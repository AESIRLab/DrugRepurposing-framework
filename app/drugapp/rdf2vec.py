# -*- coding: utf-8 -*-
"""
Module for turning graph into RDF, add drug similarity edges,
embedding the graph, and predicting new DTI

Created on Sat May 21 17:33:54 2022

@author: Carmen Reep
"""

import csv
from rdflib import Graph, Literal, URIRef
from rdflib.namespace import RDF, RDFS
import datetime
from pyrdf2vec.rdf2vec import RDF2VecTransformer
from pyrdf2vec.graphs import KG
from pyrdf2vec.walkers import RandomWalker
from pyrdf2vec.embedders import Word2Vec
import pandas as pd
import pickle
import numpy as np
import os
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
from xgboost import XGBClassifier
from scipy.stats import uniform
from sklearn.utils import class_weight
from IPython.display import HTML

today = datetime.date.today()

## RDF2VEC WITH DRUG SIMILARITIES

def is_valid_uri(value):
    """Return True if value is a non-empty, non-NaN string URI."""
    if pd.isna(value):
        return False
    if not isinstance(value, str):
        return False
    value = value.strip()
    if value == "" or value.lower() in {"na", "nan", "none"}:
        return False
    return True

def csv_to_rdf(edges_file_str, nodes_file_str, uri_type_dict):
    ''' This function turns a csv graph file into an RDF graph and saves it as a turtle file.
    Drug-gene interactions are removed to avoid bias in the drug-gene prediction.
    :param edges_file_str: the string of the edges csv of the graph to turn into rdf
    :param nodes_file_str: the string of the nodes csv of the graph to turn into rdf    
    :param uri_type_dict: a dictionary of node types in the graph to uri links of these types
    :return: RDF graph in turtle format
    '''
    # open the graph files
    graph_nodes_file = pd.read_csv(nodes_file_str)
    graph_nodes_file = graph_nodes_file.dropna(subset=["id", "uri"])
    
    # initialise the graph
    output_graph = Graph()
    
    input_file = csv.DictReader(open(edges_file_str))

    missing_nodes_count = 0
    
    for row in input_file:
        # convert row from an OrderedDict to a regular dict
        row = dict(row)
    
        subject_id = row['subject_id']
        property_uri = row['property_uri']
        object_id = row['object_id']
        property_label = row['property_label']
    
        #if link is from drug to gene, remove link in rdf graph to avoid bias while predicting
        #subject_type = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'semantic_groups'].iloc[0] 
        #object_type = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'semantic_groups'].iloc[0] 
        # Skip edges where subject or object not present in nodes table
        if subject_id not in graph_nodes_file['id'].values or object_id not in graph_nodes_file['id'].values:
            continue  # avoid out-of-bounds

        # if link is from drug to gene, remove link in rdf graph to avoid bias while predicting
        subject_type = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'semantic_groups'].iloc[0]
        object_type = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'semantic_groups'].iloc[0]

        if subject_id not in graph_nodes_file['id'].values or object_id not in graph_nodes_file['id'].values:
            missing_nodes_count += 1
            continue

        if (subject_type == 'drug') and (object_type == 'gene'):
            continue
        
        # get the uris for each node
        subject_uri = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'uri'].iloc[0]
        object_uri = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'uri'].iloc[0]
        
        # get subject and object type and turn into uri using uri_type_dict
        subject_type = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'semantic_groups'].iloc[0]
        try:
            subject_type_uri = uri_type_dict[subject_type]
        except:
            continue
        object_type = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'semantic_groups'].iloc[0]
        try:
            object_type_uri = uri_type_dict[object_type]
        except:
            continue
    
        # get subject and object label
        subject_label = graph_nodes_file.loc[graph_nodes_file['id'] == subject_id, 'preflabel'].iloc[0]
        object_label = graph_nodes_file.loc[graph_nodes_file['id'] == object_id, 'preflabel'].iloc[0]
        
     	# add types to the graph
        if is_valid_uri(subject_uri) and is_valid_uri(subject_type_uri):
            output_graph.add((URIRef(subject_uri), RDF.type, URIRef(subject_type_uri)))
        if is_valid_uri(object_uri) and is_valid_uri(object_type_uri):
            output_graph.add((URIRef(object_uri), RDF.type, URIRef(object_type_uri)))

        # add labels to the graph
        if is_valid_uri(subject_uri):
            output_graph.add((URIRef(subject_uri), RDFS.label, Literal(subject_label)))
        if is_valid_uri(object_uri):
            output_graph.add((URIRef(object_uri), RDFS.label, Literal(object_label)))

        # add properties to the graph
        if is_valid_uri(property_uri):
            output_graph.add((URIRef(property_uri), RDF.type, RDF.Property))
            output_graph.add((URIRef(property_uri), RDFS.comment, Literal(property_label)))

        # add triples to the graph
        if (
            is_valid_uri(subject_uri)
            and is_valid_uri(property_uri)
            and is_valid_uri(object_uri)
        ):
            output_graph.add((URIRef(subject_uri), URIRef(property_uri), URIRef(object_uri)))

    print(f"Skipped {missing_nodes_count} edges due to missing subject/object nodes.")
    
    # serialize and save the graph as turtle format
    output_graph.serialize(destination='my_graph_removed.ttl', format='ttl', encoding="utf-8")
    
    return output_graph.serialize(format="ttl")

def rdf_to_vec(nodes_file_str, drug_sim_str, rdf_graph, semantic_group='gene'):
    """This function turns an RDF graph into vectors, one for each entity in the specified semantic group,
    and saves it as a dictionary to a file."""
    
    print('rdf2vec running')
    # Open the graph files
    graph_nodes_file = pd.read_csv(nodes_file_str)
    drug_sim_file = pd.read_csv(drug_sim_str)
    print('files opened')
    
    # Entities we want to classify
    if semantic_group == 'drug':
        subj = drug_sim_file['subject_id'].tolist()
        obj = drug_sim_file['object_id'].tolist()
        alldrugs = list(set(subj + obj))
        entities = graph_nodes_file.loc[graph_nodes_file['id'].isin(alldrugs), 'uri'].tolist()
    else:
        entities = graph_nodes_file.loc[
            graph_nodes_file['semantic_groups'] == semantic_group, 'uri'
        ].tolist()

    print(f'entities found: {len(entities)}')
    #entities = [
    #    e for e in entities if isinstance(e, str) and e.strip() != "" and e.lower() not in {"nan", "none"}
    #]

<<<<<<< HEAD
    #if len(entities) == 0:
    #    raise ValueError(f"No valid {semantic_group} URIs found in {nodes_file_str}")

    # Initialize the transformer
    transformer = RDF2VecTransformer(
=======
    with open(os.path.join(os.getcwd(), 'extracted_entities.csv'), 'w') as f:
        for e in entities:
            f.write(f"{e},\n")
    
    if len(entities) == 0:
        raise ValueError(f"No valid {semantic_group} URIs found in {nodes_file_str}")

    # Initialize the transformer
    transformer = RDF2VecTransformer(
        Word2Vec(),
>>>>>>> 6c7d80776fdcc54cd34da8d6e6f98affd13a59f9
        walkers=[RandomWalker(max_depth=4, max_walks=10, n_jobs=2, random_state=22)],
        verbose=1,
    )

    # Build knowledge graph
    kg = KG(rdf_graph, fmt='turtle')

    print(f"[rdf2vec.py] [ref_to_vec] number of entities in KG: {len(kg._entities)}")
    with open(os.path.join(os.getcwd(), 'kg_entities.csv'), 'w') as f:
        for e in kg._entities:
            f.write(f"{e.name},\n")
    
    # 🔹 NEW STEP: filter out entities not in the RDF graph (isolated URIs)
<<<<<<< HEAD
    valid_entities = list(set(entities))
    print(f"{len(valid_entities)} entities after filtering disconnected ones (from {len(entities)})")

    if len(valid_entities) == 0:
        raise ValueError(f"No valid {semantic_group} URIs remain after filtering disconnected ones.")
    
    embeddings, literals = transformer.fit_transform(kg, valid_entities)

    # 🔹 NEW: defensive retry to catch internal PyRDF2Vec indexing issues
    #try:
    #    embeddings, literals = transformer.fit_transform(kg, valid_entities)
    #except IndexError as e:
    #    print("⚠️ PyRDF2Vec walk alignment issue detected; attempting recovery...")
    #    # Second attempt: filter entities that have at least one triple in the RDF
    #    from rdflib import Graph, URIRef
    #    rdf_g = Graph()
    #    rdf_g.parse(rdf_graph, format="turtle")

    #    connected_entities = []
    #    for e in valid_entities:
    ##        uri = URIRef(e)
    #        if (uri, None, None) in rdf_g or (None, None, uri) in rdf_g:
    #            connected_entities.append(e)

    #    print(f"Retrying with {len(connected_entities)} connected entities out of {len(valid_entities)}.")
    #    embeddings, literals = transformer.fit_transform(kg, connected_entities)
    #    valid_entities = connected_entities
=======
    # valid_entities = entities
    # print(f"{len(valid_entities)} entities after filtering disconnected ones (from {len(entities)})")

    # if len(valid_entities) == 0:
    #     raise ValueError(f"No valid {semantic_group} URIs remain after filtering disconnected ones.")

    # 🔹 NEW: defensive retry to catch internal PyRDF2Vec indexing issues
    # try:
    valid_entities = list(set(entities))
    print(f"[rdf2vec.py] [ref_to_vec] unique entities: {len(valid_entities)}")
    embeddings, literals = transformer.fit_transform(kg, valid_entities)
    # except IndexError as e:
    #     print(e)
    #     print("⚠️ PyRDF2Vec walk alignment issue detected; attempting recovery...")
    #     # Second attempt: filter entities that have at least one triple in the RDF
    #     from rdflib import Graph, URIRef
    #     rdf_g = Graph()
    #     rdf_g.parse(rdf_graph, format="turtle")

    #     connected_entities = []
    #     for e in valid_entities:
    #         uri = URIRef(e)
    #         if (uri, None, None) in rdf_g or (None, None, uri) in rdf_g:
    #             connected_entities.append(e)

    #     print(f"Retrying with {len(connected_entities)} connected entities out of {len(valid_entities)}.")
    #     embeddings, literals = transformer.fit_transform(kg, connected_entities)
    #     valid_entities = connected_entities
>>>>>>> 6c7d80776fdcc54cd34da8d6e6f98affd13a59f9

    #disconnected = len(valid_entities) - len(connected_entities)
    #print(f"Removed {disconnected} disconnected entities (no triples in RDF).")

    print('embeddings created')

    # Save embeddings as dictionary
    entity_embedding_dict = dict(zip(valid_entities, embeddings))
    with open(f"{semantic_group}_embedding_dict.pkl", "wb") as f:
        pickle.dump(entity_embedding_dict, f)

    return entity_embedding_dict


def fuse_embeddings(gene_embedding_dict, drug_embedding_dict, drug_edges_file_str, graph_nodes_file_str, genes_of_interest):
    ''' This function fuses the embeddings for the genes and drugs that have a known link,
    and adds the same amount of drug-gene pairs of random non-known links
    :param gene_embedding_dict: all gene embeddings (dict: key is gene, value is embedding)
    :param drug_embedding_dict: all drug embeddings (dict: key is drug, value is embedding)
    :param drug_edges_file_str: the string of the edges csv of the drug graph
    :param graph_nodes_file_str: the string of the nodes csv of the whole graph
    :param genes_of_interest: the genes for which predictions are wanted
    :return: training df of drug-gene pair with fused embeddings and class label, and prediction df
    '''
    array_drugs = np.array(list(drug_embedding_dict.keys()))
    drugs_emb = np.array(list(drug_embedding_dict.values()))
    array_genes = np.array(list(gene_embedding_dict.keys()))
    
    dgidb_edges_df = pd.read_csv(drug_edges_file_str)
    graph_nodes_df = pd.read_csv(graph_nodes_file_str)

    # Create a gene-drug interaction dataframe
    gene_drug_df_total = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
    
    array_drugs_len = len(array_drugs)
    for gene in array_genes:
        gene_emb = gene_embedding_dict[gene]
        gene_iri = gene
        gene_id = graph_nodes_df.loc[graph_nodes_df['uri'] == gene_iri, 'id'].iloc[0]
        fused_emb = np.multiply(gene_emb, drugs_emb)  # Hadamard operation
        gene_drug_df = pd.DataFrame({'fused_embedding': fused_emb.tolist()})
        gene_drug_df.insert(loc=0, column='drug', value=array_drugs)
        gene_drug_df.insert(loc=0, column='gene', value=[gene_iri] * array_drugs_len)
        
        # Add class column
        if gene_id not in dgidb_edges_df['object_id'].unique():
            class_lst = [0] * array_drugs_len  # No interaction
        else:
            class_lst = []
            for drug_iri in array_drugs:
                drug_id = graph_nodes_df.loc[graph_nodes_df['uri'] == drug_iri, 'id'].iloc[0]
                if len(dgidb_edges_df.loc[(dgidb_edges_df['subject_id'] == drug_id) & (dgidb_edges_df['object_id'] == gene_id), 'property_id']) != 0:
                    class_lst.append(dgidb_edges_df.loc[(dgidb_edges_df['subject_id'] == drug_id) & (dgidb_edges_df['object_id'] == gene_id), 'property_id'].item())
                else:
                    class_lst.append(0)  # No interaction
        gene_drug_df['class'] = class_lst
        
        gene_drug_df_total = pd.concat([gene_drug_df_total, gene_drug_df], ignore_index=True)
    
    print(f'[rdf2vec.py][fuse_embeddings] gene_drug_df_total head is:{gene_drug_df_total.head()}')
    # Split on class (known and unknown)
    gene_drug_df_known = gene_drug_df_total[gene_drug_df_total['class'] != 0]
    gene_drug_df_unknown = gene_drug_df_total[gene_drug_df_total['class'] == 0]
    
    # Split unknown into predict and train
    gene_drug_df_unknown_to_predict = gene_drug_df_unknown[gene_drug_df_unknown['gene'].isin(genes_of_interest)]
    gene_drug_df_unknown_not_predict = gene_drug_df_unknown[~gene_drug_df_unknown['gene'].isin(genes_of_interest)]

    if gene_drug_df_known.empty and gene_drug_df_unknown_not_predict.empty:
        print("No known or unknown DTIs; returning empty predict_df.")
        return pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class']), pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
    
    k = min(len(gene_drug_df_unknown_not_predict), 1000) if not gene_drug_df_known.empty else min(len(gene_drug_df_unknown_not_predict), 100)
    print(f'(rdf2vec.py)(fuse_embeddings) value of k sample of data to train of unknown_not_predict df: {k}')
    
    if len(gene_drug_df_unknown_not_predict) > 0:
        gene_drug_df_unknown_to_train = gene_drug_df_unknown_not_predict.sample(n=min(k, len(gene_drug_df_unknown_not_predict)))
        gene_drug_df_unknown_to_predict = gene_drug_df_unknown_not_predict.drop(gene_drug_df_unknown_to_train.index)
        predict_df = gene_drug_df_unknown_to_predict[['gene', 'drug', 'fused_embedding', 'class']].copy()
    else:
        predict_df = pd.DataFrame(columns=['gene', 'drug', 'fused_embedding', 'class'])
    
    # Combine known and unknown for training
    train_df = pd.concat([gene_drug_df_known, gene_drug_df_unknown_to_train], ignore_index=True)
    
    # Save output files
    path = os.getcwd() + '/embedding'
    if not os.path.isdir(path): os.makedirs(path)
    train_df.to_csv(f'{path}/genetarget_train_v{today}.csv', index=False)
    if not predict_df.empty:
        predict_df.to_csv(f'{path}/genetarget_predict_v{today}.csv', index=False)
    
    return train_df, predict_df

def get_genes_of_interest(edges_df, nodes_df, phenotype_of_interest):
    ''' This function finds all genes that cause the phenotype of interest.
    :param edges_df: the dataframe of the monarch edges
    :param nodes_df: the dataframe of the monarch nodes
    :param phenotype_of_interest: the id (string) of the phenotype of interest (e.g. 'HP:0002072' for chorea)
    :return: a list of iris of all genes of interest
    '''
    genes = set()
    # Original: Gene -> disease (causal)
    causal_edges = edges_df[edges_df['property_id'] == 'RO:0003303']
    genes.update(causal_edges[causal_edges['object_id'] == phenotype_of_interest]['subject_id'])
    
    # Enhanced: Pheno -> disease -> genes, plus direct gene-phenotype links
    if phenotype_of_interest.startswith('HP:'):
        # Pheno -> disease
        pheno_disease_edges = edges_df[edges_df['subject_id'] == phenotype_of_interest]
        diseases = pheno_disease_edges['object_id'].unique()
        for disease in diseases:
            disease_genes = edges_df[(edges_df['object_id'] == disease) & (edges_df['property_id'] == 'RO:0003303')]['subject_id']
            genes.update(disease_genes)
        
        # Direct gene -> phenotype (e.g., HP:0000726)
        direct_edges = edges_df[(edges_df['object_id'] == phenotype_of_interest) & (edges_df['property_id'] == 'RO:0002200')]
        genes.update(direct_edges['subject_id'])
    
    print(f'Genes of interest: {len(genes)} from symptom: {phenotype_of_interest}')
    return list(genes)          

def ml_prediction(train_df, predict_df, date):
    ''' This function builds a ML model using the training data and predicts the interactions of interest
    :param train_df: the gene, drug, embedding, class dataframe used for training
    :param predict_df: gene, drug, embedding, class dataframe used for prediction
    :param date: the date of the creation of the graph file
    :return: the predictions for the interactions in predict_df
    '''
    print('ML model is being trained, be patient')
    
    # X is the embeddings, turn list of embeddings into columns
    emb_col = train_df['fused_embedding']
    X = pd.DataFrame.from_dict(dict(zip(emb_col.index, emb_col.values))).T
    X = X.astype(float)
    
    # y is the labels, map to integers (0 for unknown, 1+ for known)
    y = train_df['class'].copy()
    class_mapping = {str(val): idx + 1 for idx, val in enumerate(train_df['class'].unique()) if str(val) != '0'}
    class_mapping['0'] = 0  # Unknown class is 0
    y = y.map(class_mapping).astype(int)

    unique_classes = np.unique(y)
    if len(unique_classes) == 1:
        print("Only one class detected; adding dummy positives for training.")
        mask = np.random.random(len(y)) < 0.01
        y = y.copy()
        y[mask] = 1  # Use integer 1 for dummy positives
        objective = 'binary:logistic'
    else:
        objective = 'multi:softmax' if len(unique_classes) > 2 else 'binary:logistic'
    
    # Define parameters to be tuned
    parameters = {
        'min_child_weight': [2, 3, 5, 8, 13, 20, 30],
        'gamma': [0, 0.2, 0.5, 0.8, 1.2, 1.6, 2.5, 4, 6],
        'reg_alpha': [0, 0.5, 1, 3, 5, 10],
        'reg_lambda': [0, 0.5, 1, 3, 5, 10],
        'subsample': uniform(0.5, 0.5),
        'colsample_bytree': uniform(0.2, 0.8),
        'max_depth': [4, 6, 8, 10, 12, 14, 16],
        'n_estimators': [35, 45, 50, 70, 80, 90, 100],
        'learning_rate': uniform(0, 0.3),
    }
    
    # Define the model
    xgb_model_hyp = XGBClassifier(objective=objective, use_label_encoder=False)

    # Find the best hyperparameters using repeated stratified k-fold and randomized search
    rskf = RepeatedStratifiedKFold(n_splits=10, n_repeats=5, random_state=123)
    randomized_search = RandomizedSearchCV(xgb_model_hyp, param_distributions=parameters,
                                          scoring='f1_weighted', n_iter=20, n_jobs=-1,
                                          error_score='raise', cv=rskf.split(X, y), verbose=3,
                                          refit=True)
    
    # Compute sample weights
    weight = class_weight.compute_sample_weight('balanced', y)
    randomized_search.fit(X, y, sample_weight=weight)

    # Find best params and score
    print('best found hyperparameters:', randomized_search.best_params_)
    print('score of the model:', randomized_search.best_score_)
    
    # Predict
    if predict_df.empty:
        print("No predictions; returning empty DF.")
        return pd.DataFrame(columns=['drug', 'gene', 'predicted_interaction', 'prob'])
    
    emb_col_pred = predict_df['fused_embedding']
    X_pred = pd.DataFrame.from_dict(dict(zip(emb_col_pred.index, emb_col_pred.values))).T
    X_pred = X_pred.astype(float)

    predictions = randomized_search.predict(X_pred)
    predictions_prob = randomized_search.predict_proba(X_pred)
    
    # Map predictions back to original class labels for readability
    reverse_mapping = {v: k for k, v in class_mapping.items()}
    predictions = [reverse_mapping.get(pred, '0') for pred in predictions]
    
    interaction_predictions_df = pd.DataFrame({
        'drug': predict_df['drug'],
        'gene': predict_df['gene'],
        'predicted_interaction': predictions,
        'prob': np.max(predictions_prob, axis=1)
    })
    
    # Add labels to drugs and genes for readability
    graph_fname_nodes = f'./graph/graph_nodes_v{date}.csv'
    nodes_df = pd.read_csv(graph_fname_nodes)
    interaction_predictions_df['drug'] = interaction_predictions_df['drug'].apply(
        lambda drug_uri: f'<a href="{drug_uri}">{nodes_df.loc[nodes_df["uri"] == drug_uri, "preflabel"].iloc[0]}</a>'
        if not nodes_df.loc[nodes_df["uri"] == drug_uri].empty else drug_uri
    )
    interaction_predictions_df['gene'] = interaction_predictions_df['gene'].apply(
        lambda gene_uri: f'<a href="{gene_uri}">{nodes_df.loc[nodes_df["uri"] == gene_uri, "preflabel"].iloc[0]}</a>'
        if not nodes_df.loc[nodes_df["uri"] == gene_uri].empty else gene_uri
    )
    
    # Save output files
    path = os.getcwd() + '/predictions'
    if not os.path.isdir(path): os.makedirs(path)
    interaction_predictions_df.to_csv(f'{path}/drug_gene_pred_v{today}.csv', index=False)
    
    return interaction_predictions_df

def drug_list_rank(symptom_user_input, disease_name, interaction_predictions_df, min_prob=0.5):
    ''' This function ranks the found drugs according to the number of genes
    this drug interacts with and sorts according to prediction confidence.
    :param interaction_predictions_df: the predictions for the interactions of interest
    :return: html file of ranked list of drugs with interactions scores
    '''

    # drop 0 predictions
    df_predict_nozero=interaction_predictions_df[interaction_predictions_df['predicted_interaction'] != '0']
    # drop interactions with min_prob of belonging in that class
    df_predict_highprob=df_predict_nozero[df_predict_nozero['prob'] > min_prob]
    df_predict_highprob['geneinteraction']=list(zip(df_predict_highprob.gene, df_predict_highprob.predicted_interaction, df_predict_highprob.prob))
    groupedby = df_predict_highprob.groupby('drug')['geneinteraction'].apply(lambda x: x.values.tolist()).to_frame()
    # add how many interactions each drug has
    groupedby['nr_interactions'] = groupedby.apply(lambda row: len(row.geneinteraction), axis=1)
    # add sum of score
    groupedby['scoresum'] = groupedby.apply(lambda row: sum(n for _,_,n in row.geneinteraction), axis=1)
    
    # sort in descending order
    groupedby_sorted = groupedby.sort_values(by=['nr_interactions', 'scoresum'], ascending=False)
    groupedby_sorted = groupedby_sorted.reset_index()
    groupedby_sorted = groupedby_sorted[['drug', 'geneinteraction']]

    # unstack geneinteraction column
    groupedby_sorted_unstacked= groupedby_sorted.set_index('drug', append=True).geneinteraction.apply(pd.Series).stack().reset_index(level=[0, 2], drop=True).reset_index()
    # create columns of tuples
    groupedby_sorted_unstacked[['gene', 'interaction', 'score']] = pd.DataFrame(groupedby_sorted_unstacked[0].tolist(), index=groupedby_sorted_unstacked.index)
    # drop tuple columns
    groupedby_sorted_unstacked = groupedby_sorted_unstacked.drop(0, 1)

    # create list of tuples for the multi-index
    groupedby_sorted_unstacked['druggene'] = list(zip(groupedby_sorted_unstacked.drug, groupedby_sorted_unstacked.gene))
    tuples_druggene = groupedby_sorted_unstacked['druggene'].tolist()
    index = pd.MultiIndex.from_tuples(tuples_druggene, names=('Drug', 'Gene'))
    # create the multi-index dataframe
    df_multiindex = pd.DataFrame(groupedby_sorted_unstacked[['interaction', 'score']].to_numpy(), index=index,
                                 columns=['Interaction type', 'Confidence'])
    
    # turn interaction type ID back to string for easier read
    df_multiindex.loc[df_multiindex['Interaction type'] == 'RO:0002408', 'Interaction type'] = 'inhibits'
    df_multiindex.loc[df_multiindex['Interaction type'] == 'RO:0002406', 'Interaction type'] = 'activates'
    df_multiindex.loc[df_multiindex['Interaction type'] == 'RO:0011002', 'Interaction type'] = 'regulates'
    
    # round predictions to 3 decimal places
    df_multiindex['Confidence'] = df_multiindex['Confidence'].apply(lambda x: round(x, 3))
    
    #change directory back to drugapp
    os.chdir('..')
    os.chdir('..')
    os.chdir('..')

    pd.set_option('colheader_justify', 'left')   # FOR TABLE <th>
    
    html_string = '''
    {{% extends "layout.html" %}}
    {{% block content %}}
    <html>
      <h2>Predictions </h2>
      Disease:<h4> {diseasename} </h4>
      Symptom: <h4>{symptom} </h4>
      Disease graph created on {diseasedate}. <br>
      Click on the drug or gene for information about that drug or gene.
      
      <body>
        {table}
      </body>
    </html>.
    {{% endblock content %}}
    '''
    with open('./drugapp/templates/{}; {}.html'.format(disease_name, symptom_user_input), 'w') as f:
        f.write(html_string.format(diseasename = disease_name[:-13], diseasedate = disease_name[-11:-1], symptom = symptom_user_input, table=df_multiindex.to_html(escape=False,classes='mystyle')))
    
    return df_multiindex

def rdf2vec_general(symptom_user_input, date, disease_name_date):
    """
    This function runs the rdf2vec script with the specified symptom by the user.
    :param symptom_user_input: the input specified symptom
    :param date: the date of creation of the disease graph
    :return: html file of the resulting predictions table
    """
    today = datetime.date.today()

    # Base path for this disease
<<<<<<< HEAD
    #base_path = f'./drugapp/data/{disease_name_date}/monarch'
    base_path = os.path.join(os.getcwd(), 'monarch')
    
=======
    print(f'[rdf2vec.py] [rdf2vec_general] os.getcwd(): ${os.getcwd()}')
    base_path = os.path.join(os.getcwd(), 'monarch')# f'./drugapp/data/{disease_name_date}/monarch'
>>>>>>> 6c7d80776fdcc54cd34da8d6e6f98affd13a59f9
    
    monarch_edges_dis_file = '{}/monarch_edges_disease_v{}.csv'.format(base_path, date) 
    monarch_edges_symp_file = '{}/monarch_edges_symptom_v{}.csv'.format(base_path, today)    
    monarch_nodes_dis_file = '{}/monarch_nodes_disease_v{}.csv'.format(base_path, date)
    monarch_nodes_symp_file = '{}/monarch_nodes_symptom_v{}.csv'.format(base_path, today)
    # open csv files
    edges_dis_df = pd.read_csv(monarch_edges_dis_file)
    edges_symp_df = pd.read_csv(monarch_edges_symp_file)    
    nodes_dis_df = pd.read_csv(monarch_nodes_dis_file)
    nodes_symp_df = pd.read_csv(monarch_nodes_symp_file)
    # concatenate dfs
    combined_monarch_nodes_df = pd.concat([nodes_dis_df, nodes_symp_df])
    combined_monarch_edges_df = pd.concat([edges_dis_df, edges_symp_df])
    
    symptom_id = combined_monarch_nodes_df[['id']].loc[combined_monarch_nodes_df['preflabel'] == symptom_user_input].values[0][0]

    uriType_dict = {
        'anatomical entity': 'http://identifiers.org/sio:SIO_001262', 
        'disease': 'http://identifiers.org/sio:SIO_010299', 
        'molecular function': 'https://identifiers.org/GO:0003674', 
        'biological process': 'https://identifiers.org/GO:0008150', 
        'genotype': 'http://identifiers.org/sio:SIO_001079', 
        'gene': 'http://identifiers.org/sio:SIO_010035', 
        'phenotype': 'http://identifiers.org/sio:SIO_010056', 
        'variant': 'http://identifiers.org/sio:SIO_001381', 
        'pathway': 'http://identifiers.org/sio:SIO_001107',
        'cellular component': 'https://identifiers.org/GO:0005575', 
        'drug': 'http://identifiers.org/sio:SIO_010038', 
        }

    # create the rdf graph (where drug-gene interactions are not included)
    graph_edges_str = './graph/graph_edges_v{}.csv'.format(today)
    graph_nodes_str = './graph/graph_nodes_v{}.csv'.format(today)
    # create the RDF graph save as .ttl file
    csv_to_rdf(graph_edges_str, graph_nodes_str, uriType_dict)
    
    # embed the RDF graph focusing on genes and drugs
    drugsim_file_str = './similaritygraph/drugdrugsim_v{}.csv'.format(today)
    gene_embedding_dict = rdf_to_vec(graph_nodes_str, drugsim_file_str, 'my_graph_removed.ttl', semantic_group = 'gene')
    drug_embedding_dict = rdf_to_vec(graph_nodes_str, drugsim_file_str, 'my_graph_removed.ttl', semantic_group = 'drug')
    
    # find genes of interest (these will not be in training set)    
    genes_of_interest = get_genes_of_interest(combined_monarch_edges_df, combined_monarch_nodes_df, symptom_id)
    
    # fuse embeddings of drug-gene pairs and create training and prediction df
    
    drug_edges_file_str = './DGIdb/DGIdb_edges_v{}.csv'.format(today)
    graph_nodes_file_str = './graph/graph_nodes_v{}.csv'.format(today)
    train_df, predict_df = fuse_embeddings(gene_embedding_dict, drug_embedding_dict, drug_edges_file_str, graph_nodes_file_str, genes_of_interest)

    # predict the unkown interactions of interest
    interaction_predictions_df = ml_prediction(train_df, predict_df, today)
    # rank the found drugs

    interaction_predictions_df = pd.read_csv('./predictions/drug_gene_pred_v{}.csv'.format(today))
    drug_list_rank(symptom_user_input, disease_name_date, interaction_predictions_df,min_prob=0.9)
    
    # delete the files not needed anymore (everything except the original disease graph)
    # these extra files take up a lot of space
    os.chdir('./drugapp/data/{}'.format(disease_name_date)) # go into the disease folder
    os.remove('./monarch/monarch_edges_symptom_v{}.csv'.format(today))
    os.remove('./monarch/monarch_nodes_symptom_v{}.csv'.format(today))
    os.remove('./monarch/monarch_orthopeno_network_symptom_v{}.csv'.format(today))
    os.remove('./DGIdb/DGIdb_edges_v{}.csv'.format(today))
    os.remove('./DGIdb/DGIdb_nodes_v{}.csv'.format(today))
    os.remove('./DGIdb/DGIdb_network_v{}.csv'.format(today))
    os.remove('./similaritygraph/drugdrugsim_v{}.csv'.format(today))
    os.remove('./graph/graph_edges_v{}.csv'.format(today)) # keep for neo4j
    os.remove('./graph/graph_nodes_v{}.csv'.format(today)) # keep for neo4j
    os.remove('drug_embedding_dict.pkl')
    os.remove('gene_embedding_dict.pkl')
    os.remove('./embedding/genetarget_predict_v{}.csv'.format(today))
    os.remove('./embedding/genetarget_train_v{}.csv'.format(today))
    os.remove('my_graph_removed.ttl')
    
  
    
if __name__ == '__main__':
    
    #If a previous graph is selected, make sure current directory is that folder
    # and monarch and dgidb and graph files are the latest in the folder
    if _userinput_:
        diseasename = 'Huntington disease (2022-06-27)' # the user input
        symptoms_name_lst, date, input_number= symptom_list_specified_folder(input_folder=diseasename)
        ## ASK USER INPUT ##
        symptom_user_input= 'Anxiety' # the user input
        #run_monarch_symptom(input_symptom, date)
        #run_dgidb(date)
        #run_drugsimilarity()
        #run_combine_graphs(date)
        rdf2vec_general(symptom_user_input, date, diseasename)

    else: # a new graph is created with phenotype MIM number
        symptoms_name_lst, date, diseasename = symptom_list_today()
        ## ASK USER INPUT ##
        symptom_user_input= 'Anxiety' # the user input
        #run_monarch_symptom(input_symptom, date)
        #run_dgidb(date)
        #run_drugsimilarity()
        #run_combine_graphs(date)
        rdf2vec_general(symptom_user_input, date, diseasename)
    

# -*- coding: utf-8 -*-
"""
Module for DGIdb network preparation and management

Created on Wed Apr  6 10:49:08 2022

@author: Carmen Reep
"""

import requests
import datetime
import pandas as pd
from biothings_client import get_client
import ast
import os

# timestamp
today = datetime.date.today()


def get_genes(nodes_df):
    """
    This functions finds all genes in Monarch. It returns a dict of id, name pairs of the gene IDs  
    :param monarch_nodes_csv: csv file of monarch nodes
    :return: dict of gene ids/name
    """
    print('get_genes() is running')
    
    # keep rows with 'gene' as semantic label
    df_genes = nodes_df.loc[nodes_df['semantic_groups'] == 'gene']    
    #get names and ids of genes
    gene_name_list = df_genes['preflabel'].to_list()
    #dict of id and names
    gene_id_list = df_genes['id'].to_list()
    id_name_dict = dict(zip(gene_id_list, gene_name_list))
    
    return id_name_dict

def normalize_genes_to_graph(id_name_dict):
    """
    Normalizes gene IDs to Entrez IDs for DGIdb input.
    :param id_name_dict: dict of gene IDs/names
    :return: dict mapping original IDs to Entrez IDs
    """
    print('\nMapping Monarch gene ontologies to Entrez genes...')
    # final dict
    concept_dct = dict()
    symbol_dct = dict()

    # input
    zfin = list()
    ensembl = list()
    hgnc = list()
    mgi = list()
    wormbase = list()
    flybase = list()
    
    for idx in id_name_dict:
        if ':' in idx:
            prefix = idx.split(':')[0].lower()
            id_val = idx.split(':')[1]
            if prefix == 'flybase':
                flybase.append(id_val)
            elif prefix == 'wormbase':
                wormbase.append(id_val)
            elif prefix == 'mgi':
                mgi.append(idx)  # Full MGI:ID
            elif prefix == 'hgnc':
                hgnc.append(id_val)
            elif prefix == 'ensembl':
                ensembl.append(id_val)
            elif prefix == 'zfin':
                zfin.append(id_val)
    
    mg = get_client('gene')
    # Queries with lowercase scopes, full IDs
    if mgi:
        print(f"Querying {len(mgi)} MGI IDs...")
        df_mgi = mg.querymany(mgi, scopes='mgi', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_mgi.columns:
            df_mgi['entrezgene'] = df_mgi['entrezgene'].astype('Int64')
            concept_dct.update(dict(zip(df_mgi.index, df_mgi['entrezgene'])))
        else:
            print('no MGI IDs mapped to Entrez IDs')
        if 'symbol' in df_mgi.columns:
            symbol_dct.update(dict(zip(df_mgi.index, df_mgi['symbol'])))
    
    if flybase:
        print(f"Querying {len(flybase)} FlyBase IDs...")
        df_flybase = mg.querymany(flybase, scopes='flybase', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_flybase.columns:
            df_flybase['entrezgene'] = df_flybase['entrezgene'].astype('Int64')
            concept_dct.update({f'FlyBase:{k}': v for k, v in dict(zip(df_flybase.index, df_flybase['entrezgene'])).items()})
        else:
            print('no flybase IDs mapped to Entrez IDs')
        if 'symbol' in df_flybase.columns:
            symbol_dct.update({f'FlyBase:{k}': v for k, v in dict(zip(df_flybase.index, df_flybase['symbol'])).items()})
    
    if hgnc:
        print(f"Querying {len(hgnc)} HGNC IDs...")
        df_hgnc = mg.querymany(hgnc, scopes='hgnc', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_hgnc.columns:
            df_hgnc['entrezgene'] = df_hgnc['entrezgene'].astype('Int64')
            concept_dct.update({f'HGNC:{k}': v for k, v in dict(zip(df_hgnc.index, df_hgnc['entrezgene'])).items()})
        else:
            print('no HGNC IDs mapped to Entrez IDs')
        if 'symbol' in df_hgnc.columns:
            symbol_dct.update({f'HGNC:{k}': v for k, v in dict(zip(df_hgnc.index, df_hgnc['symbol'])).items()})
    
    if wormbase:
        print(f"Querying {len(wormbase)} WormBase IDs...")
        df_wormbase = mg.querymany(wormbase, scopes='wormbase', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_wormbase.columns:
            df_wormbase['entrezgene'] = df_wormbase['entrezgene'].astype('Int64')
            concept_dct.update({f'WormBase:{k}': v for k, v in dict(zip(df_wormbase.index, df_wormbase['entrezgene'])).items()})
        else:
            print('no WormBase IDs mapped to Entrez IDs')
        if 'symbol' in df_wormbase.columns:
            symbol_dct.update({f'WormBase:{k}': v for k, v in dict(zip(df_wormbase.index, df_wormbase['symbol'])).items()})
    
    if zfin:
        print(f"Querying {len(zfin)} ZFIN IDs...")
        df_zfin = mg.querymany(zfin, scopes='zfin', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_zfin.columns:
            df_zfin['entrezgene'] = df_zfin['entrezgene'].astype('Int64')
            concept_dct.update({f'ZFIN:{k}': v for k, v in dict(zip(df_zfin.index, df_zfin['entrezgene'])).items()})
        else:
            print('no ZFIN IDs mapped to Entrez IDs')
        if 'symbol' in df_zfin.columns:
            symbol_dct.update({f'ZFIN:{k}': v for k, v in dict(zip(df_zfin.index, df_zfin['symbol'])).items()})
    
    if ensembl:
        print(f"Querying {len(ensembl)} Ensembl IDs...")
        df_ensembl = mg.querymany(ensembl, scopes='ensembl.gene', fields='entrezgene,symbol', size=1, as_dataframe=True)
        if 'entrezgene' in df_ensembl.columns:
            df_ensembl['entrezgene'] = df_ensembl['entrezgene'].astype('Int64')
            concept_dct.update({f'ENSEMBL:{k}': v for k, v in dict(zip(df_ensembl.index, df_ensembl['entrezgene'])).items()})
        else:
            print('no ensembl IDs mapped to Entrez IDs')
        if 'symbol' in df_ensembl.columns:
            symbol_dct.update({f'ENSEMBL:{k}': v for k, v in dict(zip(df_ensembl.index, df_ensembl['symbol'])).items()})
    
    # Filter non-NA
    concept_dct = {k: v for k, v in concept_dct.items() if pd.notna(v)}
    symbol_dct = {k: v for k, v in symbol_dct.items() if pd.notna(v)}
    print(f"Mapped {len(concept_dct)} genes to Entrez IDs")
    print(f"Mapped {len(symbol_dct)} genes to symbols")
    return concept_dct, symbol_dct

def hit_dgidb_api(genes_list, rows=2000):
    """
    Performs GraphQL API calls to DGIdb v5 for drug-gene interactions.
    :param genes_list: list of gene symbols (strings)
    :param rows: max results
    :return: API response dict
    """
    if not genes_list:
        print("No genes to query; skipping DGIdb API.")
        return {'data': {'genes': {'nodes': []}}}
    
    print('\nThe function "hit_dgidb_api()" is running...This may take a while')
    endpoint = 'https://dgidb.org/api/graphql'
    query = """
    query Interactions($genes: [String!]) {
      genes(names: $genes) {
        nodes {
          name
          conceptId
          interactions {
            drug { name, conceptId }
            interactionScore
            interactionTypes { type }
            sources { sourceDbName }
            publications { pmid }
          }
        }
      }
    }
    """
    # genes_list is symbols
    print(f"Querying DGIdb with {len(genes_list)} gene names: {genes_list[:5]}...")  # Debug print
    variables = {'genes': genes_list}
    response = requests.post(endpoint, json={'query': query, 'variables': variables})
    
    if response.status_code != 200:
        print(f"DGIdb API error: {response.status_code}, text: {response.text}")
        return {'data': {'genes': {'nodes': []}}}
    
    json_data = response.json()
    print(f"DGIdb returned {len(json_data.get('data', {}).get('genes', {}).get('nodes', []))} gene nodes")  # Debug
    print('\nFinished "hit_dgidb_api()".')
    return json_data

def get_edges_objects(data):
    """
    Prepares DGIdb GraphQL response.
    :param data: DGIdb API response dict
    :return: gene_drug_dict {gene: [interactions]}
    """
    drug_dict = dict()
    gene_nodes = data.get('data', {}).get('genes', {}).get('nodes', [])
    
    for gene in gene_nodes:
        gene_id = gene.get('conceptId', 'NA')
        gene_name = gene.get('name', 'NA')
        interactions = gene.get('interactions', [])
        
        for interaction in interactions:
            if gene_id not in drug_dict:
                drug_dict[gene_id] = []
            drug_dict[gene_id].append({
                'drugName': interaction.get('drug', {}).get('name', 'NA'),
                'drugConceptId': interaction.get('drug', {}).get('conceptId', 'NA'),
                'interactionTypes': [t.get('type', 'NA') for t in interaction.get('interactionTypes', [])],
                'pmids': [p.get('pmid', 'NA') for p in interaction.get('publications', [])],
                'sources': [s.get('sourceDbName', 'NA') for s in interaction.get('sources', [])],
                'score': interaction.get('interactionScore', 'NA')
            })
    
    print(f"DGIdb interactions mapped for {len(drug_dict)} genes")
    return drug_dict

def create_csv(id_name_dict, drug_dict, nodes_df, filename='DGIdb_network'):
    """
    This function saves the DGIdb network as a csv file
    :param id_name_dict: dict of gene ids/name to convert
    :param drug_dict: dict of searchterm/drugs 
    :param nodes_df: the dataframe of all monarch nodes 
    :param filename: the name of the csv file to save
    :return: location of the saved csv file
    """   
    
    df = pd.DataFrame(columns = ['target_id', 'target_name', 'target_iri', 'target_category', 'interaction_id', 'interaction_types', 'interaction_iri', 'drug_id', 
                                    'drug_name', 'drug_iri', 'drug_category', 'pm_ids', 'sources', 'score'] )
    
    ## create dict to turn activationtype into ols ontology
    interaction_dict = {}
    
    #inhibits
    interaction_dict['RO:0002408'] = ['antagonist', 'antibody', 'antisense oligonucleotide', 'blocker', 'cleavage',
                                      'inhibitor', 'inhibitory allosteric modulator', 'inverse agonist', 
                                      'negative modulator', 'partial antagonist', 'suppressor'] 
    #activates
    interaction_dict['RO:0002406'] = ['activator', 'agonist', 'chaperone', 'cofactor', 'inducer', 'partial agonist',
                                      'positive modulator', 'stimulator', 'vaccine'] 
    #regulates
    interaction_dict['RO:0011002'] = ['NA', 'None', 'n/a', 'other/unknown', 'adduct', 'allosteric modulator',
                                      'binder', 'ligand', 'modulator', 'multitarget', 'potentiator',
                                      'product of', 'substrate']      
    # turn into value:key
    interaction_dict_conv = {}
    for k,vs in interaction_dict.items():
        for v in vs:
            interaction_dict_conv[v] = k
            
    for searchterm, drug_l in drug_dict.items():
        interactionId_col = []
        interactionTypes_col = []
        interactionIri_col = []
        drugId_col = []
        drugName_col = []
        drugIri_col = []
        drugCategory_col = []
        pmids_col = []
        sources_col = []
        score_col = []
        
        for i in range(len(drug_l)):
            interactionId_col.append(drug_l[i]['interactionTypes']) # first type, later change to id
            interactionTypes_col.append(drug_l[i]['interactionTypes'])
            interactionIri_col.append(drug_l[i]['interactionTypes']) # first type, later change to iri
            drugId_col.append(drug_l[i]['drugConceptId'])
            drugName_col.append(drug_l[i]['drugName'])
            drugIri_col.append(drug_l[i]['drugConceptId']) # first id, later change to iri
            drugCategory_col.append('drug')
            pmids_col.append(drug_l[i]['pmids'])
            sources_col.append(drug_l[i]['sources'])
            score_col.append(drug_l[i]['score'])
            
        
        df_one_gene = pd.DataFrame(list(zip(interactionId_col, interactionTypes_col, interactionIri_col, drugId_col, 
                                   drugName_col, drugIri_col, drugCategory_col, pmids_col, sources_col, score_col)), 
                          columns =['interaction_id', 'interaction_types', 'interaction_iri', 'drug_id', 
                                    'drug_name', 'drug_iri', 'drug_category', 'pm_ids', 'sources', 'score'])
        
        #turn interaction id into ontology
        df_one_gene['interaction_id'] = df_one_gene['interaction_id'].apply(lambda x: x[0] if len(x)>0 else 'NA') # keep first 
        df_one_gene=df_one_gene.replace({"interaction_id": interaction_dict_conv})
        
        #create interaction iri column
        df_one_gene['interaction_iri'] = df_one_gene['interaction_id'].values #copy all ids from the just created id column
        df_one_gene['interaction_iri'] = ['http://purl.obolibrary.org/obo/' + idx.replace(':', '_') for idx in df_one_gene['interaction_iri']]

        #create drug iri column 
        iricol=[]
        for idx in df_one_gene['drug_iri']:
            iri = 'https://identifiers.org/'+idx
            iricol.append(iri)
        df_one_gene['drug_iri'] = iricol
        
        target_id_col = [searchterm]*len(df_one_gene.index) # create column of search term gene id, using the id_name_dict of genes
        target_name_col = [id_name_dict.get(searchterm, 'NA')]*len(df_one_gene.index)
        target_iri_col = [nodes_df.loc[nodes_df['id'] == searchterm, 'uri'].iloc[0] if not nodes_df.loc[nodes_df['id'] == searchterm].empty else 'NA']*len(df_one_gene.index) # create column of search term iri using nodes_df
        target_category_col = ['gene']*len(df_one_gene.index)
        df_one_gene.insert(loc=0, column='target_category', value=target_category_col)
        df_one_gene.insert(loc=0, column='target_iri', value=target_iri_col)       
        df_one_gene.insert(loc=0, column='target_name', value=target_name_col)
        df_one_gene.insert(loc=0, column='target_id', value=target_id_col)
        
        df_one_gene.interaction_types = df_one_gene.interaction_types.apply(lambda y: 'NA' if len(y)==0 else y) #replace empty lists in the column interaction_types with 'None'
        
        df = pd.concat([df, df_one_gene])

    # print output file
    path = os.getcwd() + '/DGIdb'
    os.makedirs(path, exist_ok=True)
    pd.DataFrame(df).fillna('NA').to_csv('{}/{}_v{}.csv'.format(path, filename, today), index=False)

    return print("\nSaving DGIdb network at: '{}/{}_v{}.csv'...\n".format(path, filename, today))

def read_connections(filename):
    """
    This function reads DBIdb CSV file.
    :param filename: the csv file string
    :return: dataframe
    """

    path = os.getcwd() + '/DGIdb'
    csv_path = path + '/' + filename

    network_df = pd.read_csv('{}'.format(csv_path))
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))

    return network_df

def build_edges(edges_df):
    """
    Builds edges network with graph schema.
    :param edges_df: network dataframe
    :return: graph edges DataFrame
    """
    print('\nThe function "build_edges()" is running...')
    
    # Safe eval only if 'sources' exists
    if 'sources' in edges_df.columns:
        edges_df['sources'] = edges_df['sources'].apply(lambda x: ast.literal_eval(str(x)) if isinstance(x, str) else x)
    else:
        edges_df['sources'] = 'NA'
    
    uriPrefixes_dct = {'pmid': 'https://pubmed.ncbi.nlm.nih.gov/'}
    ref_text = 'This edge comes from the DGIdb 5.0.'
    ref_date = today
    
    edges_l = []
    for edge in edges_df.itertuples():
        pmid_l = str(edge.pm_ids).strip('][').split(', ')
        pmid_l = [n.strip() for n in pmid_l if n.strip()]
        ref_uri = uriPrefixes_dct['pmid'] + ','.join(pmid_l) if pmid_l and pmid_l[0] else 'NA'
        
        sub_id = 'NA' if pd.isna(edge.target_id) else edge.target_id
        rel_id = 'NA' if pd.isna(edge.interaction_id) else edge.interaction_id
        obj_id = 'NA' if pd.isna(edge.drug_id) else edge.drug_id
        rel_label = 'NA' if pd.isna(edge.interaction_types) else "|".join(ast.literal_eval(str(edge.interaction_types)) if isinstance(edge.interaction_types, str) else edge.interaction_types)
        rel_iri = 'NA' if pd.isna(edge.interaction_iri) else edge.interaction_iri
        rel_def = f'interaction score = {edge.score}' if not pd.isna(edge.score) else 'NA'
        
        edges_l.append({
            'subject_id': obj_id,  # drug
            'object_id': sub_id,  # gene
            'property_id': rel_id,
            'property_label': rel_label,
            'property_description': rel_def,
            'property_uri': rel_iri,
            'reference_uri': ref_uri,
            'reference_supporting_text': ref_text,
            'reference_date': ref_date
        })
    
    df = pd.DataFrame(edges_l)
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', 'property_label', 'property_description', 'property_uri']]
    path = os.getcwd() + '/DGIdb'
    os.makedirs(path, exist_ok=True)
    df.fillna('NA').to_csv(f'{path}/DGIdb_edges_v{today}.csv', index=False)
    
    print(f'\n* This is the size of the edges file data structure: {df.shape}')
    print(f'* These are the edges attributes: {df.columns}')
    print(f'* This is the first record:\n{df.head(1)}')
    print(f'\nThe DGIdb network edges are built and saved at: {path}/DGIdb_edges_v{today}.csv\n')
    print('\nFinished build_edges().\n')
    return df
    
def build_nodes(edges_df):
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe 
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_nodes()" is running...')
    
    # Safe eval only if 'sources' exists
    if 'sources' in edges_df.columns:
        edges_df['sources'] = edges_df['sources'].apply(lambda x: ast.literal_eval(str(x)) if isinstance(x, str) else x)
    
    # build semantic groups dictionary
    # collide concepts in a concept dict
    concept_dct = dict()

    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['target_id'] = tuple[0]
            record['target_name'] = tuple[1]
            record['target_iri'] = tuple[2]    
            record['target_category'] = tuple[3]  
            record['interaction_id'] = tuple[4]
            record['interaction_types'] = tuple[5]
            record['interaction_iri'] = tuple[6]
            record['drug_id'] = tuple[7]
            record['drug_name'] = tuple[8]
            record['drug_iri'] = tuple[9]
            record['drug_category'] = tuple[10]
            record['pm_ids'] = tuple[11]
            record['sources'] = tuple[12]
            record['score'] = tuple[13]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        sid = getattr(edge, 'target_id', 'NA')
        oid = getattr(edge, 'drug_id', 'NA')
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))

    # build concept attributes dict: id integration of sub and obj IDs in a common data structure
    concept_dct = dict()

    for edge in edges_df.itertuples():
        # id: integration of sub and obj IDs in a unique data structure
        sid = getattr(edge, 'target_id', 'NA')
        slab = getattr(edge, 'target_name', 'NA')
        sgroup = getattr(edge, 'target_category', 'NA')
        siri = getattr(edge, 'target_iri', 'NA')
        oid = getattr(edge, 'drug_id', 'NA')
        olab = getattr(edge, 'drug_name', 'NA')
        ogroup = getattr(edge, 'drug_category', 'NA')
        oiri = getattr(edge, 'drug_iri', 'NA')

        concept_dct[sid] = {'preflabel': slab,
                            'semantic_groups': sgroup,
                            'uri': siri,
                            'name': 'NA',
                            'synonyms': 'NA', 
                            'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab,
                            'semantic_groups': ogroup,
                            'uri': oiri,
                            'name': 'NA',
                            'synonyms': 'NA', 
                            'description': 'NA'}

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['uri'] = concept_dct[concept]['uri']
        node['preflabel'] = concept_dct[concept]['preflabel']
        node['name'] = concept_dct[concept]['preflabel'] 
        node['synonyms'] = 'NA'
        node['description'] = 'NA' 
        nodes_l.append(node)

    # save nodes file
    df = pd.DataFrame(nodes_l)
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/DGIdb'
    os.makedirs(path, exist_ok=True)
    df.fillna('NA').to_csv('{}/DGIdb_nodes_v{}.csv'.format(path, today), index=False)

    # print info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe DGIdb network nodes are built and saved at: {}/DGIdb_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished build_nodes().\n')

    return df

def run_dgidb(date=None, disease_folder=None):
    """
    Runs the DGIdb script and saves nodes/edges files.
    :param date: date of disease graph
    :param disease_folder: disease folder name
    :return: None (saves CSVs)
    """
    path = os.getcwd() + '/DGIdb'
    os.makedirs(path, exist_ok=True)
    
    # Adjust paths if disease_folder provided
    base_path = os.getcwd()
    if disease_folder:
        base_path = os.path.join(base_path, 'drugapp', 'data', disease_folder)
    
    monarch_nodes_dis_file = os.path.join(base_path, 'monarch', f'monarch_nodes_disease_v{date}.csv')
    monarch_nodes_symp_file = os.path.join(base_path, 'monarch', f'monarch_nodes_symptom_v{date}.csv')
    
    nodes_dis_df = pd.read_csv(monarch_nodes_dis_file)
    nodes_symp_df = pd.read_csv(monarch_nodes_symp_file) if os.path.exists(monarch_nodes_symp_file) else pd.DataFrame()
    combined_monarch_nodes_df = pd.concat([nodes_dis_df, nodes_symp_df], ignore_index=True)
    
    id_name_dict = get_genes(combined_monarch_nodes_df)
    id_to_entrez_dct, id_to_symbol_dct = normalize_genes_to_graph(id_name_dict)
    entrez_to_id_dict = {v: k for k, v in id_to_entrez_dct.items()}
    seed_symbols = list(id_to_symbol_dct.values())
    print(f"Seed gene symbols sample: {seed_symbols[:5]}")  # Debug
    
    edges = hit_dgidb_api(genes_list=seed_symbols, rows=2000)
    gene_drug_dict = get_edges_objects(edges)
    
    # Remap keys from conceptId "ncbigene:entrez" to original
    for old_key in list(gene_drug_dict.keys()):
        if 'ncbigene:' in old_key:
            entrez_str = old_key.split('ncbigene:')[1]
            try:
                entrez = float(entrez_str)
                original = entrez_to_id_dict.get(entrez, None)
                if original:
                    gene_drug_dict[original] = gene_drug_dict.pop(old_key)
            except ValueError:
                print(f"Invalid Entrez in key: {old_key}")
    
    create_csv(id_name_dict, gene_drug_dict, nodes_df=combined_monarch_nodes_df, filename='DGIdb_network')
    DGIdb_connections = read_connections(f'DGIdb_network_v{today}.csv')
    build_edges(DGIdb_connections)
    build_nodes(DGIdb_connections)


if __name__ == '__main__':
    run_dgidb()
# -*- coding: utf-8 -*-
"""
Module for graph building and management

Created on Thu Jun  9 10:50:12 2022

@author: Carmen Reep
"""
import pandas as pd
import os
import datetime
import requests
from drugapp.Monarch import curie_to_uri, is_real_uri

# VARIABLES
today = datetime.date.today()

def build_edges(dgidb, monarch_disease, monarch_symptom, drugsim, input_from_file=False):
    """
    Builds the edges graph.
    :param dgidb: path to DGIdb edges file
    :param monarch_disease: path to Monarch disease edges file
    :param monarch_symptom: path to Monarch symptom edges file
    :param drugsim: path to drug similarity edges file
    :param input_from_file: boolean to read from file
    :return: edges dataframe
    """
    print('\nThe function "build_edges()" is running...')
    print(f'dgidb file is : {dgidb}')
    edges_l = []
    
    if input_from_file:
        print('\nPreparing networks...')
        # Read DGIdb edges
        if dgidb and os.path.exists(dgidb):
            dgidb_df = pd.read_csv(dgidb)
            edges_l.append(dgidb_df)
        else:
            print(f"Warning: DGIdb edges file not found at {dgidb}")
        
        # Read Monarch disease edges
        if monarch_disease and os.path.exists(monarch_disease):
            monarch_dis_df = pd.read_csv(monarch_disease)
            edges_l.append(monarch_dis_df)
        else:
            print(f"Warning: Monarch disease edges file not found at {monarch_disease}")
        
        # Read Monarch symptom edges
        if monarch_symptom and os.path.exists(monarch_symptom):
            monarch_symp_df = pd.read_csv(monarch_symptom)
            edges_l.append(monarch_symp_df)
        else:
            print(f"Warning: Monarch symptom edges file not found at {monarch_symptom}")
        
        # Read drug similarity edges
        if drugsim and os.path.exists(drugsim):
            drugsim_df = pd.read_csv(drugsim)
            edges_l.append(drugsim_df)
        else:
            print(f"Warning: Drug similarity edges file not found at {drugsim}")
    else:
        monarch_dis_df = pd.DataFrame(monarch_disease)
        monarch_symp_df = pd.DataFrame(monarch_symptom)
        dgidb_df = pd.DataFrame(dgidb)
        drugsim_df = pd.DataFrame(drugsim)
        edges_l = [monarch_dis_df, monarch_symp_df, dgidb_df, drugsim_df]
    
    print('Monarch:')
    print(f'Disease edges: {monarch_dis_df.shape}, Columns: {monarch_dis_df.columns}')
    print(f'Symptom edges: {monarch_symp_df.shape}, Columns: {monarch_symp_df.columns}')
    print('DGIdb:')
    print(f'Edges: {dgidb_df.shape}, Columns: {dgidb_df.columns}')
    print('Drug similarity edges:')
    print(f'Edges: {drugsim_df.shape}, Columns: {drugsim_df.columns}')
    
    # Concatenate edges
    print('\nConcatenating into a graph...')
    if edges_l:
        statements = pd.concat(edges_l, ignore_index=True, join="inner")
    else:
        statements = pd.DataFrame(columns=['subject_id', 'property_id', 'object_id', 'reference_uri',
                                           'reference_supporting_text', 'reference_date', 'property_label',
                                           'property_description', 'property_uri'])
    
    print(f'Combined shape: {statements.shape}')
    
    # Drop duplicates
    print('\nDrop duplicated rows...')
    statements.drop_duplicates(keep='first', inplace=True)
    print(f'After dropping duplicates: {statements.shape}')
    
    # Add property_uri for CURIEs
    curie_dct = {
        'ro': 'http://purl.obolibrary.org/obo/',
        'bfo': 'http://purl.obolibrary.org/obo/',
        'geno': 'http://purl.obolibrary.org/obo/',
        'dc': 'http://purl.org/dc/elements/1.1/',
        'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
        'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
        'skos': 'http://www.w3.org/2004/02/skos/core#',
        'pato': 'http://purl.obolibrary.org/obo/',
        'sio': 'http://semanticscience.org/resource/',
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',
        'encode': 'https://www.encodeproject.org/search/?searchTerm='
    }
    
    for i, row in statements.iterrows():
        if ':' in str(row['property_uri']):
            property_uri = row['property_uri']
        elif ':' in str(row['property_id']) and str(row['property_id']).split(':')[0].lower() == 'skos':
            property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].split(':')[1]
        elif ':' in str(row['property_id']):
            try:
                property_uri = curie_dct[row['property_id'].split(':')[0].lower()] + row['property_id'].replace(':', '_')
            except KeyError:
                property_uri = None
                #print(f'There is a reference curie with an unrecognized namespace: {row["property_id"]}')
        else:
            property_uri = None
        statements.at[i, 'property_uri'] = property_uri
    
    # Save graph
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    os.makedirs(path, exist_ok=True)
    statements = statements[['subject_id', 'property_id', 'object_id', 'reference_uri',
                             'reference_supporting_text', 'reference_date', 'property_label',
                             'property_description', 'property_uri']]
    statements.fillna('NA').to_csv(f'{path}/graph_edges_v{today}.csv', index=False)
    
    print(f'\n* This is the size of the edges file data structure: {statements.shape}')
    print(f'* These are the edges attributes: {statements.columns}')
    print(f'* This is the first record:\n{statements.head(1)}')
    print(f'\nThe knowledge graph edges are built and saved at: {path}/graph_edges_v{today}.csv\n')
    print('\nFinished build_edges().\n')
    
    return statements

def build_nodes(statements, dgidb, monarch_disease, monarch_symptom, input_from_file=False):
    """
    Builds the nodes graph.
    :param statements: graph edges dataframe
    :param dgidb: path to DGIdb nodes file
    :param monarch_disease: path to Monarch disease nodes file
    :param monarch_symptom: path to Monarch symptom nodes file
    :param input_from_file: boolean to read from file
    :return: nodes dataframe
    """
    print('\nThe function "build_nodes()" is running...')
    
    if input_from_file:
        print('\nPreparing networks...')
        if dgidb and os.path.exists(dgidb):
            dgidb_df = pd.read_csv(dgidb)
        else:
            print(f"Warning: DGIdb nodes file not found at {dgidb}")
            dgidb_df = pd.DataFrame()
        
        if monarch_disease and os.path.exists(monarch_disease):
            monarch_dis_df = pd.read_csv(monarch_disease)
        else:
            print(f"Warning: Monarch disease nodes file not found at {monarch_disease}")
            monarch_dis_df = pd.DataFrame()
        
        if monarch_symptom and os.path.exists(monarch_symptom):
            monarch_symp_df = pd.read_csv(monarch_symptom)
        else:
            print(f"Warning: Monarch symptom nodes file not found at {monarch_symptom}")
            monarch_symp_df = pd.DataFrame()
    else:
        monarch_dis_df = pd.DataFrame(monarch_disease)
        monarch_symp_df = pd.DataFrame(monarch_symptom)
        dgidb_df = pd.DataFrame(dgidb)
    
    print('Monarch:')
    print(f'Disease nodes: {monarch_dis_df.shape}, Columns: {monarch_dis_df.columns}')
    print(f'Symptom nodes: {monarch_symp_df.shape}, Columns: {monarch_symp_df.columns}')
    print('DGIdb:')
    print(f'Nodes: {dgidb_df.shape}, Columns: {dgidb_df.columns}')
    
    # Extract nodes from edges
    print('\nAnnotating nodes in the graph...')
    st_nodes_l = pd.concat([statements.subject_id, statements.object_id], ignore_index=True)
    st_nodes_l.drop_duplicates(inplace=True)
    st_nodes_df = pd.DataFrame({'id': st_nodes_l})
    print(f'Graph nodes from edges: {st_nodes_df.shape}')
    
    # Annotate nodes
    monarch_dis_nodes = pd.merge(monarch_dis_df, st_nodes_df, how='inner', on='id')
    monarch_symp_nodes = pd.merge(monarch_symp_df, st_nodes_df, how='inner', on='id')
    dgidb_nodes = pd.merge(dgidb_df, st_nodes_df, how='inner', on='id')
    
    print('Annotation check:')
    print(f'Monarch disease nodes: {monarch_dis_nodes.shape}')
    print(f'Monarch symptom nodes: {monarch_symp_nodes.shape}')
    print(f'DGIdb nodes: {dgidb_nodes.shape}')
    
    # Concatenate all nodes
    print('\nConcatenating all nodes...')
    nodes = pd.concat([monarch_dis_nodes, monarch_symp_nodes, dgidb_nodes], ignore_index=True, join="inner")
    print(f'Combined nodes: {nodes.shape}')
    diff = set(st_nodes_df.id) - set(nodes.id)
    print(f'Nodes in edges but not annotated: {diff}')
    
    # Handle synonyms and drop duplicates
    print('\nDrop duplicated rows...')
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    nodes.drop_duplicates(keep='first', inplace=True)
    print(f'After dropping duplicates: {nodes.shape}')
    
    print('\nDrop duplicated nodes...')
    nodes.drop_duplicates(subset=['id'], keep='first', inplace=True)
    print(f'After dropping duplicate nodes: {nodes.shape}')
    
    # Check annotation completeness
    if len(set(st_nodes_df.id)) != len(set(nodes.id)):
        print('\nThere is a problem in the annotation of nodes.')
        print(f'Monarch disease nodes not in graph: {set(monarch_dis_df.id) - set(monarch_dis_nodes.id)}')
        print(f'Monarch symptom nodes not in graph: {set(monarch_symp_df.id) - set(monarch_symp_nodes.id)}')
        print(f'DGIdb nodes not in graph: {set(dgidb_df.id) - set(dgidb_nodes.id)}')
    else:
        print('\nAll graph nodes are annotated.')
    
    # Save nodes
    print('\nSaving final graph...')
    path = os.getcwd() + "/graph"
    os.makedirs(path, exist_ok=True)
    nodes = nodes[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    nodes['synonyms'] = nodes.synonyms.apply(lambda x: str('|'.join(x)) if isinstance(x, list) else x)
    # added for synthesizing.
    nodes['uri'] = nodes.apply(lambda row: row['uri'] if is_real_uri(row['uri']) else (curie_to_uri(row['id']) or 'NA'), axis=1)
       
    nodes.fillna('NA').to_csv(f'{path}/graph_nodes_v{today}.csv', index=False)
    
    print(f'\n* This is the size of the nodes file data structure: {nodes.shape}')
    print(f'* These are the nodes attributes: {nodes.columns}')
    print(f'* This is the first record:\n{nodes.head(1)}')
    print(f'\nThe knowledge graph nodes are built and saved at: {path}/graph_nodes_v{today}.csv\n')
    print('\nFinished build_nodes().\n')
    
    return nodes

def run_combine_graphs(date, disease_folder=None):
    """
    Runs the combine graphs script and saves nodes and edges files.
    :param date: date of the graph
    :param disease_folder: disease folder name (optional)
    :return: None (saves CSVs)
    """
    # Adjust base path for disease-specific folder
    base_path = os.getcwd()
    # if disease_folder:
    #     base_path = os.path.join(base_path, 'drugapp', 'data', disease_folder)
    
    # File paths
    drugsim_file = os.path.join(os.getcwd(), 'similaritygraph', f'drugdrugsim_v{today}.csv')
    dgidb_edges_file = os.path.join(os.getcwd(), 'DGIdb', f'DGIdb_edges_v{today}.csv')
    dgidb_nodes_file = os.path.join(os.getcwd(), 'DGIdb', f'DGIdb_nodes_v{today}.csv')
    monarch_edges_file = os.path.join(base_path, 'monarch', f'monarch_edges_disease_v{date}.csv')
    monarch_nodes_file = os.path.join(base_path, 'monarch', f'monarch_nodes_disease_v{date}.csv')
    monarch_edges_symptom_file = os.path.join(base_path, 'monarch', f'monarch_edges_symptom_v{date}.csv')
    monarch_nodes_symptom_file = os.path.join(base_path, 'monarch', f'monarch_nodes_symptom_v{date}.csv')
    
    # Build network
    print("Preparing networks...")
    statements = build_edges(
        dgidb=dgidb_edges_file,
        monarch_disease=monarch_edges_file,
        monarch_symptom=monarch_edges_symptom_file,
        drugsim=drugsim_file,
        input_from_file=True
    )
    build_nodes(
        statements,
        dgidb=dgidb_nodes_file,
        monarch_disease=monarch_nodes_file,
        monarch_symptom=monarch_nodes_symptom_file,
        input_from_file=True
    )

if __name__ == '__main__':
    run_combine_graphs(date=str(today))
    
import requests
import os
import datetime
import pandas as pd
from biothings_client import get_client
from tqdm import tqdm
import math

# timestamp
today = datetime.date.today()
curr_year = int(str(today)[:4])

def is_real_uri(value):
    """Return True if value is a real-looking URI string (not 'NA', 'nan', empty, or NaN)."""
    if value is None:
        return False
    if isinstance(value, float) and math.isnan(value):
        return False
    if not isinstance(value, str):
        return False
    v = value.strip()
    if v == "" or v.lower() in {"na", "nan", "none", "n/a"}:
        return False
    # Simple heuristic: must contain a scheme or at least a ':' (like http: or identifiers.org: or http)
    if v.startswith("http") or ":" in v:
        return True
    return False

def curie_to_uri(curie: str) -> str:
    """Convert a CURIE like 'HGNC:1234' to a valid URI string, or return None if impossible."""
    if not isinstance(curie, str) or ':' not in curie:
        return None
    prefix, local = curie.split(':', 1)
    prefix_l = prefix.lower()
    local = local.strip()
    if not local:
        return None
    if prefix_l == 'hgnc':
        return f'https://identifiers.org/HGNC:{local}'
    if prefix_l == 'mondo':
        return f'http://purl.obolibrary.org/obo/MONDO_{local}'
    if prefix_l in {'hp', 'hpo'}:
        return f'http://purl.obolibrary.org/obo/HP_{local}'
    if prefix_l in {'go', 'pato', 'geno'}:
        return f'http://purl.obolibrary.org/obo/{prefix.upper()}_{local}'
    # fallback to identifiers.org for any other CURIE prefix
    return f'https://identifiers.org/{prefix}:{local}'

def read_connections(filename):
    """
    This function reads monarch_connections CSV file.
    :param filename: complete path to the monarch connections csv file string
    :return: monarch edges dataframe
    """

    # monarch network
    path = os.getcwd() + '/monarch'
    csv_path = path + '/' + filename
    network_df = pd.read_csv(csv_path)
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))

    return network_df

def hit_monarch_api(node='HGNC:4851', rows=500, category=None, category_in=None, use_entity=True):
    """
    This function performs API calls to Monarch to retrieve out and in edges from a query node.
    If use_entity=True and node is a disease, uses entity API for phenotypes/genes, else uses association API.
    :param node: node id to query (string). Default: 'HGNC:4851' (=HTT gene).
    :param rows: the maximum number of results to return (integer). Default: 500.
    :param category: association category for out (string). Default: None.
    :param category_in: association category for in (string). Default: None.
    :param use_entity: whether to use entity API for disease nodes (boolean). Default: True.
    :return: two API response objects: 'out' and 'in' response objects
    """
    if use_entity and node.startswith('MONDO:'):
        # Use entity API for disease nodes
        entity_url = f'https://api.monarchinitiative.org/v3/api/entity/{node}'
        r_entity = requests.get(entity_url)
        if r_entity.status_code != 200:
            print(f"Entity API failed for {node}: {r_entity.status_code}")
            return {'items': []}, {'items': []}
        json_data = r_entity.json()
        
        # r_out: Disease → Phenotypes (has_phenotype)
        r_out_items = []
        phenotypes = json_data.get('has_phenotype', [])
        phenotype_labels = json_data.get('has_phenotype_label', [])
        disease_name = json_data.get('name', 'NA')
        for i, pheno_id in enumerate(phenotypes):
            r_out_items.append({
                'subject': node,
                'subject_label': disease_name,
                'object': pheno_id,
                'object_label': phenotype_labels[i] if i < len(phenotype_labels) else 'NA',
                'predicate': 'RO:0002200',
                'subject_category': 'biolink:Disease',
                'object_category': 'biolink:PhenotypicFeature',
                'publications': ['NA']
            })
        
        # r_in: Genes → Disease (causal_gene)
        r_in_items = []
        causal_genes = json_data.get('causal_gene', [])
        for gene in causal_genes:
            r_in_items.append({
                'subject': gene.get('id', 'NA'),
                'subject_label': gene.get('name', 'NA'),
                'object': node,
                'object_label': disease_name,
                'predicate': 'RO:0003303',
                'subject_category': 'biolink:Gene',
                'object_category': 'biolink:Disease',
                'publications': ['NA']
            })
        
        # Debug: Print number of phenotypes and genes found
        print(f"Entity API for {node}: {len(r_out_items)} phenotypes, {len(r_in_items)} genes")
        return {'items': r_out_items}, {'items': r_in_items}
    
    # Fallback to association API
    biolink = 'https://api.monarchinitiative.org/v3/api/association'
    parameters = {'limit': rows}
    if category:
        parameters['category'] = category
    
    out_params = {'subject': node, **parameters}
    r_out = requests.get(biolink, params=out_params)
    
    parameters_in = {'limit': rows}
    if category_in:
        parameters_in['category'] = category_in
    in_params = {'object': node, **parameters_in}
    r_in = requests.get(biolink, params=in_params)
    
    r_out_json = r_out.json() if r_out.status_code == 200 else {'items': []}
    r_in_json = r_in.json() if r_in.status_code == 200 else {'items': []}
    print(f"Association API for {node}: {len(r_out_json.get('items', []))} out, {len(r_in_json.get('items', []))} in")
    return r_out_json, r_in_json

def get_disease_name_id(disease_input_ID='OMIM:143100'):
    """
    This function finds the disease name and ID (MONDO preferred, fallback to OMIM) using v3 search/entity APIs.
    :param disease_input_ID: The input ID of the disease
    :return: two strings, name and ID respectively
    """
    search_url = 'https://api.monarchinitiative.org/v3/api/search'
    params = {'q': disease_input_ID, 'category': 'biolink:Disease', 'limit': 1}
    r_search = requests.get(search_url, params=params)
    items = r_search.json().get('items', [])
    
    if not items:
        return 'NA', disease_input_ID
    
    mondo_id = items[0].get('id', 'NA')
    entity_url = f'https://api.monarchinitiative.org/v3/api/entity/{mondo_id}'
    r_entity = requests.get(entity_url)
    json_data = r_entity.json()
    disease_name = json_data.get('name', 'NA')
    disease_id = json_data.get('id', disease_input_ID)
    return disease_name, disease_id

def get_edges_objects(r_out, r_in):
    """
    This function prepares the API object responses from Monarch.
    :param r_out: BioLink API 'out' response object
    :param r_in: BioLink API 'in' response object
    :return: subjects, relations, objects, and references lists
    """
    sub_l, rel_l, obj_l, ref_l = [], [], [], []
    
    for assoc_list in [r_out.get('items', []), r_in.get('items', [])]:
        for association in assoc_list:
            if not isinstance(association, dict):
                continue
            sub_l.append({
                'id': association.get('subject', 'NA'),
                'label': association.get('subject_label', 'NA'),
                'iri': 'NA',
                'category': [association.get('subject_category', 'NA')]
            })
            rel_l.append({
                'id': association.get('predicate', 'NA'),
                'label': association.get('predicate', 'NA').split(':')[-1].replace('_', ' '),
                'iri': 'NA'
            })
            obj_l.append({
                'id': association.get('object', 'NA'),
                'label': association.get('object_label', 'NA'),
                'iri': 'NA',
                'category': [association.get('object_category', 'NA')]
            })
            pubs = association.get('publications', [])
            ref_l.append('|'.join(pubs if pubs else ['NA']))
    
    return sub_l, rel_l, obj_l, ref_l

def get_edges(sub_l, rel_l, obj_l, ref_l, attribute='id'):
    edges = set()
    for i in range(len(sub_l)):
        sub = sub_l[i][attribute]
        rel = rel_l[i][attribute]
        obj = obj_l[i][attribute]
        ref = ref_l[i]
        edges.add((sub, rel, obj, ref))
    return edges

def add_attributes(sub_l, rel_l, obj_l, edges):
    metaedges = set()
    for (sub_id, rel_id, obj_id, refs) in edges:
        for i in range(len(sub_l)):
            if sub_l[i]['id'] == sub_id and rel_l[i]['id'] == rel_id and obj_l[i]['id'] == obj_id:
                metaedges.add((
                    sub_id,
                    sub_l[i]['label'],
                    'NA',
                    sub_l[i]['category'][0],
                    rel_id,
                    rel_l[i]['label'],
                    'NA',
                    obj_id,
                    obj_l[i]['label'],
                    'NA',
                    obj_l[i]['category'][0],
                    refs
                ))
                break
    return metaedges

def keep_node_type(edges, seed, nodeType='ortho'):
    propertyList = ['RO:HOM0000017', 'RO:HOM0000020'] if nodeType == 'ortho' else ['RO:0002200', 'RO:0002607', 'RO:0002326', 'GENO:0000840']
    keep = set()
    for (sub, rel, obj, ref) in edges:
        if rel in propertyList:
            if sub not in seed:
                keep.add(sub)
            if obj not in seed:
                keep.add(obj)
    return keep

def get_connections(nodes):
    keep = set()
    for node in tqdm(nodes, desc="Processing nodes"):
        try:
            r_out, r_in = hit_monarch_api(node)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            meta_edges = add_attributes(sub_l, rel_l, obj_l, edges)
            keep.update(meta_edges)
        except Exception as e:
            print(f"Error processing node {node}: {e}")
    print(f"Total edges extracted: {len(keep)}")
    return keep
    
def get_neighbours_list(seed_list):
    neighbours = set()
    for node in tqdm(seed_list):
        r_out, r_in = hit_monarch_api(node, category='biolink:DiseaseToPhenotypicFeatureAssociation', category_in='biolink:GeneToDiseaseAssociation', use_entity=True)
        sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
        edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
        neighbours.update([obj for (sub, rel, obj, ref) in edges])
    return list(neighbours)

def get_orthopheno_list(seed_list):
    orthopheno = set()
    for node in tqdm(seed_list):
        r_out, r_in = hit_monarch_api(node, category='biolink:orthologous_to', use_entity=False)
        sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
        edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
        orthologs = keep_node_type(edges, seed_list, 'ortho')
        for ortholog in orthologs:
            r_out, r_in = hit_monarch_api(ortholog, category='biolink:has_phenotype', use_entity=False)
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            phenotypes = keep_node_type(edges, [ortholog], 'pheno')
            orthopheno.update(phenotypes)
    return list(orthopheno)

def extract_edges(gene_list):
    keep = set()
    for node in tqdm(gene_list, desc="Processing nodes"):
        try:
            r_out, r_in = hit_monarch_api(node, use_entity=node.startswith('MONDO:'))
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            meta_edges = add_attributes(sub_l, rel_l, obj_l, edges)
            keep.update(meta_edges)
        except Exception as e:
            print(f"Error processing node {node}: {e}")
    print(f"Total edges extracted: {len(keep)}")
    if not keep:
        print("Warning: No edges found for provided nodes. Network is empty.")
    return keep

def print_network(network, filename):
    if not network:
        print(f"!!!THIS IS PRINT NETWORK. Warning: No edges to save for {filename}_v{today}.csv. DataFrame is empty.")
        return None
    edges = [
        {
            'subject_id': tup[0],
            'subject_label': tup[1],
            'subject_uri': tup[2],
            'subject_category': tup[3],
            'relation_id': tup[4],
            'relation_label': tup[5],
            'relation_uri': tup[6],
            'object_id': tup[7],
            'object_label': tup[8],
            'object_uri': tup[9],
            'object_category': tup[10],
            'reference_id_list': tup[11]
        }
        for tup in network
    ]
    df = pd.DataFrame(edges).fillna('NA')
    path = os.getcwd() + '/monarch'
    os.makedirs(path, exist_ok=True)
    df.to_csv(f'{path}/{filename}_v{today}.csv', index=False)
    return df

# BUILD NETWORK

def build_edges(edges_df, edges_fname):
    """
    This function builds the edges network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph edges object as a list of dictionaries, where every dictionary is a record
    """

    print('\nThe function "build_edges()" is running...')
    ## variables
    # if edges_df is not a df, it is a set of tuples, then convert to a df
    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['subject_uri'] = tuple[2]
            record['subject_category'] = tuple[3]
            record['relation_id'] = tuple[4]
            record['relation_label'] = tuple[5]
            record['relation_uri'] = tuple[6]
            record['object_id'] = tuple[7]
            record['object_label'] = tuple[8]
            record['object_uri'] = tuple[9]
            record['object_category'] = tuple[10]
            record['reference_id_list'] = tuple[11]
            connections_l.append(record)

        edges_df = pd.DataFrame(connections_l)


    # generate static variable: uriPrefixes_dct (url references)
    uriPrefixes_dct = {
        'pmid': 'https://www.ncbi.nlm.nih.gov/pubmed/',  
        'react': 'http://reactome.org/content/detail/',  
        'zfin': 'http://zfin.org/',
        'go_ref': 'http://purl.obolibrary.org/obo/go/references/',  
        'mgi': 'http://www.informatics.jax.org/accession/MGI:',  
        'flybase': 'http://flybase.org/reports/',
        'wormbase': 'http://www.wormbase.org/resources/paper/',
        'hpo': 'http://compbio.charite.de/hpoweb/showterm?id=HP:',
        'isbn-10': 'ISBN-10:',
        'isbn-13': 'ISBN-13:',
        'mondo': 'http://purl.obolibrary.org/obo/MONDO_',  
        'rgd': 'https://rgd.mcw.edu/rgdweb/report/reference/main.html?id=', 
        'omim': 'http://omim.org/entry/',  
        'sgd_ref': 'https://db.yeastgenome.org/reference/',  
        'genereviews': 'https://www.ncbi.nlm.nih.gov/books/',  
        'omia': 'http://omia.angis.org.au/',  
        'hgnc': 'https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:', 
    }

    # generate static variable: dbPrefixes_dct (source/database references)
    dbPrefixes_dct = {
        'na': 'NA',
        'nan': 'NA',
        'mgi': 'http://www.informatics.jax.org/',
        'fb': 'http://flybase.org/',
        'rgd': 'http://rgd.mcw.edu/',
        'zfin': 'http://zfin.org/',
        'sgd': 'https://www.yeastgenome.org/',
        'hgnc': 'https://www.genenames.org/',
        'xenbase': 'http://www.xenbase.org/'
    }

    # provenance variables
    ref_text = 'This edge comes from the Monarch Initiative {}.'.format(curr_year) 
    ref_date = today

    ## build graph schema network edges data structure and save edges file
    # prepare dataframe = [ {} ... {} ], where every row = {} = concept
    edges_l = list()
    for edge in edges_df.itertuples():
        # edge or row is a tuple (named and ordered attributes)
        # edge.reference_id_list >> can be 1) np.nan (float type) or 2) str without "|" 3) str with "|"
        ref_s = str(edge.reference_id_list)

        # prepare reference_uri_list attribute
        ref_uri_l = list()
        # expand to uri or NA
        pmid_l = list()
        # reference_id list iteration
        for ref in ref_s.strip().split('|'):
            # NA or database
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("Warning:")
                    print('Detected a new reference database in Monarch not yet implemented in this module. '
                          'The new source should be added to the dictionary of databases.'
                          'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                    print("In the build_edges() method, update 'dbPrefixes_dct' dictionary with '{}'".format(ref))
                    print('The edge that includes this new reference database is {}'.format(edge))
                    print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                          'until the dictionary is updated.')
                ref_uri_l.append(ref_uri)
            # publication uri: pubmed_id or url
            else:
                pref, uriId = ref.split(':')
                # separate pmid from non pmid and detect:
                # pubmed_id
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                # url
                elif ref.lower().startswith('http'):
                    ref_uri_l.append(ref)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("Warning:")
                        print('Detected a new reference source in Monarch not yet implemented in this module. '
                              'The new source should be added to the dictionary of sources.'
                              'Otherwise, the source CURIE cannot be translated to the corresponding URI.')
                        print("In the build_edges() method, update 'uriPrefixes_dct' dictionary with '{}'".format(pref))
                        print('The edge that includes this new reference source is {}'.format(edge))
                        print('The method will continue to run without problem, writing the CURIE instead of the URI,'
                              'until the dictionary is updated.')
                    ref_uri_l.append(ref_uri)
        # create multi-term pubmed url
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)

        # prepare edge attributes: sub_id, obj_id, rel_id, rel_label, rel_def, rel_uri
        sub_id = 'NA' if edge.subject_id is None or str(edge.subject_id) == 'nan' or str(edge.subject_id) == 'None' else edge.subject_id
        rel_id = 'NA' if edge.relation_id is None or str(edge.relation_id) == 'nan'  or str(edge.relation_id) == 'None' else edge.relation_id
        obj_id = 'NA' if edge.object_id is None or str(edge.object_id) == 'nan'  or str(edge.object_id) == 'None' else edge.object_id
        rel_label = 'NA' if edge.relation_label is None or str(edge.relation_label) == 'nan'  or str(edge.relation_label) == 'None' else edge.relation_label
        rel_uri = 'NA' if edge.relation_uri is None or str(edge.relation_uri) == 'nan'  or str(edge.relation_uri) == 'None' else edge.relation_uri
        rel_def = 'NA'
        
        # make sure there are no 'NA' relation ids
        if rel_id == 'NA':
            if edge.object_category == 'genotype': # if object is genotype, relation is 'has genotype'
                rel_id = 'GENO:0000222'
                rel_label = 'has_genotype'
                rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000222'
            else:
                rel_id = 'RO:0002610' # else (gene, model objects), relation is 'correlated with'
                rel_label = 'correlated with'
                rel_uri = 'http://purl.obolibrary.org/obo/RO_0002610'
        
        # change every 'contributes to' to 'contributes to condition'
        if rel_id == 'RO:0002326':
            rel_id = 'RO:0003304'
            rel_label = 'contributes to condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003304'
        
        # change every 'is causal germline mutation in' and 'is causal germline mutation partially giving rise to' and 'pathogenic for condition' to 'causes condition'
        if rel_id == 'RO:0004013' or rel_id == 'RO:0004016' or rel_id == 'GENO:0000840':
            rel_id = 'RO:0003303'
            rel_label = 'causes condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003303'

        # change every 'in orthology relationship with' to 'in 1 to 1 orthology relationship with'
        if rel_id == 'RO:HOM0000017':
            rel_id = 'RO:HOM0000020'
            rel_label = 'in 1 to 1 orthology relationship with'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_HOM0000020'  
        
        # if 'genotype' -->'has phenotype' --> 'disease', change 'has phenotype' to 'has role in modelling'
        if edge.subject_category == 'genotype' and edge.object_category == 'disease' and rel_id == 'RO:0002200':
            rel_id = 'RO:0003301'
            rel_label = 'has role in modeling'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003301'
        
        # change every 'is reference allele of' to 'is allele of'
        if rel_id == 'GENO:0000610':
            rel_id = 'GENO:0000408'
            rel_label = 'is_allele_of'
            rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000408'
        
        # build the data structure = list of edges as list of dict, where a dict is an edge
        edge = dict()
        edge['subject_id'] = sub_id
        edge['object_id'] = obj_id
        edge['property_id'] = rel_id
        edge['property_label'] = rel_label
        edge['property_description'] = rel_def
        edge['property_uri'] = rel_uri
        edge['reference_uri'] = ref_uri_list
        edge['reference_supporting_text'] = ref_text
        edge['reference_date'] = ref_date
        edges_l.append(edge)

    # save edges file
    df = pd.DataFrame(edges_l)
    print('df',df.shape)
    path = os.getcwd() + '/monarch'
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', \
             'property_label', 'property_description', 'property_uri']]
    df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,edges_fname,today), index=False)

    # print info
    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe Monarch network edges are built and saved at: {}/monarch_edges_v{}.csv\n'.format(path,today))
    print('\nFinished build_edges().\n')

    return df

def build_nodes(edges_df, nodes_fname):
    """
    This function builds the nodes network with the graph schema.
    :param edges_df: network dataframe from the extract_edges() function
    :return: graph nodes object as a list of dictionaries, where every dictionary is a record
    """
    print('\nThe function "build_nodes()" is running...')
    concept_dct = dict()

    if isinstance(edges_df, set):
        connections_l = list()
        for tuple in edges_df:
            record = dict()
            record['subject_id'] = tuple[0]
            record['subject_label'] = tuple[1]
            record['subject_uri'] = tuple[2]
            record['subject_category'] = tuple[3]
            record['relation_id'] = tuple[4]
            record['relation_label'] = tuple[5]
            record['relation_uri'] = tuple[6]
            record['object_id'] = tuple[7]
            record['object_label'] = tuple[8]
            record['object_uri'] = tuple[9]
            record['object_category'] = tuple[10]
            record['reference_id_list'] = tuple[11]
            connections_l.append(record)
        edges_df = pd.DataFrame(connections_l)

    for edge in edges_df.itertuples():
        sid = edge.subject_id 
        oid = edge.object_id 
        concept_dct[sid] = {}
        concept_dct[oid] = {}
    print('Number of concepts: {}'.format(len(concept_dct.keys())))

    concept_dct = dict()

    biolink_to_semantic = {
        'biolink:AnatomicalEntity': 'anatomical entity',
        'biolink:Disease': 'disease',
        'biolink:MolecularActivity': 'molecular function',
        'biolink:BiologicalProcess': 'biological process',
        'biolink:Genotype': 'genotype',
        'biolink:Gene': 'gene',
        'biolink:PhenotypicFeature': 'phenotype',
        'biolink:SequenceVariant': 'variant',
        'biolink:Pathway': 'pathway',
        'biolink:CellularComponent': 'cellular component',
    }

    for edge in edges_df.itertuples():
        sid = edge.subject_id 
        slab = edge.subject_label 
        sgroup = edge.subject_category
        suri = edge.subject_uri
        oid = edge.object_id 
        olab = edge.object_label 
        ogroup = edge.object_category
        ouri = edge.object_uri

        # Handle NaN
        sgroup = 'NA' if pd.isna(sgroup) else biolink_to_semantic.get(sgroup, sgroup.lower().replace('biolink:', ''))
        ogroup = 'NA' if pd.isna(ogroup) else biolink_to_semantic.get(ogroup, ogroup.lower().replace('biolink:', ''))
        
        if sgroup == 'genotype' and ogroup == 'gene' and edge.relation_id == 'GENO:0000408':
            sgroup = 'variant'

        # Normalize incoming URIs/placeholders and synthesize when missing
        # use suri/ouri only if they look like real URIs; otherwise synthesize from CURIE (id)
        subject_uri_val = suri if is_real_uri(suri) else curie_to_uri(sid)
        object_uri_val  = ouri if is_real_uri(ouri) else curie_to_uri(oid)
        
        concept_dct[sid] = {'preflabel': slab or 'NA',
                            'semantic_groups': sgroup,
                            'uri': subject_uri_val if subject_uri_val is not None else 'NA' ,
                            'synonyms': 'NA', 
                            'description': 'NA'}
        concept_dct[oid] = {'preflabel': olab or 'NA',
                            'semantic_groups': ogroup,
                            'uri': object_uri_val if object_uri_val is not None else 'NA',
                            'synonyms': 'NA', 
                            'description': 'NA'}

    # ... (rest unchanged)

    # biothings: annotate name,synonyms,description to genes
    print('\nAdding BioThings annotation: name, synonyms, description...')
    # input: (preflabel) symbol,alias
    symbols = list()
    for concept in concept_dct:
        if isinstance(concept_dct[concept]['semantic_groups'], list):
            for label in concept_dct[concept]['semantic_groups']:
                if 'gene' in label:
                    symbols.append(concept_dct[concept]['preflabel'])
        else:
            if 'gene' in concept_dct[concept]['semantic_groups']:
                symbols.append(concept_dct[concept]['preflabel'])
    
    print('symbols:', len(symbols))

    # query biothings
    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name,alias,summary', size=1, as_dataframe=True)

    # dictionary: {symbol:name}
    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
    ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
    monarch_s2n = dict(zip(ids.symbol, ids.name))
    monarch_s2s = dict(zip(ids.symbol, ids.synonyms))
    monarch_s2d = dict(zip(ids.symbol, ids.description))

    # prepare data structure = [ {} ... {} ], where every {} = concept = row
    nodes_l = list()
    for concept in concept_dct:
        # define nodes (rows) for the data structure
        preflabel = concept_dct[concept]['preflabel']
        concept_dct[concept]['synonyms'] = monarch_s2s[preflabel] if preflabel in monarch_s2s.keys() else 'NA'
        node = dict()
        node['id'] = concept
        node['semantic_groups'] = concept_dct[concept]['semantic_groups']
        node['uri'] = concept_dct[concept]['uri']
        node['preflabel'] = preflabel
        node['name'] = monarch_s2n[preflabel] if preflabel in monarch_s2n.keys() else preflabel
        node['synonyms'] = '|'.join(list(concept_dct[concept]['synonyms'])) if isinstance(
            concept_dct[concept]['synonyms'], list) else concept_dct[concept]['synonyms']
        node['description'] = monarch_s2d[preflabel] if preflabel in monarch_s2d.keys() else 'NA'
        nodes_l.append(node)

    # save nodes file
    df = pd.DataFrame(nodes_l)
    # Ensure 'uri' is a string and replace None / nan-like with 'NA'
    df['uri'] = df['uri'].apply(lambda x: x if is_real_uri(x) else 'NA')
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/monarch'
    df.fillna('NA').to_csv('{}/{}_v{}.csv'.format(path,nodes_fname, today), index=False)

    # print info
    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe Monarch network nodes are built and saved at: {}/monarch_nodes_v{}.csv\n'.format(path,today))
    print('\nFinished build_nodes().\n')

    return df

def get_symptoms_disease(disease_id, monarch_edges_csv, monarch_nodes_csv):
    """
    This function finds all symptoms of a disease.
    :param disease_id: id of the disease of which symptoms and drugs should be found (e.g. 'MONDO:0007739' for HD)  
    :param monarch_edges_csv: name of the edges csv file from Monarch
    :param monarch_nodes_csv: name of the nodes csv file from Monarch
    :return: a list of symptom names
    """
    # open csv file
    edges_df = pd.read_csv(monarch_edges_csv)
    nodes_df = pd.read_csv(monarch_nodes_csv)

    # find all nodes that are one step away from disease and have as relation 'has phenotype' (='RO:0002200')
    df_symptoms = edges_df.loc[(edges_df['subject_id'] == disease_id) & (edges_df['property_id'] == 'RO:0002200')]
    print(f"Symptom edges for {disease_id}: {len(df_symptoms)}")

    symptoms_id_lst = df_symptoms['object_id'].tolist()
    #get names of symptoms
    symptoms_name_lst = nodes_df.loc[nodes_df['id'].isin(symptoms_id_lst), 'preflabel'].to_list()
    print(f"Symptoms found: {symptoms_name_lst}")
    return symptoms_name_lst

def symptom_list_from_folder(input_folder):
    """
    This function finds the symptom list of the disease in the specified folder.
    :param input_folder: folder name containing disease_id.txt and monarch CSVs (e.g., 'Alzheimer disease type 1 (2025-10-15)')
    :return: symptom name list, date of creation, disease name with date
    """
    import os
    date = input_folder.split('(')[-1].rstrip(')')
    
    # Change to the specified folder
    base_path = os.getcwd()
    folder_path = os.path.join(base_path, 'drugapp', 'data', input_folder)
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"Folder {folder_path} does not exist")
    os.chdir(folder_path)
    
    # Read disease_id.txt
    with open("disease_id.txt", 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')
    
    disease_name_date = f"{disease_name} ({date})"
    monarch_fname_edges = f'./monarch/monarch_edges_disease_v{date}.csv'
    monarch_fname_nodes = f'./monarch/monarch_nodes_disease_v{date}.csv'
    
    # Get list of symptoms
    symptoms_name_lst = get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)
    
    # Restore original working directory
    os.chdir(base_path)
    
    return symptoms_name_lst, date, disease_name_date

def symptom_list_today():
    """
    This function finds the symptom list of the disease specified by user on 
    that same day.
    :return: symptom name list, date of creation of files required, disease name
    """
    date=today
    
    # find the disease ID that was previously written in a txt file
    with open(os.getcwd()+"/disease_id.txt", 'r') as file:
        file_text = file.read()
        input_number, disease_id, disease_name = file_text.split(';')

    disease_name_date = disease_name+' ({})'.format(date)
    
    monarch_fname_edges = './monarch/monarch_edges_disease_v{}.csv'.format(date)
    monarch_fname_nodes = './monarch/monarch_nodes_disease_v{}.csv'.format(date)
    
    # get list of symptoms for disease to ask user
    symptoms_name_lst=get_symptoms_disease(disease_id, monarch_fname_edges, monarch_fname_nodes)
    
    return symptoms_name_lst, date, disease_name_date

def run_monarch(input_number='104300'):
    input_id = 'OMIM:' + input_number
    disease_name, disease_id = get_disease_name_id(input_id)
    new_path = os.getcwd() + '/drugapp/data/' + disease_name + ' ({})'.format(today)
    os.makedirs(new_path, exist_ok=True)
    os.chdir(new_path)
    text_file = open("disease_id.txt", "w")
    text_file.write(input_number + ';' + disease_id + ';' + disease_name)
    text_file.close()
    
    seed_list = [disease_id]
    neighbour_list = get_neighbours_list(seed_list)
    orthopheno_list = get_orthopheno_list(seed_list)
    gene_list = seed_list + neighbour_list + orthopheno_list
    network = extract_edges(gene_list)
    df = print_network(network, 'monarch_orthopeno_network_disease')
    if df is None:
        print("Error: No network data to build edges and nodes.")
        return
    build_edges(df, 'monarch_edges_disease')
    build_nodes(df, 'monarch_nodes_disease')

def run_monarch_symptom(input_symptom, date, disease_folder):
    """
    This function runs the Monarch script using the symptom ID as the seed and saves nodes/edges files.
    :param input_symptom: The input symptom name (e.g., 'Dementia')
    :param date: Date string of the disease folder (e.g., '2025-10-15')
    :param disease_folder: Folder name containing disease_id.txt and CSVs (e.g., 'Alzheimer disease type 1 (2025-10-15)')
    :return: None (saves CSVs in monarch folder)
    """
    # Change to disease folder
    base_path = os.getcwd()
    folder_path = os.path.join(base_path, )# 'drugapp', 'data', disease_folder)
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"Folder {folder_path} does not exist")
    os.chdir(folder_path)
    
    # Get symptom ID from nodes CSV
    monarch_fname_nodes = f'./monarch/monarch_nodes_disease_v{date}.csv'
    try:
        nodes_df = pd.read_csv(monarch_fname_nodes)
        symptom_id_row = nodes_df[nodes_df['preflabel'] == input_symptom]
        if symptom_id_row.empty:
            raise ValueError(f"Symptom '{input_symptom}' not found in {monarch_fname_nodes}")
        symptom_id = symptom_id_row['id'].iloc[0]
    except Exception as e:
        os.chdir(base_path)
        raise Exception(f"Error finding symptom ID for '{input_symptom}': {e}")
    
    # Build Monarch network with symptom as seed
    seed_list = [symptom_id]
    neighbour_list = get_neighbours_list(seed_list)
    orthopheno_list = get_orthopheno_list(seed_list)
    gene_list = seed_list + neighbour_list + orthopheno_list
    network = extract_edges(gene_list)

    print(f"Building symptom KG for symptom: {input_symptom} (ID: {symptom_id})")
    print(f"Neighbours: {len(neighbour_list)}, Orthophenotypes: {len(orthopheno_list)}")
    print(f"Total edges: {len(network)}")
    
    # Save network
    df = print_network(network, 'monarch_orthopeno_network_symptom')
    if df is None:
        os.chdir(base_path)
        raise ValueError("No edges generated for symptom network")
    
    # Build and save edges/nodes
    build_edges(df, edges_fname='monarch_edges_symptom')
    build_nodes(df, nodes_fname='monarch_nodes_symptom')
    
    # Restore original working directory
    os.chdir(base_path)

#run_monarch()
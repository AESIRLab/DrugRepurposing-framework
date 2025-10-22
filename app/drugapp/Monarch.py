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
    return f'https://identifiers.org/{prefix}:{local}'

def read_connections(filename):
    """
    This function reads monarch_connections CSV file.
    :param filename: complete path to the monarch connections csv file string
    :return: monarch edges dataframe
    """
    path = os.getcwd() + '/monarch'
    csv_path = path + '/' + filename
    network_df = pd.read_csv(csv_path)
    print('\n* This is the size of the data structure: {}'.format(network_df.shape))
    print('* These are the attributes: {}'.format(network_df.columns))
    print('* This is the first record:\n{}'.format(network_df.head(1)))
    return network_df

def hit_monarch_api(node='HGNC:4851', rows=500, category=None, category_in=None, use_entity=True):
    r_out_items = []
    r_in_items = []
    r_out_json, r_in_json = {'items': []}, {'items': []}

    if use_entity and node.startswith('MONDO:'):
        entity_url = f'https://api.monarchinitiative.org/v3/api/entity/{node}'
        r_entity = requests.get(entity_url)
        if r_entity.status_code != 200:
            print(f"Entity API failed for {node}: {r_entity.status_code}")
            return {'items': []}, {'items': []}
        json_data = r_entity.json()
        
        disease_name = json_data.get('name', 'NA')
        for pheno_id, pheno_label in zip(json_data.get('has_phenotype', []), json_data.get('has_phenotype_label', [])):
            r_out_items.append({
                'subject': node,
                'subject_label': disease_name,
                'object': pheno_id,
                'object_label': pheno_label,
                'predicate': 'biolink:has_phenotype',
                'subject_category': 'biolink:Disease',
                'object_category': 'biolink:PhenotypicFeature',
                'publications': ['NA']
            })
        
        # Restore original category for MONDO genes
        assoc_url = f'https://api.monarchinitiative.org/v3/api/association?object={node}&category=biolink:CausalGeneToDiseaseAssociation&limit={rows}'
        r_assoc = requests.get(assoc_url)
        if r_assoc.status_code == 200:
            assoc_data = r_assoc.json().get('items', [])
            for assoc in assoc_data:
                r_in_items.append({
                    'subject': assoc.get('subject', 'NA'),
                    'subject_label': assoc.get('subject_label', 'NA'),
                    'object': node,
                    'object_label': disease_name,
                    'predicate': 'biolink:causes_or_contributes_to',
                    'subject_category': assoc.get('subject_category', 'biolink:Gene'),
                    'object_category': 'biolink:Disease',
                    'publications': assoc.get('publications_links', ['NA'])
                })
        
        print(f"Entity API for {node}: {len(r_out_items)} phenotypes, {len(r_in_items)} genes")
    else:
        # For phenotype seeds (e.g., HP:), directly query incoming genes
        if node.startswith(('HP:', 'MP:', 'ZP:')):
            assoc_url = f'https://api.monarchinitiative.org/v3/api/association?object={node}&category=biolink:GeneToPhenotypicFeatureAssociation&limit={rows}'
            r_assoc = requests.get(assoc_url)
            if r_assoc.status_code == 200:
                assoc_data = r_assoc.json().get('items', [])
                for assoc in assoc_data:
                    r_in_items.append({
                        'subject': assoc.get('subject', 'NA'),  # Gene
                        'subject_label': assoc.get('subject_label', 'NA'),
                        'object': node,
                        'object_label': assoc.get('object_label', 'NA'),
                        'predicate': assoc.get('predicate', 'biolink:causes_or_contributes_to'),
                        'subject_category': assoc.get('subject_category', 'biolink:Gene'),
                        'object_category': 'biolink:PhenotypicFeature',
                        'publications': assoc.get('publications_links', ['NA'])
                    })
            print(f"Association API for phenotype {node}: {len(r_in_items)} gene links")

        # Default association API calls (unchanged)
        base_url = 'https://api.monarchinitiative.org/v3/api/association'
        params_out = {'subject': node, 'limit': rows}
        if category: params_out['category'] = category
        r_out = requests.get(base_url, params=params_out)
        params_in = {'object': node, 'limit': rows}
        if category_in: params_in['category'] = category_in
        r_in = requests.get(base_url, params=params_in)
        
        r_out_json = r_out.json() if r_out.status_code == 200 else {'items': []}
        r_in_json = r_in.json() if r_in.status_code == 200 else {'items': []}
        r_out_items.extend(r_out_json.get('items', []))
        r_in_items.extend(r_in_json.get('items', []))
        print(f"API for {node}: {len(r_out_items)} out, {len(r_in_items)} in (cat out: {category}, in: {category_in})")

    return {'items': r_out_items}, {'items': r_in_items}

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
    biolink_to_ro = {
        'biolink:has_phenotype': 'RO:0002200',
        'biolink:causes_or_contributes_to': 'RO:0003303',
        'biolink:orthologous_to': 'RO:HOM0000017',
        'biolink:gene_to_phenotypic_feature_association': 'RO:0002327'
    }
    for assoc_list in [r_out.get('items', []), r_in.get('items', [])]:
        for association in assoc_list:
            if not isinstance(association, dict): continue
            pred = association.get('predicate', 'NA')
            mapped_pred = biolink_to_ro.get(pred, pred)
            sub_l.append({
                'id': association.get('subject', 'NA'),
                'label': association.get('subject_label', 'NA'),
                'iri': 'NA',
                'category': [association.get('subject_category', 'NA')]
            })
            rel_l.append({
                'id': mapped_pred,
                'label': mapped_pred.split(':')[-1].replace('_', ' '),
                'iri': 'NA'
            })
            obj_l.append({
                'id': association.get('object', 'NA'),
                'label': association.get('object_label', 'NA'),
                'iri': 'NA',
                'category': [association.get('object_category', 'NA')]
            })
            # Handle publications safely
            pubs = association.get('publications_links', association.get('publications', ['NA']))
            if isinstance(pubs, dict):
                pubs = [pubs.get('id', 'NA')] if 'id' in pubs else ['NA']  # Extract 'id' from single dict
            elif isinstance(pubs, list):
                pubs = [item.get('id', 'NA') if isinstance(item, dict) and 'id' in item else str(item) for item in pubs]  # Handle list of dicts or mixed
            else:
                pubs = [str(pubs)] if pubs else ['NA']
            ref_l.append('|'.join(pubs) if pubs else 'NA')
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
                # Create 12-tuple with NA for missing fields
                metaedges.add((
                    sub_id,
                    sub_l[i].get('label', 'NA'),
                    sub_l[i].get('iri', 'NA'),
                    sub_l[i].get('category', ['NA'])[0],
                    rel_id,
                    rel_l[i].get('label', rel_id.split(':')[-1].replace('_', ' ')),
                    rel_l[i].get('iri', 'NA'),
                    obj_id,
                    obj_l[i].get('label', 'NA'),
                    obj_l[i].get('iri', 'NA'),
                    obj_l[i].get('category', ['NA'])[0],
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
    print(f"Processing all {len(nodes)} nodes with edge filtering")
    # Add this mapping (biolink → RO, from original build_edges logic + Monarch docs)
    biolink_to_ro = {
        'biolink:has_phenotype': 'RO:0002200',
        'biolink:causes_or_contributes_to': 'RO:0003303',  # Or RO:0003304 for contributes
        'biolink:orthologous_to': 'RO:HOM0000017',  # 1:1 ortho
        # Add more as needed, e.g., 'biolink:gene_to_phenotypic_feature_association' → 'RO:0002327' (has phenotype qualifier)
    }
    all_edges = []  # Accumulate before filtering
    for node in tqdm(nodes, desc="Processing nodes"):
        try:
            r_out, r_in = hit_monarch_api(node, use_entity=node.startswith('MONDO:'))
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            # Temp: Collect all, then filter post-loop for efficiency
            for (sub, rel, obj, ref) in edges:
                mapped_rel = biolink_to_ro.get(rel, rel)  # Fallback to raw if no map
                all_edges.append((sub, mapped_rel, obj, ref))
        except Exception as e:
            print(f"Error processing node {node}: {e}")

    print(f'[Monarch.py][get_connections], node in the get_connections: {nodes[:5]}, {len(nodes)}')
    print(f'[Monarch.py][get_connections], printing the All_edges before meta_edges: {all_edges[:5]}, {len(all_edges)}')
    
    # Filter post-collection: Direct rels only (RO equivs), cap at 1000
    rel_filter = ['RO:0002200', 'RO:0003303', 'RO:HOM0000017']  # has_phenotype, causes_condition, orthologous
    meta_edges = [(sub, rel, obj, ref) for (sub, rel, obj, ref) in all_edges if rel in rel_filter][:1000]
    print(f'[Monarch.py][get_connections], meta_edges: {meta_edges[:5]}, {len(meta_edges)}')
    
    # Build 12-tuples (add_attributes style, but inline for simplicity)
    for (sub_id, rel_id, obj_id, refs) in meta_edges:
        # Find matching sub/rel/obj attrs (simplified; assumes sub_l/rel_l/obj_l from last iter—improve by passing full lists if needed)
        # For now, placeholder attrs; enhance with global lookup if edges stay low
        sub_label, sub_iri, sub_cat = 'NA', curie_to_uri(sub_id) or 'NA', 'biolink:Gene' if sub_id.startswith('HGNC:') else 'biolink:PhenotypicFeature'
        obj_label, obj_iri, obj_cat = 'NA', curie_to_uri(obj_id) or 'NA', 'biolink:PhenotypicFeature' if obj_id.startswith(('HP:', 'MP:', 'ZP:')) else 'biolink:Disease'
        rel_label = rel_id.split(':')[-1].replace('_', ' ')
        rel_iri = 'http://purl.obolibrary.org/obo/' + rel_id.replace(':', '_')
        keep.add((
            sub_id, sub_label, sub_iri, sub_cat,
            rel_id, rel_label, rel_iri,
            obj_id, obj_label, obj_iri, obj_cat,
            refs
        ))
    
    keep = list(keep)[:1000]  # Final cap
    print(f"Total edges extracted: {len(keep)}")
    if not keep:
        print("Warning: No edges found for provided nodes. Network is empty.")
    return keep

def get_neighbours_list(seed_list, max_depth=1, max_nodes_per_layer=20):  # Keep your BFS, but broad calls
    neighbours = set(seed_list)
    current_layer = set(seed_list)
    for depth in range(max_depth):
        next_layer = set()
        for node in tqdm(current_layer, desc=f"Depth {depth+1}"):
            r_out, r_in = hit_monarch_api(node, rows=500)  # Broad! No category
            sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
            edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
            # Original-style filter: Add keep_nodes here for junk removal
            next_layer.update([obj for (sub, rel, obj, ref) in edges if rel in ['RO:0002200', 'RO:0003303'] and obj not in neighbours])
        next_layer -= neighbours
        if len(next_layer) > max_nodes_per_layer: next_layer = set(list(next_layer)[:max_nodes_per_layer])
        neighbours.update(next_layer)
        current_layer = next_layer
        if not current_layer: break
    return list(neighbours - set(seed_list))

def get_orthologs(gene_list):
    """
    Fetches orthologs for a list of genes using GeneToGeneHomologyAssociation.
    :param gene_list: List of gene IDs (e.g., ['HGNC:620']).
    :return: List of ortholog IDs.
    """
    orthologs = set()
    for gene in tqdm(gene_list, desc="Fetching orthologs"):
        if gene.startswith('HGNC:'):
            assoc_url = f'https://api.monarchinitiative.org/v3/api/association?subject={gene}&category=biolink:GeneToGeneHomologyAssociation&limit=50'  # Cap per gene
            r_assoc = requests.get(assoc_url)
            if r_assoc.status_code == 200:
                assoc_data = r_assoc.json().get('items', [])
                orthologs.update([assoc.get('object', 'NA') for assoc in assoc_data if assoc.get('object', '').startswith(('MGI:', 'ZFIN:', 'FB:'))])  # Model organisms only
    print(f"Found {len(orthologs)} orthologs for {len(gene_list)} genes")
    return list(orthologs)

def get_orthopheno_list(seed_list):
    orthopheno = set()
    genes = set()
    for node in seed_list:
        # Broad for disease seed (no category_in—original style)
        r_out, r_in = hit_monarch_api(node)  # Drop category_in here!
        sub_l, rel_l, obj_l, ref_l = get_edges_objects(r_out, r_in)
        edges = get_edges(sub_l, rel_l, obj_l, ref_l, 'id')
        genes.update([sub for (sub, rel, obj, ref) in edges if rel == 'RO:0003303'])  # Post-filter causes
    print(f"Extracted {len(genes)} genes from seed: {list(genes)[:5]}...")
    
    if not genes:
        print("No genes found; skipping ortholog expansion.")
        return []
    
    orthologs = get_orthologs(list(genes))
    for ortholog in tqdm(orthologs, desc="Fetching orthophenotypes"):
        r_out_pheno, r_in_pheno = hit_monarch_api(ortholog, category='biolink:GeneToPhenotypicFeatureAssociation')
        print(f"API response for {ortholog}: {len(r_out_pheno['items']) + len(r_in_pheno['items'])} items")  # Debug
        sub_l_pheno, rel_l_pheno, obj_l_pheno, ref_l_pheno = get_edges_objects(r_out_pheno, r_in_pheno)
        edges_pheno = get_edges(sub_l_pheno, rel_l_pheno, obj_l_pheno, ref_l_pheno, 'id')
        phenotypes = [obj for (sub, rel, obj, ref) in edges_pheno if rel == 'RO:0002200' and obj.startswith(('HP:', 'MP:', 'ZP:'))]  # RO post-map
        phenotypes = phenotypes[:10]  # Tighter cap
        orthopheno.update(phenotypes)
        print(f"Added {len(phenotypes)} phenotypes from {ortholog}")
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

    ref_text = 'This edge comes from the Monarch Initiative {}.'.format(curr_year) 
    ref_date = today

    edges_l = list()
    for edge in edges_df.itertuples():
        ref_s = str(edge.reference_id_list)
        ref_uri_l = list()
        pmid_l = list()
        for ref in ref_s.strip().split('|'):
            if ':' not in ref:
                try:
                    ref_uri = dbPrefixes_dct[ref.lower()]
                except KeyError:
                    print("Warning: New reference database '{}'. Add to dbPrefixes_dct.".format(ref))
                    ref_uri = ref
                ref_uri_l.append(ref_uri)
            else:
                pref, uriId = ref.split(':')
                if ref.startswith('PMID'):
                    pmid_l.append(uriId)
                elif ref.lower().startswith('http'):
                    ref_uri_l.append(ref)
                else:
                    try:
                        ref_uri = uriPrefixes_dct[pref.lower()] + uriId
                    except KeyError:
                        print("Warning: New reference source '{}'. Add to uriPrefixes_dct.".format(pref))
                        ref_uri = ref
                    ref_uri_l.append(ref_uri)
        if len(pmid_l):
            pmid_s = ','.join(pmid_l)
            ref_uri = uriPrefixes_dct['pmid'] + pmid_s
            ref_uri_l.append(ref_uri)
        ref_uri_list = '|'.join(ref_uri_l)

        sub_id = 'NA' if pd.isna(edge.subject_id) else edge.subject_id
        rel_id = 'NA' if pd.isna(edge.relation_id) else edge.relation_id
        obj_id = 'NA' if pd.isna(edge.object_id) else edge.object_id
        rel_label = 'NA' if pd.isna(edge.relation_label) else edge.relation_label
        rel_uri = 'NA' if pd.isna(edge.relation_uri) else edge.relation_uri
        rel_def = 'NA'
        
        if rel_id == 'NA':
            if edge.object_category == 'genotype':
                rel_id = 'GENO:0000222'
                rel_label = 'has_genotype'
                rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000222'
            else:
                rel_id = 'RO:0002610'
                rel_label = 'correlated with'
                rel_uri = 'http://purl.obolibrary.org/obo/RO_0002610'
        
        # Relation mappings (from original)
        if rel_id == 'RO:0002326':
            rel_id = 'RO:0003304'
            rel_label = 'contributes to condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003304'
        if rel_id in ['RO:0004013', 'RO:0004016', 'GENO:0000840']:
            rel_id = 'RO:0003303'
            rel_label = 'causes condition'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003303'
        if rel_id == 'RO:HOM0000017':
            rel_id = 'RO:HOM0000020'
            rel_label = 'in 1 to 1 orthology relationship with'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_HOM0000020'
        if edge.subject_category == 'genotype' and edge.object_category == 'disease' and rel_id == 'RO:0002200':
            rel_id = 'RO:0003301'
            rel_label = 'has role in modeling'
            rel_uri = 'http://purl.obolibrary.org/obo/RO_0003301'
        if rel_id == 'GENO:0000610':
            rel_id = 'GENO:0000408'
            rel_label = 'is_allele_of'
            rel_uri = 'http://purl.obolibrary.org/obo/GENO_0000408'

        edge_dict = {
            'subject_id': sub_id,
            'object_id': obj_id,
            'property_id': rel_id,
            'property_label': rel_label,
            'property_description': rel_def,
            'property_uri': rel_uri,
            'reference_uri': ref_uri_list,
            'reference_supporting_text': ref_text,
            'reference_date': ref_date
        }
        edges_l.append(edge_dict)

    df = pd.DataFrame(edges_l)
    print('df', df.shape)
    path = os.getcwd() + '/monarch'
    df = df[['subject_id', 'property_id', 'object_id', 'reference_uri', 'reference_supporting_text', 'reference_date', 
             'property_label', 'property_description', 'property_uri']]
    df.fillna('NA').to_csv(f'{path}/{edges_fname}_v{today}.csv', index=False)

    print('\n* This is the size of the edges file data structure: {}'.format(pd.DataFrame(edges_l).shape))
    print('* These are the edges attributes: {}'.format(pd.DataFrame(edges_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(edges_l).head(1)))
    print('\nThe Monarch network edges are built and saved at: {}/monarch_edges_v{}.csv\n'.format(path, today))
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
        slab = edge.subject_label if edge.subject_label != 'NA' else edge.subject_id 
        sgroup = edge.subject_category
        suri = edge.subject_uri
        oid = edge.object_id 
        olab = edge.object_label if edge.object_label != 'NA' else edge.object_id 
        ogroup = edge.object_category
        ouri = edge.object_uri

        # Handle NaN
        sgroup = 'NA' if pd.isna(sgroup) else biolink_to_semantic.get(sgroup, sgroup.lower().replace('biolink:', ''))
        ogroup = 'NA' if pd.isna(ogroup) else biolink_to_semantic.get(ogroup, ogroup.lower().replace('biolink:', ''))
        
        if sgroup == 'genotype' and ogroup == 'gene' and edge.relation_id == 'GENO:0000408':
            sgroup = 'variant'

        # Normalize incoming URIs/placeholders and synthesize when missing
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

    # biothings: annotate name,synonyms,description to genes
    print('\nAdding BioThings annotation: name, synonyms, description...')
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

    mg = get_client('gene')
    df = mg.querymany(symbols, scopes='symbol,alias', fields='name,alias,summary', size=1, as_dataframe=True)

    ids = (df.reset_index().rename(columns={'query': 'symbol'}))
    monarch_s2n = {}
    monarch_s2s = {}
    monarch_s2d = {}
    if 'alias' in ids.columns:
        ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
    else:
        ids['synonyms'] = 'NA'
    if 'summary' in ids.columns:
        ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
    else:
        ids['description'] = 'NA'
    if not ids.empty:
        if 'name' in ids.columns:
            monarch_s2n = dict(zip(ids.symbol, ids.name))
        if 'alias' in ids.columns:
            ids['synonyms'] = ids.alias.apply(lambda x: x if str(x) != 'nan' else 'NA')
            monarch_s2s = dict(zip(ids.symbol, ids.synonyms))
        if 'summary' in ids.columns:
            ids['description'] = ids.summary.apply(lambda x: x if str(x) != 'nan' else 'NA')
            monarch_s2d = dict(zip(ids.symbol, ids.summary))

    nodes_l = list()
    for concept in concept_dct:
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

    df = pd.DataFrame(nodes_l)
    df['uri'] = df['uri'].apply(lambda x: x if is_real_uri(x) else 'NA')
    df = df[['id', 'semantic_groups', 'uri', 'preflabel', 'name', 'synonyms', 'description']]
    path = os.getcwd() + '/monarch'
    df.fillna('NA').to_csv(f'{path}/{nodes_fname}_v{today}.csv', index=False)

    print('\n* This is the size of the nodes file data structure: {}'.format(pd.DataFrame(nodes_l).shape))
    print('* These are the nodes attributes: {}'.format(pd.DataFrame(nodes_l).columns))
    print('* This is the first record:\n{}'.format(pd.DataFrame(nodes_l).head(1)))
    print('\nThe Monarch network nodes are built and saved at: {}/monarch_nodes_v{}.csv\n'.format(path, today))
    print('\nFinished build_nodes().\n')

    return df

def get_pheno_label(pheno_id):
    """Quick label fetch for pheno IDs."""
    if not pheno_id.startswith('HP:'): return pheno_id
    url = f'https://api.monarchinitiative.org/v3/api/entity/{pheno_id}'
    r = requests.get(url)
    if r.status_code == 200:
        return r.json().get('name', pheno_id)
    return pheno_id

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
    symptoms_name_lst = [get_pheno_label(id) for id in symptoms_id_lst]  # Fetch labels!
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
    with open("disease_id.txt", "w") as text_file:
        text_file.write(input_number + ';' + disease_id + ';' + disease_name)
    
    seed_list = [disease_id]
    print(f'[Monarch.py][run_monarch] seed_list: {seed_list}')
    neighbour_list = get_neighbours_list(seed_list, max_depth=1, max_nodes_per_layer=20)
    orthopheno_list = get_orthopheno_list(seed_list)
    gene_list = seed_list + neighbour_list + orthopheno_list
    print(f"Full gene_list size: {len(gene_list)} (seed: {len(seed_list)}, neighbors: {len(neighbour_list)}, orthopheno: {len(orthopheno_list)})")
    network = get_connections(gene_list)
    df = print_network(network, 'monarch_orthopheno_network_disease')
    if df is None:
        print("Error: No network data to build edges and nodes.")
        return
    build_edges(df, 'monarch_edges_disease')
    build_nodes(df, 'monarch_nodes_disease')

def run_monarch_symptom(input_symptom, date, disease_folder):
    base_path = os.getcwd()
<<<<<<< HEAD
    folder_path = os.path.join(base_path) #'drugapp', 'data', disease_folder)
=======
    folder_path = os.path.join(base_path, )# 'drugapp', 'data', disease_folder)
>>>>>>> 6c7d80776fdcc54cd34da8d6e6f98affd13a59f9
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"Folder {folder_path} does not exist")
    os.chdir(folder_path)
    
    monarch_fname_nodes = f'./monarch/monarch_nodes_disease_v{date}.csv'
    try:
        nodes_df = pd.read_csv(monarch_fname_nodes)
        symptom_row = nodes_df[nodes_df['preflabel'] == input_symptom]
        if symptom_row.empty:
            # Fallback: Search IDs for label match (if user picked name, reverse)
            id_map = {get_pheno_label(row['id']): row['id'] for _, row in nodes_df.iterrows()}
            symptom_id = id_map.get(input_symptom, input_symptom)  # ID if name not found
        else:
            symptom_id = symptom_row['id'].iloc[0]
    except Exception as e:
        os.chdir(base_path)
        raise Exception(f"Error finding symptom ID for '{input_symptom}': {e}")
    
    seed_list = [symptom_id]
    neighbour_list = get_neighbours_list(seed_list)
    orthopheno_list = get_orthopheno_list(seed_list)
    gene_list = seed_list + neighbour_list + orthopheno_list
    network = extract_edges(gene_list)

    print(f"Building symptom KG for symptom: {input_symptom} (ID: {symptom_id})")
    print(f"Neighbours: {len(neighbour_list)}, Orthophenotypes: {len(orthopheno_list)}")
    print(f"Total edges: {len(network)}")
    
    df = print_network(network, 'monarch_orthopheno_network_symptom')
    if df is None:
        os.chdir(base_path)
        raise ValueError("No edges generated for symptom network")
    
    build_edges(df, edges_fname='monarch_edges_symptom')
    build_nodes(df, nodes_fname='monarch_nodes_symptom')
    
    os.chdir(base_path)

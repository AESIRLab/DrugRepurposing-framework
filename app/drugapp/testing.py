import pandas as pd
date = "2025-10-15"
fn = f'./graph/graph_nodes_v{date}.csv'
df = pd.read_csv(fn)
print("Total nodes:", len(df))
print(df.semantic_groups.value_counts(dropna=False))
print("URIs null / NaN / 'NA':")
print(df.uri.isna().sum(), (df.uri.astype(str).str.lower() == 'na').sum())
print("\nSample gene nodes (first 20):")
print(df[df.semantic_groups=='gene'].head(20))
print("\nGenes without valid URI (sample 50):")
bad_genes = df[(df.semantic_groups=='gene') & (~df.uri.astype(str).str.strip().astype(bool) | df.uri.isna() | (df.uri.astype(str).str.lower().isin(['na','nan','none'])))]
print(len(bad_genes))
print(bad_genes[['id','preflabel','uri']].head(50))

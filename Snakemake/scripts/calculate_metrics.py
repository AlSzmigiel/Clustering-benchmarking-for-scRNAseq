import scanpy as sc
from sklearn.metrics import silhouette_score, adjusted_rand_score, adjusted_mutual_info_score, calinski_harabasz_score
import genieclust
import pandas as pd 
import re
from sklearn import preprocessing 

dataset = sc.read(str(snakemake.input))
partitions = snakemake.params.partition #this will allow to
print(partitions)
#use regex to match observations in adata and retrieve different clusterings.
#this has to be done because in datasets there is different metadata, stored in different ways, 
#deoending if it was Robject originally or not 
pattern = re.compile('|'.join(partitions))
print(pattern)
expression_matches = [entry for entry in list(dataset.obs.columns) if pattern.match(entry)]
label_key = snakemake.params.label_key
label_encoder = preprocessing.LabelEncoder() 
label_true = label_encoder.fit_transform(dataset.obs[label_key]) 
scores = {}
for column_name, clustering_column in dataset.obs[expression_matches].items():
    eval = genieclust.compare_partitions.compare_partitions(label_true, clustering_column)
    #eval['modularity'] = igraph.Graph.Weighted_Adjacency(dataset.obsp["connectivities"]).modularity(dataset.obs[column_name].cat.codes)
    scores[column_name] = eval


#silhouette_scores= {column_name: silhouette_score(dataset.obsp['distances'],clustering_column) 
              #for column_name, clustering_column in dataset.obs[expression_matches].items()}
scores_df = pd.DataFrame.from_dict(scores, orient='index')
#silhouette_df = pd.DataFrame.from_dict(silhouette_scores, orient='index', columns=['Silhouette Score'])



#merged_df = pd.concat([adjusted_df, silhouette_df], axis=1)

scores_df.to_csv(str(snakemake.output))

#f = open(snakemake.output.summary_metrics, 'w')
#writer = csv.writer(f)
#writer.writerow([str(snakemake.wildcards.dataset),silhouette])





import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster


# Read the data from the phenotypes.tsv file
phenotypes_data = pd.read_csv('phenotypes.tsv', delimiter='\t')
# Create a DataFrame where each row is a gene and its associated phenotypes
genes_phenotypes = phenotypes_data.groupby('Gene')['Phenotype Description'].apply(list).reset_index()
# Populate genelist
gene_names = genes_phenotypes['Gene'].tolist()

# Get a sorted list of unique phenotypes
unique_phenotypes = sorted(phenotypes_data['Phenotype Description'].unique())

# Create a dictionary to map phenotype to index
phenotype_to_index = {phenotype: index for index, phenotype in enumerate(unique_phenotypes)}

# Initialize the adjacency matrix
adjacency_matrix = np.zeros((len(unique_phenotypes), len(unique_phenotypes)), dtype=int)

# Populate the adjacency matrix
for phenotypes in genes_phenotypes['Phenotype Description']:
    for i in range(len(phenotypes)):
        for j in range(i + 1, len(phenotypes)):
            index_i = phenotype_to_index[phenotypes[i]]
            index_j = phenotype_to_index[phenotypes[j]]
            # Since it's undirected, we add to both [i, j] and [j, i]
            adjacency_matrix[index_i][index_j] += 1
            adjacency_matrix[index_j][index_i] += 1

# Create the heatmap
plt.figure(figsize=(12, 12))
sns.heatmap(adjacency_matrix, cmap='viridis')
plt.title('Heatmap of Phenotype Co-occurrence')
plt.xlabel('Phenotypes')
plt.ylabel('Phenotypes')
plt.xticks(ticks=np.arange(len(unique_phenotypes)), labels=unique_phenotypes, rotation=90)
plt.yticks(ticks=np.arange(len(unique_phenotypes)), labels=unique_phenotypes, rotation=0)
plt.savefig('hpo_heatmap.png')

# Perform hierarchical clustering
row_linkage = linkage(adjacency_matrix, method='ward')
col_linkage = linkage(adjacency_matrix.T, method='ward')

# Create a clustered heatmap with seaborn
sns.clustermap(adjacency_matrix, cmap='viridis',
               row_linkage=row_linkage, col_linkage=col_linkage,
               figsize=(12, 12))

# Save the clustered heatmap
plt.savefig('clustered_heatmap.png')

# Perform hierarchical clustering
linkage_matrix = linkage(adjacency_matrix, method='ward')

# Determine the clusters, here we use a threshold to cut the dendrogram
# The 'maxclust' parameter can be set to specify the maximum number of clusters
cluster_ids = fcluster(linkage_matrix, t=10, criterion='maxclust')

# Map each gene to its cluster
gene_to_cluster = {gene: cluster_id for gene, cluster_id in zip(gene_names, cluster_ids)}

# Create a DataFrame from this mapping
clusters_df = pd.DataFrame(list(gene_to_cluster.items()), columns=['Gene', 'ClusterID'])

# Group by 'ClusterID' and count the number of genes in each cluster
cluster_sizes = clusters_df.groupby('ClusterID').size().reset_index(name='Count')

# Sort clusters by size in descending order
sorted_clusters = cluster_sizes.sort_values(by='Count', ascending=False)

# Get the top 10 largest clusters
top_clusters = sorted_clusters.head(10)['ClusterID'].tolist()

# Now, we extract the genes for each of the top clusters and write to a text file
with open('top_clusters_genes.txt', 'w') as f:
    for cluster_id in top_clusters:
        # Get the genes in this cluster
        genes_in_cluster = clusters_df[clusters_df['ClusterID'] == cluster_id]['Gene'].tolist()
        
        # Write the cluster id and the genes to the file
        f.write(f"Cluster {cluster_id} consists of genes: {', '.join(genes_in_cluster)}\n\n")

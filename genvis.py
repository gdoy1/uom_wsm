import pandas as pd
import igraph as ig
import numpy as np
import csv
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from scipy.spatial.distance import squareform
import requests
from collections import Counter

# Load the graph from the GraphML file
file_path = 'epilepsy.graphml'
graph = ig.Graph.Read_GraphML(file_path)

# Calculate degree centrality for each gene
degrees = graph.degree()

# Combine gene names with their degrees
gene_names = graph.vs['display name']  # Grabs gene symbol from 'name' GraphML object
gene_with_degrees = zip(gene_names, degrees)

# Sort genes by their degree centrality
sorted_genes = sorted(gene_with_degrees, key=lambda x: x[1], reverse=True)

# Print the top 10 most important genes
print("Top 10 most important genes based on degree centrality:")
for gene, degree in sorted_genes[:10]:
    print(f"{gene}: {degree}")

# Calculate and set vertex size based on degree centrality:
vertex_sizes = [vertex.degree() * 10 for vertex in graph.vs]

# Visualizing the graph
visual_style = {
    "vertex_size": vertex_sizes,
    "bbox": (800, 800),
    "margin": 100
}
igraph_plot = ig.plot(graph, **visual_style)
igraph_plot.save('gene_network.png')  # Saves the plot to a file

# Generating the adjacency matrix from the graph
adjacency_matrix = np.array(graph.get_adjacency().data)

# Creating a heatmap from the adjacency matrix
plt.figure(figsize=(10, 10))
sns.heatmap(adjacency_matrix, cmap='viridis')
plt.title('Heatmap of Gene Interactions')
plt.xlabel('Genes')
plt.ylabel('Genes')
plt.savefig('gene_heatmap.png')  # Saves the heatmap to a file

# Squaring the adjacency matrix
squared_adjacency_matrix = np.matmul(adjacency_matrix, adjacency_matrix)

# Ensure the squared matrix is symmetric
squared_adjacency_matrix = (squared_adjacency_matrix + squared_adjacency_matrix.T) / 2

# Normalize the squared adjacency matrix
max_value = np.max(squared_adjacency_matrix)
normalized_squared_matrix = squared_adjacency_matrix / max_value

# Creating a heatmap from the normalized squared adjacency matrix
plt.figure(figsize=(10, 10))
sns.heatmap(normalized_squared_matrix, cmap='viridis', vmin=0, vmax=1)
plt.title('Normalized Heatmap of Gene Interactions (Squared Adjacency Matrix)')
plt.xlabel('Genes')
plt.ylabel('Genes')
plt.savefig('normalized_squared_gene_heatmap.png')  # Saves the heatmap to a file

# Sum the values in each row of the squared adjacency matrix
two_step_paths = np.sum(squared_adjacency_matrix, axis=1)

# Combine gene names with their two-step path counts
gene_with_two_step_paths = zip(gene_names, two_step_paths)

# Sort genes by their two-step path counts
sorted_genes_by_two_step_paths = sorted(gene_with_two_step_paths, key=lambda x: x[1], reverse=True)

# Print the top 10 most centrally connected genes based on two-step paths
print("Top 10 most centrally connected genes based on two-step paths:")
for gene, count in sorted_genes_by_two_step_paths[:10]:
    print(f"{gene}: {count}")

# Convert the normalized squared adjacency matrix to a distance matrix
distance_matrix = 1 - normalized_squared_matrix

# Set the diagonal elements to zero
np.fill_diagonal(distance_matrix, 0)

# Convert to a condensed distance matrix (required for clustering)
condensed_distance_matrix = squareform(distance_matrix)

# Perform hierarchical clustering
linkage_matrix = linkage(condensed_distance_matrix, 'ward')

# Plot the dendrogram
plt.figure(figsize=(20, 10))
dendrogram(linkage_matrix, labels=gene_names)
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('Gene')
plt.ylabel('Distance')
plt.savefig('gene_clustering_dendrogram.png')

# Plot degree distribution
degrees = graph.degree()
plt.figure()
plt.hist(degrees, bins=30)
plt.title('Degree Distribution')
plt.xlabel('Degree')
plt.ylabel('Frequency')
plt.savefig('degree_distribution.png')

# Check if the graph is directed, and if so, convert it to undirected for community detection
if graph.is_directed():
    undirected_graph = graph.as_undirected()
else:
    undirected_graph = graph

# Community detection using the leading eigenvector method
communities = graph.community_leading_eigenvector()
membership = communities.membership

# Map each gene to its community
gene_to_community = {gene: community for gene, community in zip(gene_names, membership)}

# Calculate the average degree centrality for each community
community_avg_centrality = {}
for gene, community in gene_to_community.items():
    if community not in community_avg_centrality:
        community_avg_centrality[community] = []
    community_avg_centrality[community].append(degrees[gene_names.index(gene)])

community_avg_centrality = {community: np.mean(centralities) for community, centralities in community_avg_centrality.items()}

# Identify the community with the highest average centrality
most_central_community = max(community_avg_centrality, key=community_avg_centrality.get)

# Extract the genes from the most central community
most_central_genes = [gene for gene, community in gene_to_community.items() if community == most_central_community]

print(f"Most central community (Community {most_central_community}) contains the following genes:")
print(most_central_genes)

with open('most_central_community_genes.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Gene'])  # Header
    for gene in most_central_genes:
        writer.writerow([gene])

print("Gene list of the most central community saved to 'most_central_community_genes.csv'")

# Get the order of leaves in the dendrogram
dendrogram_leaves = leaves_list(linkage_matrix)

# Reorder the distance matrix according to the leaves of the dendrogram
reordered_distance_matrix = distance_matrix[dendrogram_leaves, :][:, dendrogram_leaves]

# Creating a clustered heatmap
plt.figure(figsize=(20, 20))
sns.heatmap(reordered_distance_matrix, cmap='viridis', 
            xticklabels=np.array(gene_names)[dendrogram_leaves], 
            yticklabels=np.array(gene_names)[dendrogram_leaves])
plt.title('Clustered Heatmap of Gene Distance Matrix')
plt.xlabel('Genes')
plt.ylabel('Genes')
plt.savefig('clustered_gene_distance_heatmap.png')

# Calculate shortest paths
path_lengths = graph.shortest_paths_dijkstra()

# Flatten the path lengths array while filtering out infinite values
path_lengths_flattened = [path for sublist in path_lengths for path in sublist if path != 0 and np.isfinite(path)]

# Now, create the histogram plot for path lengths
plt.figure()
plt.hist(path_lengths_flattened, bins=30)
plt.title('Path Length Distribution')
plt.xlabel('Path Length')
plt.ylabel('Frequency')
plt.savefig('path_length_distribution.png')

plt.figure()
sns.heatmap([degrees], cmap='viridis')
plt.title('Node Importance Heatmap')
plt.xlabel('Genes')
plt.ylabel('Importance')
plt.savefig('node_importance_heatmap.png')

def get_phenotype_descriptions_for_gene(gene_symbol, all_descriptions):
    try:
        url = f"https://grch37.rest.ensembl.org/phenotype/gene/homo_sapiens/{gene_symbol}?content-type=application/json"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            phenotype_descriptions = [entry['description'] for entry in data if 'description' in entry]
            all_descriptions.update(phenotype_descriptions)
        else:
            print(f"Error: Received status code {response.status_code}")
    except requests.RequestException as e:
        print(f"An error occurred: {e}")

description_count = Counter()  # Initialize a Counter object to count descriptions

for gene in most_central_genes:
    print(f"Gene: {gene}")
    # Pass the Counter object as the second argument
    get_phenotype_descriptions_for_gene(gene, description_count)
    if description_count:
        print("Associated Phenotype Descriptions:")
        for description in description_count:
            print(description)
    else:
        print("No phenotype description available")
    print("\n")

def save_phenotypes_to_file(description_count, filename='phenotypes.txt'):
    with open(filename, 'w') as file:
        for description, count in description_count.items():
            file.write(f"{description}: {count}\n")
    print(f"Phenotype descriptions saved to {filename}")

# Accumulate phenotype descriptions for each gene
for gene in most_central_genes:
    get_phenotype_descriptions_for_gene(gene, description_count)

# Save the complete list of phenotype descriptions to a file
save_phenotypes_to_file(description_count)
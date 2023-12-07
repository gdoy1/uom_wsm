import requests
import csv
import pygtex
import json
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from collections import Counter
import re
import os


def get_genes_from_panel(panel_id):
    url = f"https://panelapp.genomicsengland.co.uk/api/v1/panels/{panel_id}/"
    try:
        response = requests.get(url)
        response.raise_for_status()
        panel_data = response.json()
        panel_genes = {gene['gene_data']['hgnc_symbol']: gene['confidence_level'] for gene in panel_data['genes']}
        return panel_genes
    except requests.HTTPError as http_err:
        print(f"HTTP error occurred: {http_err}")
    except Exception as err:
        print(f"An error occurred: {err}")

def get_kegg_pathways_for_gene(gene):
    url = f"https://rest.kegg.jp/find/genes/{gene}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.split('\n')
    else:
        return []

def extract_hsa_pathway_ids(pathway_data):
    pathway_ids = []
    for line in pathway_data:
        match = re.search(r"hsa:\d+", line)
        if match:
            pathway_id = match.group()
            pathway_ids.append(pathway_id)
    return pathway_ids


def read_gene_list(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        return [row[0] for row in reader]

def compare_genes_with_panel(panel_genes, most_central_genes):
    gene_status = {}
    for gene in most_central_genes:
        status = panel_genes.get(gene, 'Not Present')
        gene_status[gene] = status
    return gene_status

def get_gencode_id_for_gene(gene_symbol):
    try:
        gModel = pygtex.GeneModel(gene_symbol)
        return gModel.getGencodeId()
    except Exception as e:
        print(f"Error getting GENCODE ID for {gene_symbol}: {e}")
        return None

def get_eqtl_data(gencode_id, tissue):
    try:
        response = requests.get('https://gtexportal.org/api/v2/association/singleTissueEqtl',
                                params={"gencodeId": gencode_id, "tissueSiteDetailId": tissue, "datasetId": "gtex_v8"})
        response.raise_for_status()
        return json.loads(response.text)['data']
    except requests.HTTPError as http_err:
        print(f"HTTP error occurred while querying eQTL data for {gencode_id}: {http_err}")
    except Exception as err:
        print(f"An error occurred while querying eQTL data for {gencode_id}: {err}")

def write_to_csv(gene_status, gene_pvalue_dict, filename='updated_gene_status.csv'):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Gene', 'Status', 'Average P-Value'])
        for gene, status in gene_status.items():
            average_pvalue = gene_pvalue_dict.get(gene)
            writer.writerow([gene, status, average_pvalue])
    print(f"Updated gene status with average p-values saved to {filename}")

def calculate_average_pvalue(eqtl_data):
    if not eqtl_data:
        return None
    p_values = [item['pValue'] for item in eqtl_data if 'pValue' in item]
    return sum(p_values) / len(p_values) if p_values else None

def read_gene_status(filename):
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header
        return {row[0]: row[1] for row in reader}

def save_snp_list(eqtl_data, filename='snp_list.tsv'):
    with open(filename, 'a') as file:
        for item in eqtl_data:
            snp_id = item.get('snpId')
            p_value = item.get('pValue')
            if snp_id and p_value is not None:
                file.write(f"{snp_id}\t{p_value}\n")

# Example usage
panel_id = 'R59'
input_filename = 'most_central_community_genes.csv'
status_filename = 'updated_gene_status.csv'
tissue = 'Brain_Frontal_Cortex_BA9'
groups = ['1', '2', '3', 'Not Present']

most_central_genes = read_gene_list(input_filename)
panel_genes = get_genes_from_panel(panel_id)
gene_status = compare_genes_with_panel(panel_genes, most_central_genes)


gene_pvalue_dict = {}
gene_status_dict = read_gene_status(status_filename)
pathway_counter = Counter()

color_map = {
    "3": "green",
    "2": "#FFBF00",  # Hex code for amber
    "1": "red",
    "Not Present": "black"
}

for gene in most_central_genes:
    print('Analysing ' + gene)
    # Grab the network data
    pathway_data = get_kegg_pathways_for_gene(gene)
    pathways = extract_hsa_pathway_ids(pathway_data)
    pathway_counter.update(pathways)
    # Get the gencode ID
    gencode_id = get_gencode_id_for_gene(gene)
    print(gencode_id)
    if gencode_id:
        eqtl_data = get_eqtl_data(gencode_id, tissue)
        if eqtl_data:
            save_snp_list(eqtl_data)
        average_pvalue = calculate_average_pvalue(eqtl_data)
        if average_pvalue is not None:
            gene_pvalue_dict[gene] = average_pvalue

# Saving to pathway CSV
with open('pathways.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(['Pathway ID', 'Count'])
    for pathway, count in pathway_counter.items():
        writer.writerow([pathway, count])

# Plotting
for gene, pvalue in gene_pvalue_dict.items():
    color = color_map.get(gene_status_dict.get(gene, "Not Present"), "black")
    plt.scatter(gene, pvalue, color=color)

plt.xlabel('Genes')
plt.ylabel('Average P-value')
plt.title('Average P-values of Genes in ' + tissue)
plt.xticks(rotation=90) # Rotating gene names for readability
plt.subplots_adjust(bottom=0.4)
plt.savefig('average_pvals_colored.png')  # Save the plot with color coding

# Create a DataFrame for the boxplot
boxplot_data = pd.DataFrame(columns=['Gene Status', 'P-Value'])

# Populate the DataFrame
for gene, status in gene_status.items():
    if gene in gene_pvalue_dict:
        boxplot_data = boxplot_data.append({'Gene Status': status, 'P-Value': gene_pvalue_dict[gene]}, ignore_index=True)

# Plotting the boxplot
plt.figure(figsize=(10, 6))
sns.boxplot(x='Gene Status', y='P-Value', data=boxplot_data, palette=color_map)
plt.title('P-value Distribution Across Gene Status Groups in ' + tissue)
plt.savefig('pvalue_boxplot.png')  # Saves the boxplot

# Update the gene status CSV with average p-values
write_to_csv(gene_status, gene_pvalue_dict)
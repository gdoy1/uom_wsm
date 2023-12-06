import requests
import csv
import pygtex
import json
import matplotlib.pyplot as plt
import numpy as np


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

def write_to_csv(gene_status, filename='updated_gene_status.csv'):
    with open(filename, 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['Gene', 'Status'])
        for gene, status in gene_status.items():
            writer.writerow([gene, status])
    print(f"Updated gene status saved to {filename}")

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

# Example usage
panel_id = 'R59'
input_filename = 'most_central_community_genes.csv'
status_filename = 'updated_gene_status.csv'
tissue = 'Brain_Frontal_Cortex_BA9'

most_central_genes = read_gene_list(input_filename)
panel_genes = get_genes_from_panel(panel_id)
gene_status = compare_genes_with_panel(panel_genes, most_central_genes)
write_to_csv(gene_status)

gene_pvalue_dict = {}
gene_status_dict = read_gene_status(status_filename)

color_map = {
    "3": "green",
    "2": "#FFBF00",  # Hex code for amber
    "1": "red",
    "Not Present": "black"
}

for gene in most_central_genes:
    gencode_id = get_gencode_id_for_gene(gene)
    if gencode_id:
        eqtl_data = get_eqtl_data(gencode_id, tissue)
        average_pvalue = calculate_average_pvalue(eqtl_data)
        if average_pvalue is not None:
            gene_pvalue_dict[gene] = average_pvalue

# Plotting
for gene, pvalue in gene_pvalue_dict.items():
    color = color_map.get(gene_status_dict.get(gene, "Not Present"), "black")
    plt.scatter(gene, pvalue, color=color)

plt.xlabel('Genes')
plt.ylabel('Average P-value')
plt.title('Average P-values of Genes in ' + tissue)
plt.xticks(rotation=90) # Rotating gene names for readability
plt.savefig('average_pvals_colored.png')  # Save the plot with color coding

# Calculate average p-values for each status group
group_pvalues = {"1": [], "2": [], "3": [], "Not Present": []}
for gene, status in gene_status.items():
    if gene in gene_pvalue_dict:
        group_pvalues[status].append(gene_pvalue_dict[gene])

# Compute average p-values for each group
avg_pvalues = {group: np.mean(pvalues) if pvalues else None for group, pvalues in group_pvalues.items()}

# Prepare data for bar chart
groups = ['1', '2', '3', 'Not Present']
avg_pvalues_list = [avg_pvalues[group] for group in groups]

# Plotting bar chart
plt.bar(groups, avg_pvalues_list, color=['red', '#FFBF00', 'green', 'black'])
plt.xlabel('Gene Status')
plt.ylabel('Average P-value')
plt.title('Average P-values by Gene Status')
plt.savefig('average_pvals_by_gene_status.png')
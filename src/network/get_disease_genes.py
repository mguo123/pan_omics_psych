# Margaret GUo
# Get the disease genes from disease_to_genes.txt


### EXTRACT out cancer related genes

genes_info_file = ''../data/Homo_sapiens.gene_info'
'

file_path = '../data/disease_to_genes.txt'

with open(file_path, 'r') as f:
	lines = f.readlines()

for line in lines:
	if "cancer" in line:
		genes_str = line.split("\"")
		print(genes_str)

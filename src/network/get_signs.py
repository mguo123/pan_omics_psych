# Margaret
# get signs and directions of the Pathways Commmon 9
import numpy as np
import os
import networkx as nx
#####################################################################################################################
################################################## Create dictionaries ######################################################
#####################################################################################################################
corr_dir = 'HSA-coexpression/'
for file in os.listdir(corr_dir):



#####################################################################################################################
################################################## LOAD FILES ######################################################
#####################################################################################################################

# data_arr = [mx3] were m is the number of edges, the 
# first column is the starting node (as a gene symbol), the second column is relationship ontology, the third column is the ending node (as a gene symbol), 
# 1546602 edges
file_name = 'PathwayCommons9.All.hgnc.sif'
with open(file_name, 'r') as f:
    x = f.readlines()

    data_arr = [line.strip().split('\t') for line in x][:10]

# get relatinship mapping dictionary of the relationship ontology from Pathway commons to the type of edge that exists
relationship_csv = 'pathwaycommons_relationship_types.csv'
with open(relationship_csv, 'r') as f:
    relationship_dict = {}
    for relationship in f.readlines():
        name, description, type_r = relationship.strip().split(',')
        # type_r == 'N': undirected, type_r == 'AB' directed at to b
        relationship_dict[name] = type_r


# get the mapping between the gene_id (entrez Id) and the gene symbol (the symbol in Pathways graph) 
id_mapping = 'Homo_sapiens.gene_info.csv'
with open(id_mapping, 'r') as g:
    # columns = tax_id GeneID  Symbol  LocusTag    Synonyms    dbXrefs chromosome  map_location    description type_of_gene   ...
    # Symbol_from_nomenclature_authority  Full_name_from_nomenclature_authority   Nomenclature_status Other_designations  Modification_date   Feature_type
    id_to_sym = {}
    sym_to_id = {}
    for idx, line in enumerate(g.readlines()):
        if idx > 0:
            line = line.strip().split(',')
            id_to_sym[line[1]] = line[2]
            sym_to_id[line[2]] = line[1]

#####################################################################################################################
################################################## Find Correlation ######################################################
#####################################################################################################################
# for every edge, see if the coexpression database has an correlation value, report that correlation

def get_correl(gene_B_id, file_corr):
    with open(file_corr, 'r') as f:
        for line in f:
            if gene_B_id in line:
                 corr = float(line.strip().split('\t')[-1])
                 return corr
    return 0

corr = np.zeros(len(data_arr))
rev_corr = np.zeros(len(data_arr))
corr_dir = 'HSA-coexpression/'
new_edges = []
for idx, line in enumerate(data_arr):
    gene_A, relationship, gene_B = line

    # find correlation
    corr = 0
    if gene_A in sym_to_id and gene_B in sym_to_id:
        gene_A_id = sym_to_id[gene_A]
        gene_B_id = sym_to_id[gene_B]
        file_corr = os.path.join(corr_dir, gene_A_id)
        corr[idx] = get_correl(gene_B_id, file_corr)
        rev_corr[idx] = get_correl(gene_A_id, file_corr)

    ## Write new file with graph information
    edge_type = 'N'
    if relationship in relationship_dict:
        edge_type = relationship_dict[relationship]
    new_edges.append((gene_A, gene_B, corr[idx]))
    if edge_type == 'N': # undirected graph
        new_edges.append((gene_B, gene_A,  rev_corr[idx])) # if undirected add the correlation in the opposite direction

print(new_edges)
# create both the multi edge graph and the undirected single edge graph, need a single edge graph to run some forms of analysis
MG=nx.MultiGraph()
MG.add_weighted_edges_from([new_edges])

GG=nx.Graph()
for n,nbrs in MG.adjacency_iter():
    for nbr,edict in nbrs.items():
        minvalue=min([d['weight'] for d in edict.values()]) 
        GG.add_edge(n,nbr, weight = minvalue)
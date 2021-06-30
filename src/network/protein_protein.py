
# coding: utf-8

# In[1]:

import csv
from collections import defaultdict
#import matplotlib.pyplot as plt
import numpy as np
import scipy
import snap
from math import log, exp
from scipy.special import gammaln
from scipy.stats import combine_pvalues
from operator import itemgetter

def AddNode(Graph, node_to_id, node_name, new_node_id):
    # This function adds a new node to the graph and
    # updates the protein name -> node id mapping.
    Graph.AddNode(new_node_id)
    node_to_id[node_name] = new_node_id
    return new_node_id + 1

def OpenCorrelationInfoFile(path):
    # Opens the edge list file at 'path',
    # construct the graph and the node name to node id mapping,
    # and return both.
    with open(path, 'r') as f:
        # Skip the first line
        next(f)
        # Read all lines, constructing the graph as we go.
        new_node_id = 1
        Graph = snap.TNEANet.New()
        node_to_id = {}
        for line in f:
            values = line.split()
            start_node_id = -1
            end_node_id = -1
            if values[0] not in node_to_id:
                AddNode(Graph, node_to_id, values[0], new_node_id)
                start_node_id = new_node_id
                new_node_id += 1
            else:
                start_node_id = node_to_id[values[0]]
            if values[1] not in node_to_id:
                AddNode(Graph, node_to_id, values[1], new_node_id)
                end_node_id = new_node_id
                new_node_id += 1
            else:
                end_node_id = node_to_id[values[1]]
            
            # Add the edge
            Graph.AddEdge(start_node_id, end_node_id)

        return Graph, node_to_id

def OpenDiseaseGeneCountFile(path, node_to_id):
    with open(path, 'r') as f:
        next(f)
        disease_genes = set()
        total_line_count = 0
        total_disease_genes_count = 0
        for line in f:
            total_line_count += 1
            values = line.split()
            if values[0] in node_to_id:
                total_disease_genes_count += 1
                disease_genes.add(node_to_id[values[0]])
        print 'Total lines=%d\tDisease Genes=%d' % (total_line_count, total_disease_genes_count)
        return disease_genes
    
def ComputeInDegreeDistribution(Graph):
    # Compute the in degree distribution and return X, Y values for plotting
    d = defaultdict(lambda: 0)
    for node in Graph.Nodes():
        d[node.GetInDeg()] += 1
    
    d = {k: float(d[k]) / Graph.GetNodes() for k in d.keys()}
    X = d.keys()
    Y = d.values()
    return X, Y

def ComputeOutDegreeDistribution(Graph):
    # Compute the in degree distribution and return X, Y values for plotting
    d = defaultdict(lambda: 0)
    for node in Graph.Nodes():
        d[node.GetOutDeg()] += 1
    
    d = {k: float(d[k]) / Graph.GetNodes() for k in d.keys()}
    X = d.keys()
    Y = d.values()
    return X, Y

def ComputeStatistics(Graph):
    # Compute some simple statistics. 
    # Node counts, edge counts, clustering coeffs, degree distributions.
    print 'Number of nodes = %d' % Graph.GetNodes()
    print 'Number of edges = %d' % Graph.GetEdges()
    print 'Average clustering coefficient = %f' % snap.GetClustCf(Graph)
    # Plot the in-deg and out-deg distributions, and fit a power law distribution to them.
    X_in, Y_in = ComputeInDegreeDistribution(Graph)
    X_out, Y_out = ComputeOutDegreeDistribution(Graph)
    plt.loglog(X_in, Y_in, 'o', color = 'r', label = 'In Degree')
    plt.loglog(X_out, Y_out, '+', color = 'b', label = 'Out Degree')
    plt.xlabel('degree k (log)')
    plt.ylabel('Degree distribution P(k) (log)')
    plt.title('In/Out Degree distribution of Protein-Protein Network')
    plt.legend()
    plt.subplots_adjust(top=0.85)
    plt.savefig('degree_dist.png', dpi = 300)
    plt.show()
    # TODO: We should also compute some other statistics, such as largest CC, largest weakly connected CC,
    # local modularity, fraction of known disease proteins, etc.


# ComputeStatistics(Graph)

# In this section I'll start implementing the DIAMOND algorithm. We start with the definition to the link
# significance as defined by the DIAMOND paper. The link significance is defined as the probability that,
# if we randomly selected s_0 seed proteins as disease proteins, what is the probability that a node with
# k out edges point to exactly k_s of those seed proteins? This is a hypergeometric distribution, where
# the total event space of the node with k out edges is of size N choose k, and the situations in which exactly
# k_s of those k edges point to the s_0 seed proteins is (s_0 choose k_s) * ((N - s_0) choose (k - k_s)).
# Note that this is not exactly drawing black/white balls out of bins, so we cannot directly use scipy.hypergeom.
# This function needs to be implemented with some care, namely:
#   * We are dealing with large numbers, so directly applying N choose k is bad. Instead we should operate in
#     log space.
#   * k! can be expressed as the Gamma function Gamma(k + 1)
#   * Thus, N choose k = Gamma(N + 1) / (Gamma(k + 1) * Gamma(N - k + 1))
#   * In log space, log( N choose k ) = LogGamma(N + 1) - (LogGamma(k + 1) + LogGamma(N - k + 1))
#   * scipy has a gammaln function for LogGamma.

def LogChoose(n, k):
    return gammaln(n + 1) - (gammaln(k + 1) + gammaln(n - k + 1))

# The probability mass function is defined as:
#   pmf(k_s, k, N, s_0) = (s_0 choose k_s) * ((N - s_0) choose (k - k_s)) / (N choose k)
# In log space:
#   ln pmf = ln(s_0 choose k_s) + ln((N - s_0) choose (k - k_s)) - ln(N choose k)
#          = LogChoose(s_0, k_s) + LogChoose(N - s_0, k - k_s) - LogChoose(N, k)
# We can exponentiate that to get the answer.
def DiseaseProteinHyperGeomPmf(k, k_s, N, s_0):
    log_pmf = LogChoose(s_0, k_s) + LogChoose(N - s_0, k - k_s) - LogChoose(N, k)
    return exp(log_pmf)

# Finally, we can compute the link signficance as the p-value of the null hypothesis that the linking
# is random. This is just the complementary CDF of the hypergeometric distribution:
#   p-value = P( X >= k_s ) = \sum_{k_i = k_s}^{k} pmf(k_i, k, N, s_0)

def CCDF_HyperGeom(k, k_s, N, s_0):
    ccdf = 0.0
    for k_i in range(k_s, k + 1):
        ccdf += DiseaseProteinHyperGeomPmf(k, k_i, N, s_0)
    return ccdf

def ComputeNumberOfLinksToSeedProteins(Graph, node_id, seed_node_ids):
    k_s = 0
    seed_neighbors = set()
    node = Graph.GetNI(node_id)
    
    # For the basic DIAMOND case, we will treat edges as undirected, and
    # sum up undirected edges to all neighboring seeded disease modules.
    for i in range(node.GetInDeg()):
        in_id = node.GetInNId(i)
        if in_id in seed_node_ids and in_id not in seed_neighbors:
            seed_neighbors.add(in_id)
            k_s += 1
    for i in range(node.GetOutDeg()):
        out_id = node.GetOutNId(i)
        if out_id in seed_node_ids and out_id not in seed_neighbors:
            seed_neighbors.add(out_id)
            k_s += 1
    return k_s

# This function is needed because we can't simply call node.GetDeg(), because for PNEANET, GetDeg() is 
# equal to the sum of the in-deg and out-deg, which is not necessarily the number of neighbors.
def ComputeNeighborCount(Graph, node_id):
    neighbors = set()
    node = Graph.GetNI(node_id)
    for i in range(node.GetDeg()):
        nid = node.GetNbrNId(i)
        if nid not in neighbors:
            neighbors.add(nid)
    return len(neighbors)
    
def ComputePValue(node_id, neighbor_count, k_s_map, N, s_0):
    k = neighbor_count[node_id]
    k_s = k_s_map[node_id]
    return CCDF_HyperGeom(k, k_s, N, s_0)

def ComputePValueMap(N, s_0, neighbor_count, k_s_map):
    p_value_map = {}
    i = 0
    for node_id in neighbor_count.keys():
        if i % 100 == 0:
            print 'Compute PValueMap: iteration %d' % i
        p_value_map[node_id] = ComputePValue(node_id, neighbor_count, k_s_map, N, s_0)
        i+= 1
    return p_value_map

def DIAMOND(Graph, seed_node_ids, n_iterations=100):
    seed_node_ids_prime = seed_node_ids.copy()
    # Compute a mapping from node to neighbors count.
    neighbor_count = {node.GetId(): ComputeNeighborCount(Graph, node.GetId())
                      for node in Graph.Nodes()}
    # Compute a mapping from each node to the number of seed proteins it is connected to.
    k_s_map = {node.GetId(): ComputeNumberOfLinksToSeedProteins(Graph, node.GetId(), seed_node_ids_prime)
               for node in Graph.Nodes()}

    N = Graph.GetNodes()
    s_0 = len(seed_node_ids)
    
    for i in range(n_iterations):
        p_value_map = ComputePValueMap(N, s_0, neighbor_count, k_s_map)
        p_values = [(k, v) for k, v in p_value_map.iteritems()]
        sorted_p_values = sorted(p_values, key=itemgetter(1))
        # Get the lowest p-value node.
        lowest_id = sorted_p_values[0][0]
        p_value = sorted_p_values[0][1]
        print 'Added node: %d\tp-value=%.10f' % (lowest_id, p_value)
        print 'k=%d\tk_s=%d\tN=%d\ts_0=%d' % (neighbor_count[lowest_id], k_s_map[lowest_id], N, s_0)
        # Add to the seed nodes.
        seed_node_ids_prime.add(lowest_id)
        # Remove from neighbor_count and k_s_map.
        del neighbor_count[lowest_id]
        del k_s_map[lowest_id]
        # Update the k_s_map. It only changes for the neighbors of lowest_id.
        neighbors = set()
        lowest_node = Graph.GetNI(lowest_id)
        for i in range(lowest_node.GetDeg()):
            nid = lowest_node.GetNbrNId(i)
            if nid not in neighbors:
                neighbors.add(nid)
        for nid in neighbors:
            k_s_map[nid] += 1


# How do we modify DIAMOND to work for signed edges? Essentially, we want to calculate two different values:
#   -- For positive edges of a node, given that it has k+ positive edges, what's the probability that k_s+ of them
#      are linked to seed nodes? It would assume the form:
#        (s_0 choose k_s+) * ((N - s_0) choose (k+ - k_s+)) / (N choose k+)
#   -- For negative edges, by symmetry we have a similar form.

def OpenThresholdedFile(path):
    # Opens the edge list file at 'path',
    # construct the graph and the node name to node id mapping,
    # and return both.
    with open(path, 'r') as f:
        # Skip the first line
        next(f)
        # Read all lines, constructing the graph as we go.
        new_node_id = 1
        Graph = snap.TNEANet.New()
        node_to_id = {}
        edge_to_sign = {}
        for line in f:
            values = line.split()
            start_node_id = -1
            end_node_id = -1
            if values[0] not in node_to_id:
                AddNode(Graph, node_to_id, values[0], new_node_id)
                start_node_id = new_node_id
                new_node_id += 1
            else:
                start_node_id = node_to_id[values[0]]
            if values[1] not in node_to_id:
                AddNode(Graph, node_to_id, values[1], new_node_id)
                end_node_id = new_node_id
                new_node_id += 1
            else:
                end_node_id = node_to_id[values[1]]
            
            # Add the edge
            Graph.AddEdge(start_node_id, end_node_id)
            # Add the edge sign
            edge_to_sign[(start_node_id, end_node_id)] = int(values[2])

        return Graph, node_to_id, edge_to_sign

def OpenSeededGenesFile(path, node_to_id):
    with open(path, 'r') as f:
        disease_genes = set()
        total_line_count = 0
        total_disease_genes_count = 0
        for line in f:
            total_line_count += 1
            if line in node_to_id:
                total_disease_genes_count += 1
                disease_genes.add(node_to_id[line])
        print 'Total lines=%d\tDisease Genes=%d' % (total_line_count, total_disease_genes_count)
        return disease_genes

def ComputeNumberOfPosNegLinksToSeedProteins(Graph, node_id, seed_node_ids, edge_to_sign, pos_or_neg):
    # This function computes k_s+ / k_s-, aka the total number of positive / negative links pointing 
    # to seed nodes, per node.
    k_s = 0
    node = Graph.GetNI(node_id)
    
    # Count the number of out edges to seed proteins that have the same sign as pos_or_neg.
    for i in range(node.GetOutDeg()):
        out_id = node.GetOutNId(i)
        if out_id in seed_node_ids and edge_to_sign[(node_id, out_id)] == pos_or_neg:
            k_s += 1
    return k_s

def ComputePosNegNeighborCount(Graph, node_id, edge_to_sign, pos_or_neg):
    # This function computes k+ / k-, aka the total number of positive / negative links per node.
    neighbor_count = 0
    node = Graph.GetNI(node_id)
    for i in range(node.GetOutDeg()):
        out_id = node.GetOutNId(i)
        if edge_to_sign[(node_id, out_id)] == pos_or_neg:
            neighbor_count += 1
    return neighbor_count

def ComputePosNegPValue(node_id, pos_neighbor_count, neg_neighbor_count, pos_k_s_map, neg_k_s_map, N, s_0):
    # If a particular node does not have any positive edges, then the CCDF_HyperGeom will come out to be 0.
    # Compute p_value+
    if node_id in pos_neighbor_count:
        k_pos = pos_neighbor_count[node_id]
    else:
        k_pos = 0
    if node_id in pos_k_s_map:
        k_s_pos = pos_k_s_map[node_id]
    else:
        k_s_pos = 0
    pos_p_value = CCDF_HyperGeom(k_pos, k_s_pos, N, s_0)
    # Compute p_value-
    if node_id in neg_neighbor_count:
        k_neg = neg_neighbor_count[node_id]
    else:
        k_neg = 0
    if node_id in neg_k_s_map:
        k_s_neg = neg_k_s_map[node_id]
    else:
        k_s_neg = 0
    neg_p_value = CCDF_HyperGeom(k_neg, k_s_neg, N, s_0)
    # Compute the chi squared statistic (fisher's method): \chi^2 = -2 * (ln p_value+ + ln p_value-).
    # The combined p-value is returned from scipy.stats.combine_pvalues[1]
    return combine_pvalues([pos_p_value, neg_p_value], method='fisher', weights=None)[1]

def ComputePosNegPValueMap(Graph, N, s_0, pos_neighbor_count, neg_neighbor_count, pos_k_s_map, neg_k_s_map, seed_node_ids):
    p_value_map = {}
    i = 0
    for node in Graph.Nodes():
        node_id = node.GetId()
        if node_id in seed_node_ids:
            continue
        #if i % 100 == 0:
        #    print 'Compute PValueMap: iteration %d' % i
        p_value_map[node_id] = ComputePosNegPValue(node_id, pos_neighbor_count, neg_neighbor_count, pos_k_s_map, neg_k_s_map, N, s_0)
        i+= 1
    return p_value_map

def PosNegDIAMOND(Graph, seed_node_ids, edge_to_sign, p_value_cutoff=0.05, n_cutoff=500):
    seed_node_ids_prime = seed_node_ids.copy()
    # Compute a mapping from node to neighbors count.
    pos_neighbor_count = {node.GetId(): ComputePosNegNeighborCount(Graph, node.GetId(), edge_to_sign, 1)
                          for node in Graph.Nodes() if node not in seed_node_ids_prime}
    print 'Total number of positive edges: %d' % np.sum([v for k, v in pos_neighbor_count.iteritems()])
    neg_neighbor_count = {node.GetId(): ComputePosNegNeighborCount(Graph, node.GetId(), edge_to_sign, 0)
                          for node in Graph.Nodes() if node not in seed_node_ids_prime}
    print 'Total number of negative edges: %d' % np.sum([v for k, v in neg_neighbor_count.iteritems()])
    # Compute a mapping from each node to the number of seed proteins it is connected to.
    pos_k_s_map = {node.GetId(): ComputeNumberOfPosNegLinksToSeedProteins(Graph, node.GetId(), seed_node_ids_prime, edge_to_sign, 1)
                   for node in Graph.Nodes()}
    neg_k_s_map = {node.GetId(): ComputeNumberOfPosNegLinksToSeedProteins(Graph, node.GetId(), seed_node_ids_prime, edge_to_sign, 0)
                   for node in Graph.Nodes()}

    # Sanity Check
    for k, v in pos_neighbor_count.iteritems():
        assert v >= pos_k_s_map[k]
    for k, v in neg_neighbor_count.iteritems():
        assert v >= neg_k_s_map[k]
                   
    N = Graph.GetNodes()
    
    last_pv = 0
    it = 0
    newly_found_disease_genes = set()
    while last_pv < p_value_cutoff and len(newly_found_disease_genes) < n_cutoff:
        s_0 = len(seed_node_ids_prime)
        print 'Running iteration %d' % it
        p_value_map = ComputePosNegPValueMap(Graph, N, s_0, pos_neighbor_count, neg_neighbor_count, pos_k_s_map, neg_k_s_map, seed_node_ids_prime)
        p_values = [(k, v) for k, v in p_value_map.iteritems()]
        sorted_p_values = sorted(p_values, key=itemgetter(1))
        # Get the lowest p-value node.
        lowest_id = sorted_p_values[0][0]
        p_value = sorted_p_values[0][1]
        print 'Added node: %d\tp-value=%.10f' % (lowest_id, p_value)
        print 'k+=%d\tk_s+=%d\tk-=%d\tk_s-=%d\tN=%d\ts_0=%d' % (pos_neighbor_count[lowest_id], pos_k_s_map[lowest_id], neg_neighbor_count[lowest_id], neg_k_s_map[lowest_id], N, s_0)
        # Add to the seed nodes.
        seed_node_ids_prime.add(lowest_id)
        newly_found_disease_genes.add(lowest_id)
        # Remove from neighbor_count and k_s_map.
        del pos_neighbor_count[lowest_id]
        del neg_neighbor_count[lowest_id]
        del pos_k_s_map[lowest_id]
        del neg_k_s_map[lowest_id]
        # Update the k_s_map. It only changes for the neighbors of lowest_id.
        neighbors = set()
        lowest_node = Graph.GetNI(lowest_id)
        for i in range(lowest_node.GetInDeg()):
            nid = lowest_node.GetInNId(i)
            if nid not in seed_node_ids_prime:
                neighbors.add(nid)
        for nid in neighbors:
            #print 'nid: %d k+=%d k-=%d k_s+=%d k_s-=%d' % (nid, pos_neighbor_count[nid], neg_neighbor_count[nid], pos_k_s_map[nid], neg_k_s_map[nid])
            if nid in pos_k_s_map and (nid, lowest_id) in edge_to_sign and edge_to_sign[(nid, lowest_id)] == 1:
                pos_k_s_map[nid] += 1
            elif nid in neg_k_s_map and (nid, lowest_id) in edge_to_sign and edge_to_sign[(nid, lowest_id)] == 0:
                neg_k_s_map[nid] += 1
            assert pos_neighbor_count[nid] >= pos_k_s_map[nid]
            assert neg_neighbor_count[nid] >= neg_k_s_map[nid]
        last_pv = p_value
        it += 1

    return seed_node_ids_prime, newly_found_disease_genes

def WriteNewDiseasesToFile(new_disease_genes, node_to_id, path):
    # First compute the reverse mapping of node_to_id
    id_to_node = {v: k for k, v in node_to_id.iteritems()}
    new_disease_gene_names = [id_to_node[i] for i in new_disease_genes]
    with open(path, 'w') as f:
        for item in new_disease_gene_names:
            f.write("%s\n" % item)
            
def OpenGenesCsvFile(path, node_to_id, col_idx):
    disease_genes = set()
    with open(path, 'r') as f:
        next(f)
        next(f)
        i = 0
        for line in f:
            cols = line.split(',')
            gene = cols[col_idx]
            if gene in node_to_id:
                disease_genes.add(node_to_id[gene])
            i+= 1
        print 'Total number of disease genes in column: %d\tNumber of disease genes found in network: %d' % (i, len(disease_genes))
        return disease_genes
            

Graph, node_to_id, edge_to_sign = OpenThresholdedFile('new_edges_t0.5_binary')

#disease_genes = OpenGenesCsvFile('seed_nodes.csv', node_to_id, 0)
disease_genes = OpenDiseaseGeneCountFile('diease_gene-count.txt', node_to_id)
# disease_genes = OpenSeededGenesFile('cancer_genes.txt', node_to_id)

disease_module, new_disease_genes = PosNegDIAMOND(Graph, disease_genes, edge_to_sign, p_value_cutoff=0.05, n_cutoff=500)
WriteNewDiseasesToFile(new_disease_genes, node_to_id, 'new_disease_names.txt')


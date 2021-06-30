# Run random validation
from validation import read_results_file, get_disease_stats, get_go_stats, write_stats_file
from collections import defaultdict
import random

NUM_ITER = 10
NUM_RANDOM = 500

node_file = '../DATA/nodes.txt'
with open(node_file, 'r') as f:
	lines = f.readlines()
	nodes = [line.strip() for line in lines]
print(len(nodes), 'number of nodes')

# storage dictionaries
disease_overall = defaultdict(list)
go_overall = defaultdict(list)

for i in range(NUM_ITER):
    print(i)

    # generate random nodes
    random_nodes = set()
    while len(random_nodes) < NUM_RANDOM:
        random_nodes.add(random.choice(nodes))
    print(len(random_nodes), 'random nodes generated')


    # Get z-scores for diseases and  frequencies for go terms
    disease_score_avg, disease_score_str = get_disease_stats(random_nodes)
    go_freq = get_go_stats(random_nodes)
    for disease in disease_score_avg:
        disease_overall[disease].append(disease_score_avg[disease])
    for go in go_freq:
        go_overall[go].append(go_freq[go])

# find average overall all iterations for disease z-scores and go frequencies
disease_overall_avg = {}
go_overall_avg = {}
for disease in disease_overall:
    disease_overall_avg[disease ] = sum(disease_overall[disease])/float(NUM_ITER)
for go in go_overall:
    go_overall_avg[go] = sum(go_overall[go])/float(NUM_ITER)


out_name_diseases_avg = "../Results/random_disease_avg.txt"
out_name_go = "../Results/random_go.txt"


write_stats_file(disease_overall_avg, out_name_diseases_avg)
write_stats_file(go_overall_avg, out_name_go)

import requests
import nltk
import os,sys
from collections import defaultdict
import time
import pickle, csv

ACCEPTABLE_ENDINGS = ['v1.txt', 'v2.txt', 'up.txt', 'dn.txt']
REF_DIR = 'reference'
if not os.path.isdir(REF_DIR):
    os.mkdir(REF_DIR)

def query(data_type, uid):
    """
    :param data_type: the Pharos data type (e.g. ligand, target ..)
    :param uid: Pharos concept identifier
    :return: JSON containing the query results
    """
    url = "https://pharos.nih.gov/idg/api/v1/{}({})?view=full".format(data_type, uid)
    try:
        resp = requests.get(url).json()
    except:
        resp = None
    return resp


def get_target_id(target_common_name):
    """
    :param target_common_name: target name (NCBI gene name)
    :return: corresponding Pharos ID for the target
    """
    try:
        resp = query("targets", target_common_name)
        target_id = resp['id']
    except:
        target_id = None

    return target_id



def target_prop_data(qID, label):
    """
    :param qID: Pharos ID for the gene target
    :param label: property of interest
    :return: info - JSON of the query results
    """

    url = "https://pharos.nih.gov/idg/api/v1/targets/{}/properties(label={})".format(qID, label)
    info = requests.get(url).json()
    return info

def target_link_data(qID, kind):
    """
    :param qID: Pharos ID for the gene target
    :param kind: link of interest
    :return:
    """
    url = "https://pharos.nih.gov/idg/api/v1/targets/{}/links(kind={})".format(qID, kind)
    info = requests.get(url).json()
    return info


def get_target_diseases(target_common_name):
    qID = get_target_id(target_common_name)
    if qID is not None:
        info = target_link_data(qID, 'ix.idg.models.Disease')
        target_disease_scores = {}
        for association in info:
            try:
                if isinstance(association, dict):
                    properties = association['properties']
                    info_dict = {}
                    for prop in properties:
                        if prop['label'] == 'IDG Disease':
                            info_dict['disease'] = prop['term']
                        elif prop['label'] == 'IDG Confidence':
                            info_dict['confidence'] = float(prop['numval'])
                        elif prop['label'] == 'IDG Z-score':
                            info_dict['z_score'] = float(prop['numval'])
                    if 'z_score' in info_dict: # if we have some quantifiable measure for the data
                        target_disease_scores[info_dict['disease']] = info_dict['z_score']
            except:
                print(association)
                raise
        # for key, value in sorted(target_disease_scores.items(), key=lambda x: x[1], reverse=True):
        #     print ("%s: %s" % (key, value))
        return target_disease_scores
    else:
        return None

def write_target_disease_score(target_disease_scores, target_common_name):
    out_file_name = os.path.join(REF_DIR, target_common_name + '_disease_scores.pkl')
    print('wrote', out_file_name)
    if not os.path.isfile(out_file_name):
        with open(out_file_name, 'wb') as f:
            pickle.dump(target_disease_scores, f)


# wrapper function that reads from file the object containing information specific for that gene from pharos, is that info is not found, pull from pharos
def read_target_disease_score(target_common_name):
    in_file_name = os.path.join(REF_DIR, target_common_name + '_disease_scores.pkl')
    if os.path.isfile(in_file_name):
        # print('loaded', in_file_name)
        pickle_in = open(in_file_name,"rb")
        target_disease_scores = pickle.load(pickle_in)
    else:
        target_disease_scores = get_target_diseases(target_common_name)
        if target_disease_scores is not None:
            write_target_disease_score(target_disease_scores, target_common_name)

    return target_disease_scores


def get_GO_terms(target_name):
    """
    :param target_name: NCBI gene name
    :return: list of GO process, function, and component (list of 3 dictionaries) terms for that drug
    """

    # Query for the Pharos ID of the Gene target
    qID = get_target_id(target_name)

    if qID is not None:

        # Get tissue information from the properties for that target
        GO_info = target_prop_data(qID, "*GO")

        # Grab GO term annotations for each gene
        GO_process = set()
        GO_function = set()
        GO_component = set()
        for prop in GO_info:
            if isinstance(prop, dict):
                if prop["label"] == 'GO Process':
                    GO_process.add(prop['term'])
                elif prop["label"] == 'GO Function':
                    GO_function.add(prop['term'])
                elif prop["label"] == 'GO Component':
                    GO_component.add(prop['term'])
        if len(GO_process) == 0:
            GO_process.add("unknown")
        if len(GO_function) == 0:
            GO_function.add("unknown")
        if len(GO_component) == 0:
            GO_component.add("unknown")

        return GO_process, GO_function, GO_component
    else:
        return None


def write_target_go_terms(GO_process, GO_function, GO_component, target_common_name):

    out_file_name = os.path.join(REF_DIR, target_common_name + '_go_terms.pkl')
    print('wrote', out_file_name)
    if not os.path.isfile(out_file_name):
        with open(out_file_name, 'wb') as f:
            pickle.dump([GO_process, GO_function, GO_component], f)


def read_target_go_terms(target_common_name):
    in_file_name = os.path.join(REF_DIR, target_common_name + '_go_terms.pkl')
    if os.path.isfile(in_file_name):
        # print('read', in_file_name)
        pickle_in = open(in_file_name,"rb")
        GO_process, GO_function, GO_component = pickle.load(pickle_in)
        return GO_process, GO_function, GO_component
    else:
        result = get_GO_terms(target_common_name)
        if result is not None:
            GO_process, GO_function, GO_component = result
            write_target_go_terms(GO_process, GO_function, GO_component, target_common_name)
            return GO_process, GO_function, GO_component
        else:
            return None
    


def read_results_file(filepath):
    print('read in', filepath)
    with open(filepath, 'r') as f:
        lines = f.readlines()

        genes = [line.strip() for line in lines]

    return genes # list of gene names

def get_disease_stats(gene_list):
    # store disease z-scores for all genes as a list of z-scores, these will be summed and the average over all of the possible genes will be found
    disease_scores_list = defaultdict(list)
    for i, gene in enumerate(gene_list):
        # print(i)        
        target_disease_scores = read_target_disease_score(gene)
        if target_disease_scores is not None:
            # print(gene, '!!!!!!!!!!!!!!!!')
            for key, value in target_disease_scores.items():
                # print(key, value)
                disease_scores_list[key].append(value)


    disease_score_str = {}

    disease_score_avg = {}
    for key, list_val in disease_scores_list.items():
        disease_score_avg[key] = sum(list_val)/float(len(gene_list))
        disease_score_str[key] = '\t'.join([str(s) for s in list_val])

    return disease_score_avg, disease_score_str


def get_go_stats(gene_list):
    num_gene = float(len(gene_list))
    go_freq = defaultdict(int)
    for i, gene in enumerate(gene_list):
        # print(i)
        GO_info = read_target_go_terms(gene)
        if GO_info is not None:
            GO_process = GO_info[0] # don't care about components or function
            for go_term in GO_process: #iterate through set
                go_freq[go_term] += 1/num_gene

    return go_freq


def write_stats_file(dict_to_write, out_file_name):
    print('wrote file', out_file_name)
    with open(out_file_name, 'w') as f:
        for key, value in sorted(dict_to_write.items(), key=lambda x: x[1], reverse=True):

            f.write ("%s\t%s\n" % (key, value))

def get_node_to_id(path):
    # Opens the edge list file at 'path',
    # construct the node name to node id mapping,
    # and return both.
    print('open threshold file:', path)
    with open(path, 'r') as f:
        # Skip the first line
        next(f)
        # Read all lines, constructing the graph as we go.
        new_node_id = 1
        node_to_id = {}
        edge_to_sign = {}
        for line in f:
            values = line.split()
            start_node_id = -1
            end_node_id = -1
            if values[0] not in node_to_id:
                node_to_id[values[0]] = new_node_id
                start_node_id = new_node_id
                new_node_id += 1
            else:
                start_node_id = node_to_id[values[0]]
            if values[1] not in node_to_id:
                node_to_id[values[1]] = new_node_id
                end_node_id = new_node_id
                new_node_id += 1
            else:
                end_node_id = node_to_id[values[1]]


        return node_to_id


def get_original_genes(path, node_to_id, col_idxs):
    #col_idxs list of integers referring to columns 
    print('open genes file', path, 'col', col_idxs)
    disease_genes = set()
    name_file_arr = []
    num_seed_genes = 0
    # get columns
    with open(path, 'rb') as f:
        lines = [row for row in csv.reader(f.read().splitlines())]

        # iterate through col indexes
        for col_idx in col_idxs:
            result=[]
            for x in lines:
                if x[col_idx] != '':
                    result.append(x[col_idx])
        
            name_file_arr.append(result[0].lower())
            num_seed_genes += int(result[1])
            i = 0
            for gene in result[2:]:
                if gene in node_to_id:
                    disease_genes.add(genes)
                i+= 1

    name_file = "_".join(name_file_arr)
    print('name_file:', name_file)
    print ('Total number of disease genes in columns: %d\tNumber of disease genes found in network: %d' % (num_seed_genes, len(disease_genes)))
    return disease_genes, name_file, num_seed_genes

def get_node_to_id(path):
    # Opens the edge list file at 'path',
    # construct the node name to node id mapping,
    # and return both.
    print('open threshold file:', path)
    with open(path, 'r') as f:
        # Skip the first line
        next(f)
        # Read all lines, constructing the graph as we go.
        new_node_id = 1
        node_to_id = {}
        edge_to_sign = {}
        for line in f:
            values = line.split()
            start_node_id = -1
            end_node_id = -1
            if values[0] not in node_to_id:
                node_to_id[values[0]] = new_node_id
                start_node_id = new_node_id
                new_node_id += 1
            else:
                start_node_id = node_to_id[values[0]]
            if values[1] not in node_to_id:
                node_to_id[values[1]] = new_node_id
                end_node_id = new_node_id
                new_node_id += 1
            else:
                end_node_id = node_to_id[values[1]]


        return node_to_id

def get_original_genes(path, node_to_id, col_idx):
    #col_idxs list of integers referring to columns 
    print('open genes file', path, 'col', col_idx)
    disease_genes = set()
    name_file_arr = []
    num_seed_genes = 0
    # get columns
    with open(path, 'r') as f:
        lines = [row for row in csv.reader(f.read().splitlines())]

        # iterate through col indexes
        name_file_arr.append(str(col_idx))

        result=[]
        for x in lines:
            if x[col_idx] != '':
                result.append(x[col_idx])
    
        name_file_arr.append(result[0].lower())
        num_seed_genes += int(result[1])
        i = 0
        for gene in result[2:]:
            if gene in node_to_id:
                disease_genes.add(gene)
            i+= 1

    name_file = "_".join(name_file_arr)
    print('name_file:', name_file)
    print ('Total number of disease genes in columns: %d\tNumber of disease genes found in network: %d' % (num_seed_genes, len(disease_genes)))
    return sorted(list(disease_genes)), name_file, num_seed_genes


if __name__ == '__main__':
    # read in the output file of protein_protein_cmd.py which contains top 500 significant genes, and creates a dictionary with all of the 
    assert len(sys.argv) == 3, "must call command as python validation_withseed.py ../DATA/new_edges_t0.5_binary ../DATA/seed_nodes.csv"

    result_path = '../Results'

    graph_path = sys.argv[1]
    node_to_id = get_node_to_id(graph_path)
    print(len(node_to_id))

    seed_path = sys.argv[2]


    files_to_run = sorted(os.listdir(result_path), reverse=True)
    for file in files_to_run:
        if any([file.endswith(end) for end in ACCEPTABLE_ENDINGS]):
            out_name = os.path.splitext(os.path.basename(file))[0] 

            number  = int(out_name.split('_')[1])



            #### DEBUG
            # gene_list = gene_list[:10]
            out_name_diseases_avg = "".join([os.path.join(result_path,out_name), "_diseases_avg", ".txt"])
            out_name_diseases_all = "".join([os.path.join(result_path,out_name), "_diseases_all", ".txt"])
            out_name_go = "".join([os.path.join(result_path,out_name), "_go", ".txt"])

        # if os.path.isfile(out_name_diseases_all):
            print(file)
            out_genes = read_results_file(os.path.join(result_path,file))

            seed_genes, name_file, num_seed_genes = get_original_genes(seed_path, node_to_id, number)

            gene_list = seed_genes + out_genes
            print(len(gene_list))

            disease_score_avg, disease_score_str = get_disease_stats(gene_list)
            # print(disease_score_avg)

            write_stats_file(disease_score_avg, out_name_diseases_avg)
            write_stats_file(disease_score_str, out_name_diseases_all)
           
            go_freq = get_go_stats(gene_list)
            # print(go_freq)
            write_stats_file(go_freq, out_name_go)


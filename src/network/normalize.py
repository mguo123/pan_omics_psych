# normalize
import os
import math

THRES_DIS = 2/math.sqrt(500) # 0.089, defined as the population wide z-score required to get within 2 z scores 
THRES_GO = THRES_DIS/2
def read_results(file_path):
    result_dict = {}
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            try:
                key, value = line.strip().split('\t')
                result_dict[key] = float(value)
            except ValueError:
                print(line, 'not read')
                raise
    return result_dict

def normalize(raw_dict, rand_dict, thres=THRES_DIS):
    norm_dict={}
    for key, value in raw_dict.items():
        if key in rand_dict and value > thres:
            norm_dict[key] = value/ rand_dict[key]
    return norm_dict

def write_results(dict_to_write, out_file_name):
    print('wrote file', out_file_name)
    with open(out_file_name, 'w') as f:
        for key, value in sorted(dict_to_write.items(), key=lambda x: x[1], reverse=True):
            f.write ("%s\t%s\n" % (key, value))



# # get random disease dictionary
# rand_disease_file = "../Results/all_disease_names_diseases_avg.txt"
# rand_go_file = "../Results/all_disease_names_go.txt"

rand_disease_file = "../Results/random_disease_avg.txt"
rand_go_file = "../Results/random_go.txt"
rand_disease_dict = read_results(rand_disease_file)
rand_go_dict = read_results(rand_go_file)


# get paths of raw files
result_disease_path = '../Results/results_seed'
result_go_path = result_disease_path#'../Results/go_results_no_seed'
files_to_run_disease = sorted(os.listdir(result_disease_path)) 
files_to_run_go = sorted(os.listdir(result_go_path))

for file in files_to_run_disease:
    if any([file.endswith('start_avg.txt'),file.endswith('diseases_avg.txt')]):
        try:
            out_name = os.path.splitext(os.path.basename(file))[0] 
            # print(out_name, '111')
            # comparison_file = os.path.join(result_disease_path, out_name[:-3] + "start_avg.txt")
            # print(comparison_file)
            # comparison_dict = read_results(comparison_file)

            out_name_norm = "".join([os.path.join(result_disease_path,out_name), "_norm", ".txt"])
            # print(out_name_norm)
            # if output file does not already exist
            if not os.path.isfile(out_name_norm):
                print(file)

                raw_dict = read_results(os.path.join(result_disease_path,file))
                norm_dict = normalize(raw_dict, rand_disease_dict)
                write_results(norm_dict, out_name_norm)
        except:
            print('failed for', out_name)
            raise

for file in files_to_run_go:
    if (file.endswith('go.txt')):
        out_name = os.path.splitext(os.path.basename(file))[0] 

        out_name_norm = "".join([os.path.join(result_go_path,out_name), "_norm", ".txt"])

        # if output file does not already exist
        if not os.path.isfile(out_name_norm):
            print(file)

            raw_dict = read_results(os.path.join(result_go_path,file))
            norm_dict = normalize(raw_dict, rand_go_dict, thres=THRES_GO)
            write_results(norm_dict, out_name_norm)


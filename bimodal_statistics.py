import argparse
import numpy as np
import gzip
import math

top_k_percent = 0.1


def round_float_for_print(in_float):
    in_float = 0 if math.isnan(in_float) else in_float
    return str(round(in_float, 3))


def process_blocks(block_file, beta_file, out_file):
    beta_array = np.fromfile(beta_file, dtype=np.uint8).reshape((-1, 2))
    new_lines = []
    num_appends = 0
    with gzip.open(block_file, 'r') as fin:
        for line in fin:
            line = line.rstrip().decode()
            tokens = line.split("\t")
            start_ind = int(tokens[3])
            end_ind = int(tokens[4])
            is_first = True
            prev = 0
            max_dif = 0
            prop_list = []

            for cur_ind in range(start_ind, end_ind):
                row = beta_array[cur_ind]
                methyl_prop = row[0] / row[1]
                if not math.isnan(methyl_prop):
                    prop_list.append(methyl_prop)
            prop_list = sorted(prop_list)
            for el in prop_list:
                if is_first:
                    prev = el
                    is_first = False
                else:
                    cur_dif = abs(el - prev)
                    prev = el
                    if cur_dif > max_dif:
                        max_dif = cur_dif
            if len(prop_list) >= 5:
                top_k_elements = int(top_k_percent * len(prop_list))
                top_k_elements = 1 if top_k_elements == 0 else top_k_elements
                bottom_k_val_avg = round_float_for_print(np.mean(prop_list[0:top_k_elements]))
                tail = len(prop_list) - top_k_elements
                top_k_val_avg = round_float_for_print(np.mean(prop_list[tail:]))
                max_dif = round_float_for_print(max_dif)
                new_line = line + "\t{}\t{}\t{}\n".format(max_dif, bottom_k_val_avg, top_k_val_avg)
            else:
                new_line = line + "\t{}\t{}\t{}\n".format("0", "0", "0")
            num_appends += 1
            new_lines.append(new_line)
            if len(new_lines) > 500:
                with open(out_file, "a") as myfile:
                    for n_line in new_lines:
                        myfile.write(n_line)
                new_lines = []
    print(num_appends)
    with open(out_file, "a") as myfile:
        for n_line in new_lines:
            myfile.write(n_line)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--block_file', default="/cs/zbio/jrosensk/block_files_2/homog/homog_global/Prostate-Epithelial-Z000000RV.homog.gz")
    parser.add_argument('--beta_file', default="/cs/cbio/jon/grail_atlas/data/Prostate-Epithelial-Z000000RV.beta")
    parser.add_argument('--out_file', default="/cs/zbio/jrosensk/block_files_2/homog/homog_global/Prostate-Epithelial-Z000000RV.statistics.homog")
    return parser.parse_args()


def main():
    args = parse_args()
    process_blocks(args.block_file, args.beta_file, args.out_file)


if __name__ == '__main__':
    main()
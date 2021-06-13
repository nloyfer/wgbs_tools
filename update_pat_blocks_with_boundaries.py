import argparse

def get_boundary(boundary_it):
    cur_start_boundary_line = boundary_it.readline()
    start_boundary_tokens = cur_start_boundary_line.split("\t")
    return int(start_boundary_tokens[1])


def find_new_boundary(boundary_it, cand_index, cur_position, prev_position):
    cur_start_bound_distance = abs(cand_index - cur_position)
    prev_start_bound_distance = abs(cand_index - prev_position)
    while cur_start_bound_distance < prev_start_bound_distance:
        prev_position = cur_position
        cur_position = get_boundary(boundary_it)
        cur_start_bound_distance = abs(cand_index - cur_position)
        prev_start_bound_distance = abs(cand_index - prev_position)
    return prev_position, cur_position


def block_boundary_iterator(block_file, boundary_file, out_file):
    with open(block_file, 'rt') as block_lines, open(boundary_file, 'rt') as start_boundaries, open(boundary_file, 'rt'
                                                                                                    ) as end_boundaries:
        prev_start_boundary = get_boundary(start_boundaries)
        cur_start_boundary = get_boundary(start_boundaries)
        prev_end_boundary = get_boundary(end_boundaries)
        cur_end_boundary = get_boundary(end_boundaries)
        has_next = True
        to_print_list = []
        counter = 0
        while has_next:
            try:
                cur_block = block_lines.readline()
                if cur_block:
                    counter += 1
                    block_tokens = cur_block.split("\t")
                    block_start = int(block_tokens[3])
                    block_end = int(block_tokens[4])
                    prev_start_boundary, cur_start_boundary = find_new_boundary(start_boundaries, block_start,
                                                                                cur_start_boundary, prev_start_boundary)
                    prev_end_boundary, cur_end_boundary = find_new_boundary(end_boundaries, block_end,
                                                                            cur_end_boundary, prev_end_boundary)
                    # if prev_end_boundary == prev_start_boundary:
                    #     actual_end_boundary += prev_end_boundary + 2
                    to_print = "{}\t{}\t{}\n".format(block_tokens[0], prev_start_boundary, prev_end_boundary)
                    to_print_list.append(to_print)
                    if len(to_print_list) > 1000:
                        with open(out_file, "a") as myfile:
                            for to_print in to_print_list:
                                myfile.write(to_print)
                        to_print_list = []
                else:
                    has_next = False
            except StopIteration as e:
                has_next = False
        with open(out_file, "a") as myfile:
            for to_print in to_print_list:
                myfile.write(to_print)
        print(counter)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--block_file')
    parser.add_argument('--boundary_file')
    parser.add_argument('--out_file')
    return parser.parse_args()


def main():
    args = parse_args()
    block_boundary_iterator(args.block_file, args.boundary_file, args.out_file)


if __name__ == '__main__':
    main()
    # block_boundary_iterator("/cs/zbio/jrosensk/block_files_2/Aorta-Endothel-Z00000422.blocks.tsv",
    #                         "/cs/cbio/jon/segmentation_files/block_boundaries.tsv",
    #                         "/cs/zbio/jrosensk/block_files_2/Aorta-Endothel-Z00000422.blocks.unified.tsv")

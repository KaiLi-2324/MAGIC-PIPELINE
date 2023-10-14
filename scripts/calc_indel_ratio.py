import os
import csv
import argparse
from collections import defaultdict


def read_in_indel(path_indel_vcf):
    print(f"Start reading from {path_indel_vcf}")
    indel_pos = defaultdict(lambda: "")
    with open(path_indel_vcf) as file:
        for line in filter(lambda k: not k.startswith("#"), file):
            items = line.strip().split("\t")
            pos = f"{items[0]}-{items[1]}"
            if len(items[3]) > len(items[4]):
                indel_pos[pos] = "del"
            elif len(items[3]) < len(items[4]):
                indel_pos[pos] = "in"
    return indel_pos


def read_in_gt(path_gt, path_indel_vcf):
    indel_pos = read_in_indel(path_indel_vcf=path_indel_vcf)
    print(f"Start reading from {path_gt}")
    sample_indel_counts = defaultdict(lambda: defaultdict(lambda: 0))
    for item in csv.DictReader(open(path_gt), delimiter="\t"):
        pos = f"{item['CHROM']}-{item['POS']}"
        indel_type = indel_pos[pos]
        if not indel_type:
            continue
        samples_with_count = dict(filter(lambda k: k[1] == "0/1" or k[1] == "1/1", item.items()))
        for each_sample in samples_with_count:
            sample_indel_counts[each_sample][indel_type] += 1
    return sample_indel_counts


def write_indel_ratio(sample_indel_counts, path_out):
    print(f"Start writing to {path_out}")
    with open(path_out, "w+") as out:
        out.write("sample\tin/del\n")
        for each_sample in sample_indel_counts:
            insertion_counts = sample_indel_counts[each_sample]["in"]
            deletion_counts = sample_indel_counts[each_sample]["del"]
            indel_ratio = insertion_counts / deletion_counts
            out.write("{}\t{}\n".format(
                each_sample,
                indel_ratio
            ))


def main():
    parser = argparse.ArgumentParser(description="This script is designed to calculate insertion/deletion ratio")
    parser.add_argument("--path_indel_vcf", help="Path to indel vcf file", action="store", type=str, required=True)
    parser.add_argument("--path_out", help="Path to the output in/del result file", action="store", type=str,
                        required=True)
    parser.add_argument("--path_gt", help="Path to the GT format file by --extract-FORMAT-info GT", action="store",
                        type=str, required=True)
    args = parser.parse_args()
    if not os.path.exists(os.path.dirname(args.path_out)):
        os.makedirs(os.path.dirname(args.path_out))
    sample_indel_counts = read_in_gt(path_gt=args.path_gt, path_indel_vcf=args.path_indel_vcf)
    write_indel_ratio(sample_indel_counts=sample_indel_counts, path_out=args.path_out)


if __name__ == '__main__':
    main()

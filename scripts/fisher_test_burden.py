import os
import csv
import time
import random
import argparse
import numpy as np
from tqdm import tqdm
import multiprocessing
import scipy.stats as stats
from itertools import groupby
from operator import itemgetter
from collections import defaultdict
import statsmodels.stats.multitest as multi

"""
@Author: Kai Li
@email: lemonsky79@gmail.com
@Description: This script is designed to perform fisher exact test and permutation for each gene
"""


class InvalidSampleTypeError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class GenesParsingError(BaseException):
    def __init__(self, genes):
        self.genes = genes

    def __str__(self):
        return "Unexpected Error occurred when parsing the following genes {}".format(",".join(self.genes))


class GenesPermNumberNotEqualError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class EmptyResultError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


class SamplesNotFoundError(BaseException):
    def __init__(self, info):
        self.info = info

    def __str__(self):
        return self.info


def require_arg(args_list):
    # warn the user to input the required args
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values not in args_list:
                msg = "Args of {} must be chosen from {}".format(f"-{self.dest}", ",".join(args_list))
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)

    return RequireArgs


def require_argl():
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if type(values) != int:
                msg = "Args of {} must be integer".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            elif values < 1:
                msg = "Args of {} must be greater than or equal to 1".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            else:
                cpu_count = multiprocessing.cpu_count()
                if values > cpu_count:
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Local cpu count: {cpu_count}, "
                          f"processes user provide: {values} Adjusting to {cpu_count}")
                    # values = cpu_count
            setattr(args, self.dest, values)

    return RequireArgs


def get_case_control_samples(path_samples):
    # format of case case control sample is as follow
    # sample    sd_left sd_right    type
    # 19BY17889 -13 -12.25  case
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_samples}")
    case_samples = []
    control_samples = []
    for each_item in csv.DictReader(open(path_samples), delimiter="\t"):
        if each_item["type"] == "case":
            case_samples.append(each_item["sample"])
        elif each_item["type"] == "control":
            control_samples.append(each_item["sample"])
        else:
            raise InvalidSampleTypeError(f"Sample types must be case/control, got {each_item['type']}")
    return case_samples, control_samples


def fisher_perm_each_gene(each_gene, this_gene_case_sample_num,
                          this_gene_control_sample_num, case_samples, control_samples, hypothesis, pathway=None):
    # in this function, each_gene stands for a gene
    try:
        variant_sample_num = this_gene_case_sample_num + this_gene_control_sample_num

        def permutation():
            # perform permutation for each gene for 1000 times
            random_samples = random.sample(case_samples + control_samples, variant_sample_num)
            case_variant_intersect_samples = list(frozenset(random_samples).intersection(case_samples))
            control_variant_intersect_samples = list(frozenset(random_samples).intersection(control_samples))
            # yield the p value
            yield stats.fisher_exact([[len(case_variant_intersect_samples),
                                       len(case_samples) - len(case_variant_intersect_samples)],
                                      [len(control_variant_intersect_samples),
                                       len(control_samples) - len(control_variant_intersect_samples)]],
                                     alternative=hypothesis)[1]

        odds_ratio, this_gene_p_value = stats.fisher_exact([[this_gene_case_sample_num,
                                                             len(case_samples) - this_gene_case_sample_num],
                                                            [this_gene_control_sample_num,
                                                             len(control_samples) - this_gene_control_sample_num]],
                                                           alternative=hypothesis)
        # if input is a pathway list, we return the median of the 1000 permutation values
        # if no pathway is specified, we return the 1000 permutation values
        if pathway:
            this_gene_p_permutation = np.median([permutation().__next__() for _ in range(1000)])
        else:
            this_gene_p_permutation = [permutation().__next__() for _ in range(1000)]
        return (each_gene, this_gene_p_value, this_gene_p_permutation, odds_ratio)
    except Exception as e:
        return (each_gene, str(e), "error")


def read_in_snp_matrix(path_matrix):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_matrix}")
    samples = []
    total_genes = []
    samples_variants_count = defaultdict(list)
    genes_variant_samples = dict()
    with open(path_matrix) as file:
        for line in file:
            items = line.strip().split("\t")
            if line.startswith("variant"):
                # in the new vcf2mat.py output, samples start from the 8th column
                samples = items[7:]
            elif not samples:
                raise SamplesNotFoundError(f"Samples not found in {path_matrix}")
            else:
                gene = items[1]
                total_genes.append(gene)
                # in the new vcf2mat.py output, samples start from the 8th column
                genotypes = items[7:]
                # for each_sample_genotype in zip(samples, items[5:]):
                #    if each_sample_genotype[1] == "1" or each_sample_genotype[1] == "2":
                #        samples_variants_count[gene].append(each_sample_genotype[0])
                for k, g in groupby(zip(samples, genotypes), itemgetter(1)):
                    if k == "1" or k == "2":
                        samples_variants_count[gene].extend(list(list(zip(*g))[0]))
    total_genes = list(set(total_genes))
    for each_gene in total_genes:
        if each_gene in samples_variants_count:
            genes_variant_samples[each_gene] = list(set(samples_variants_count[each_gene]))
        else:
            genes_variant_samples[each_gene] = []
    return genes_variant_samples


def count_pathway_variant_samples(genes_variant_samples, path_pathway):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_pathway}")
    pathway_variant_samples = defaultdict(list)
    pathway_variant_samples_dedup = dict()
    total_pathways = []
    with open(path_pathway, encoding='utf-8') as file:
        for line in file:
            items = line.strip().split("\t")
            pathway = items[0]
            #if pathway != items[1]:
                #print(f"[warning]: Pathway 1st column not equal to the 2nd column: {items[0].decode('utf-8')} {items[1].decode('utf-8')}")
                #raise ValueError(f"1st column not same to the 2nd column in {path_pathway}")
            this_pathway_genes = items[2:]
            for each_gene in this_pathway_genes:
                if each_gene in genes_variant_samples:
                    pathway_variant_samples[pathway].extend(genes_variant_samples[each_gene])
            total_pathways.append(pathway)
    total_pathways = list(set(total_pathways))
    for each_pathway in total_pathways:
        if each_pathway in pathway_variant_samples:
            pathway_variant_samples_dedup[each_pathway] = list(set(pathway_variant_samples[each_pathway]))
        else:
            pathway_variant_samples_dedup[each_pathway] = []
    return pathway_variant_samples_dedup


def parse_snp_matrix(path_matrix, case_samples, control_samples, method, path_out, num_process, progress_bar,
                     hypothesis, path_pathway=None):
    genes_variant_samples = read_in_snp_matrix(path_matrix)
    if path_pathway:
        # to simplify the code, we use genes_variant_samples to represent the returned pathway variant samples
        # note that in the following genes_variant_samples, the keys are pathway terms rather than genes, the
        # value of each key stands for samples with snp
        genes_variant_samples = count_pathway_variant_samples(genes_variant_samples, path_pathway)
    genes = list(genes_variant_samples.keys())
    gene_case_control_samples = defaultdict(dict)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start performing fisher test and permutation")

    mp_list = []
    pbar = ""
    if progress_bar:
        pbar = tqdm(desc=f"  Fisher test/Permutation", total=len(genes))

    def call_back(result):
        if progress_bar:
            pbar.update()
        mp_list.append(result)

    pool = multiprocessing.Pool(processes=int(num_process))
    for each_gene in genes:
        this_gene_variant_samples = genes_variant_samples[each_gene]
        this_gene_variant_case_samples = list(frozenset(this_gene_variant_samples).intersection(case_samples))
        this_gene_variant_control_samples = list(frozenset(this_gene_variant_samples).intersection(control_samples))
        gene_case_control_samples[each_gene]["case"] = len(this_gene_variant_case_samples)
        gene_case_control_samples[each_gene]["control"] = len(this_gene_variant_control_samples)
        pool.apply_async(func=fisher_perm_each_gene, args=(each_gene, len(this_gene_variant_case_samples),
                                                           len(this_gene_variant_control_samples),
                                                           case_samples, control_samples, hypothesis, path_pathway,),
                         callback=call_back)

    pool.close()
    pool.join()

    if progress_bar:
        pbar.close()

    genes_pvals_pairs = dict()
    genes_with_error = dict()
    total_genes_permutation = []
    if not mp_list:
        raise EmptyResultError("No results got!")
    for each_item in mp_list:
        each_gene = each_item[0]
        if each_item[-1] == "error":
            genes_with_error[each_gene] = each_item[1]
        else:
            gene_case_control_samples[each_gene]["case_wild"] = len(case_samples) - gene_case_control_samples[each_gene]["case"]
            gene_case_control_samples[each_gene]["control_wild"] = len(control_samples) - gene_case_control_samples[each_gene]["control"]
            # each_item[3] is odds ratio
            gene_case_control_samples[each_gene]["odds_ratio"] = each_item[3]
            # each_item[1] is p value
            genes_pvals_pairs[each_gene] = each_item[1]
            if path_pathway:
                # if pathway is input, the each_item[2] is a number, if not, it's a python list
                gene_case_control_samples[each_gene]["pperm"] = each_item[2]
            else:
                total_genes_permutation.extend(sorted(each_item[2]))

    if genes_with_error:
        path_error = f"{path_out}.err"
        print(f"[ERROR]: Error occurred when performing fisher/permutation for certain genes, Writing failed parsing "
              f"genes to {path_error}")
        with open(path_error, "w+") as err:
            for each_gene in genes_with_error:
                err.write(f"{each_gene}\t{genes_with_error[each_gene]}\n")
        raise GenesParsingError(genes_with_error)

    # Because genes in the Manager.dict is out of order after the multiprocessing method, we have to reorder the genes
    # Calculate adjusted p value for each gene, and sort the genes according to it's p value by ascending order
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start sorting genes and adjusting p values")
    genes_sorted = sorted(genes_pvals_pairs, key=lambda k: genes_pvals_pairs[k])
    pvals_sorted = [genes_pvals_pairs[each_gene] for each_gene in genes_sorted]
    pvals_adjusted = dict(zip(genes_sorted, multi.multipletests(pvals_sorted, method=method)[1].tolist()))

    # if pathway is not input, merge all the p perm into a sorted list, and calculate the permutation for each gene by
    # each 1000 p perms in the list, if it's input, just assign total_genes_permutation to total_genes_permutation_mean
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start calculating p permutation")
    if path_pathway:
        total_genes_permutation_mean = [gene_case_control_samples[each_gene]["pperm"] for each_gene in genes_sorted]
    else:
        total_genes_permutation_mean = list(map(np.mean, np.array_split(sorted(total_genes_permutation), len(genes_sorted))))
    if len(total_genes_permutation_mean) != len(genes_sorted):
        raise GenesPermNumberNotEqualError("Gene numbers not equal to the calculated p permutation numbers")
    genes_pperm_pairs = dict(zip(genes_sorted, total_genes_permutation_mean))

    if not os.path.exists(os.path.dirname(path_out)):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Creating  directory {os.path.dirname(path_out)}")
        os.makedirs(os.path.dirname(path_out))

    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {path_out}")
    with open(path_out, "w+") as out:
        out.write("gene\tcase\tcase_wild\tcontrol\tcontrol_wild\tp_value\todds_ratio\tp_adj\tp_permutation\n")
        for each_gene in genes_sorted:
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(each_gene,
                                                                    gene_case_control_samples[each_gene]["case"],
                                                                    gene_case_control_samples[each_gene]["case_wild"],
                                                                    gene_case_control_samples[each_gene]["control"],
                                                                    gene_case_control_samples[each_gene]["control_wild"],
                                                                    genes_pvals_pairs[each_gene],
                                                                    gene_case_control_samples[each_gene]["odds_ratio"],
                                                                    pvals_adjusted[each_gene],
                                                                    genes_pperm_pairs[each_gene]))


def main():
    parser = argparse.ArgumentParser(description="This script is designed to perform fisher exact test and permutation")
    parser.add_argument("-p", help="Path to input matrix with columns representing samples and rows representing snps",
                        action="store", dest="p", type=str, required=True)
    parser.add_argument("-c", help="Path to input file with case and control samples, must be two columns, first col "
                                   "for samples, the other for case/control", action="store", dest="c", type=str,
                        required=True)
    parser.add_argument("-o", help="Path to put output matrix with the following information: "
                                   "gene\tcase\tcase_wild\tcontrol\tcontrol_wild\tp_value\tp_adj\tp_permutation",
                        action="store", dest="o", type=str, required=True)
    parser.add_argument("-m", help="Method to use in the adjustment of p values, must be chosen from "
                                   "bonferroni/sidak/holm-sidak/holm/simes-hochberg/hommel/fdr_bh/fdr_by/fdr_tsbh"
                                   "/fdr_tsbky, default is fdr_bh",
                        action=require_arg(["bonferroni", "sidak", "holm-sidak", "holm",
                                            "simes-hochberg", "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"]),
                        type=str, default="fdr_bh")
    parser.add_argument("-l", help="Processes to use in the fisher test and permutation step, default is 10",
                        action=require_argl(), type=int, default=10)
    parser.add_argument("-g", help="Whether to show the progress bar, default is off",
                        action="store_true", default=False)
    parser.add_argument("-s", help="Alternative hypothesis for fisher exact test, must be chosen from "
                                   "two-sided/less/greater, less and greater is one-sided",
                        action=require_arg(["two-sided", "less", "greater"]), type=str, default="two-sided")
    parser.add_argument("-w", help="Path to pathway table, each row stands for a pathway, the 1st and 2nd col are "
                                   "pathway name while the other cols are genes in this pathway. e.g. "
                                   "pathway\tpathway\tgene1\tgene2...",
                        action="store", default=None)
    args = parser.parse_args()
    case_samples, control_samples = get_case_control_samples(args.c)
    parse_snp_matrix(path_matrix=args.p, case_samples=case_samples, control_samples=control_samples, method=args.m,
                     path_out=args.o, num_process=args.l, progress_bar=args.g, hypothesis=args.s, path_pathway=args.w)


if __name__ == '__main__':
    main()

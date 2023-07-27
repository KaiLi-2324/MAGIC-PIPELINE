import re
import argparse
from collections import Counter


def format_genotype(genotype, reverse=False):
    genotype_items = genotype.split(":")
    if not reverse:
        if genotype_items[0] == "0/0":
            return "0"
        elif genotype_items[0] == "1/1":
            return "2"
        elif genotype_items[0] == "0/1" or genotype_items[0] == "1/0":
            return "1"
        else:
            return "."
    else:
        if genotype_items[0] == "0/0":
            # replace 0/0 to 1/1
            return "2"
        elif genotype_items[0] == "1/1":
            # replace 1/1 to 0/0
            return "0"
        elif genotype_items[0] == "0/1" or genotype_items[0] == "1/0":
            return "1"
        else:
            return "."


def parse_vcf(path_vcf, path_out, variant_type, e, cadd, LoF, frq_low, frq_up, path_corr):
    formats = []
    sample_names = []
    snp_out_format = ["0/0", "0/1", "1/1", "1/0", "./."]
    pattern_CSQ = re.compile(r";CSQ=(.*?)$")
    pattern_AC = re.compile(r"AC=(.*?);AF=(.*?);AN=(.*?);")
    missense_types = ["inframe_deletion", "inframe_insertion", "missense_variant", "stop_lost"]
    lof_types = ["stop_gained", "frameshift_variant", "start_lost", "splice_acceptor_variant", "splice_donor_variant"]
    print(f"Start reading from {path_vcf}")
    out = open(path_out, "w+")
    corr = open(path_corr, "w+")
    with open(path_vcf) as file:
        for line in file:
            if line.startswith("#"):
                if re.search(r'Format: (.*?)">', line):
                    print(f"Start parsing format line")
                    formats = re.findall(r'Format: (.*?)">', line)[0].split("|")
                elif line.startswith("#CHROM"):
                    print(f"Start collecting sample names")
                    sample_names = line.strip().split("\t")[9:]
                    out.write("variant\tsymbol\tENSG\tENST\tsample_stat\tvariant_stat\tinfo\t" + "\t".join(sample_names) + "\n")
            elif not formats:
                raise ValueError(f"Format line not found in {path_vcf}")
            elif not sample_names:
                raise ValueError(f"Sample name not found in {path_vcf}")
            else:
                items = line.strip().split("\t")
                # print(pattern_AC.findall(items[7]))
                ac, af, an = list(map(float, pattern_AC.findall(items[7])[0]))
                samples_info_items = list(map(lambda k: k.split(":")[0], items[9:]))
                # count the frequency of "0/0", "1/1", "1/0", "0/1" and "./."
                snp_type_counts = Counter(samples_info_items)
                snp_type_out = []
                for each_snp in snp_out_format:
                    if each_snp in snp_type_counts:
                        snp_type_out.append(str(snp_type_counts[each_snp]))
                    else:
                        snp_type_out.append("0")
                # whether to output the corrected variants, if true, the variant will be output
                whether_output_corrected = False
                if af <= 0.05:
                    genotypes = list(map(format_genotype, items[9:]))
                else:
                    whether_output_corrected = True
                    genotypes = list(map(lambda k: format_genotype(k, reverse=True), items[9:]))
                transcripts = pattern_CSQ.findall(items[7])[0].split(",")
                if len(transcripts) > 1:
                    raise ValueError(f"More than one transcript detected in {items[7]}")
                this_snp_trans_type = ""
                each_transcript = transcripts[0]
                this_transcript_items = each_transcript.split("|")
                # replace the empty string in this_transcript_items with 0
                this_transcript_format_dict_rep0 = dict(zip(formats, list(map(lambda k: k if k else 0, this_transcript_items))))
                genom_af = this_transcript_format_dict_rep0["gnomAD_AF"]
                if any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, lof_types))):
                    if LoF and this_transcript_format_dict_rep0["LoF"] == LoF:
                        this_snp_trans_type = "lof"
                    else:
                        this_snp_trans_type = "lof"
                elif any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, missense_types))):
                    if (re.search(r"^deleterious\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^probably_damaging\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"]))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) > cadd:
                        this_snp_trans_type = "missense_damage"
                    if (re.search(r"^tolerated\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^benign\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"]))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) < cadd:
                        this_snp_trans_type = "missense_benign"
                if this_transcript_format_dict_rep0["Consequence"] == "synonymous_variant" and float(this_transcript_format_dict_rep0["CADD_PHRED"]) < cadd:
                    this_snp_trans_type = "synonymous"
                # if not float(this_transcript_format_dict_rep0[e]) < frequency:
                #    continue
                # variant_stat AC/AF/AN/gnomad_AF
                this_snp_variant_stat = f"{ac}/{af}/{an}/{genom_af}"
                # whether to output this variant
                whether_output = False
                if not e:
                    whether_output = True
                elif float(frq_low) < float(this_transcript_format_dict_rep0[e]) <= float(frq_up):
                    whether_output = True
                if whether_output:
                    if variant_type == this_snp_trans_type:
                        if whether_output_corrected:
                            corr.write(line)
                        out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                            f"{items[0]}_{items[1]}_{items[3]}_{items[4]}",
                            this_transcript_format_dict_rep0['SYMBOL'],
                            this_transcript_format_dict_rep0['Gene'],
                            this_transcript_format_dict_rep0['Feature'],
                            "/".join(snp_type_out),
                            this_snp_variant_stat,
                            each_transcript,
                            "\t".join(genotypes)
                        ))
    out.close()
    corr.close()


def main():
    parser = argparse.ArgumentParser(description="This script is to convert vcf to a matrix")
    parser.add_argument("-p", help="Path to reference file", action="store", dest="path_vcf", type=str, required=True)
    parser.add_argument("-o", help="Path to output file", action="store", dest="path_out", type=str, required=True)
    parser.add_argument("-e", help="Variable to filter, must be chosen from MAX_AF, MAX_AF_POPS, AFR_AF, AMR_AF, "
                                   "EAS_AF, EUR_AF, SAS_AF, AA_AF, EA_AF, gnomAD_AF, gnomAD_AFR_AF, gnomAD_AMR_AF, "
                                   "gnomAD_ASJ_AF, gnomAD_EAS_AF, gnomAD_FIN_AF, gnomAD_NFE_AF, gnomAD_OTH_AF, "
                                   "gnomAD_SAS_AF, default is not given. If not given, won't use the database to "
                                   "filter variants", action="store", dest="e", type=str, default="",
                        choices=["MAX_AF", "MAX_AF_POPS", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF", "AA_AF",
                                 "EA_AF", "gnomAD_AF", "gnomAD_AFR_AF", "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "",
                                 "gnomAD_EAS_AF", "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF", "gnomAD_SAS_AF"])
    parser.add_argument("-f", help="Lower limit for -e, a float number to filter the variants by database, default is "
                                   "negative infinity. eg. when -f 0, it means that if gnomAD_AF > 0, then the "
                                   "transcript is output", action="store", dest="frq_low", type=str, default="-inf")
    parser.add_argument("-F", help="Upper limit for -e, a float number to filter the variants by database, default is "
                                   "positive infinity. eg. when -F 0.05, it means that if gnomAD_AF <= 0.05, then the "
                                   "transcript is output", action="store", dest="frq_up", type=str, default="inf")
    parser.add_argument("-v", help="Variant types to filter the vcf, must be chosen from lof/missense_damage/"
                                   "missense_benign/synonymous/3UTR/5UTR, default is lof", action="store",
                        dest="variant", type=str, default="lof",
                        choices=["lof", "missense_damage", "missense_benign", "synonymous", "3UTR", "5UTR"])
    parser.add_argument("-l", help="Chosen from HC/LC, default is not given, if given, HC or LC WON'T be used to filter"
                                   " for lof variant", action="store", dest="lof", type=str, choices=["HC", "LC", ""],
                        default="")
    parser.add_argument("-s", help="Threshold for CADD_PHRED, used in missense_damage/missense_benign/synonymous, "
                                   "default is 15", action="store", dest="cadd", type=str, default=15)
    parser.add_argument("-c", help="Path to put corrected variants", action="store", dest="path_corr", type=str,
                        required=True)
    args = parser.parse_args()
    if args.e:
        if float(args.frq_low) > float(args.frq_up):
            raise ValueError(f"Lower limit for {args.e} can not be greater than upper limit!")
    parse_vcf(path_vcf=args.path_vcf, path_out=args.path_out, variant_type=args.variant, e=args.e, cadd=args.cadd,
              LoF=args.lof, frq_low=args.frq_low, frq_up=args.frq_up, path_corr=args.path_corr)


if __name__ == '__main__':
    main()

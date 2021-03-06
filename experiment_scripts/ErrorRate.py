#!/usr/bin/env python
"""ErrorRate.py randomly selects alignments from two experiments and tests their accuracy at a given coverage
The input to this script is a directory of PCR generated alignments and a directory of genomic DNA generated
alignments. The program then randomly selects a given number of alignments and calls VClr to make variant
calls on the sites. The script then records the number of sites called as methyl/non-methyl and compares it
to ground truth (methylated in the genomic DNA case, non-methyl in the PCR case). It then records the error
rate and repeats with more alignments (without replacement).

"""
from __future__ import print_function, division
import sys
import glob
import subprocess as sp
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from random import shuffle


def parse_reads_probs(alignment_files):
    lines = []
    assert len(alignment_files) > 0, "Didn't provide any alignment files"
    for f in alignment_files:
        for line in open(f, "r"):
            lines.append(line)
    assert len(lines) > 0, "Problem parsing files"
    return lines


def variant_call_stats(probs, read_score, strand, threshold):
    assert len(probs) > 0, "Probs length 0"
    vclr_command = ["VClr", "-tool=methyl"]
    if strand is not None:
        vclr_command.append("-strand={}".format(strand))
    if threshold is not None:
        vclr_command.append("-t={}".format(threshold))
    if read_score is not None:
        vclr_command.append("-s={}".format(read_score))

    vclr = sp.Popen(vclr_command, stdin=sp.PIPE, stdout=sp.PIPE)

    for l in probs:
        vclr.stdin.write(l)
    output = vclr.communicate()[0]
    output = output.split("\n")[1:-1]

    vc_dat = {
        "site": [],
        "call": [],
        "coverage": [],
    }

    for l in output:
        l = l.strip().split()
        vc_dat["site"].append(int(l[0]))
        vc_dat["call"].append(str(l[1]))
        vc_dat["coverage"].append(int(l[2]))

    vc_dat = pd.DataFrame(vc_dat)

    methyl_calls = len(vc_dat[(vc_dat["call"] == "E") | (vc_dat["call"] == "I")])
    canonical_calls = len(vc_dat[(vc_dat["call"] == "C") | (vc_dat["call"] == "A")])
    mean_coverage = np.mean(vc_dat["coverage"])

    return methyl_calls, canonical_calls, mean_coverage


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-pcr", action="store", dest="pcr", required=True)
        parser.add_argument("-gen", action="store", dest="genomic", required=True)
        parser.add_argument("-N", action="store", dest="N", type=int, required=False, default=10)
        parser.add_argument("-i", action="store", dest="iter", type=int, required=False, default=None)
        parser.add_argument("-t", action="store", dest="threshold", required=False, default=None)
        parser.add_argument("-strand", action="store", dest="strand", required=False, default=None)
        parser.add_argument("-s", action="store", dest="read_score", required=False, default=None)
        args = parser.parse_args()
        return args
    args = parse_args()

    print("accuracy\tsensitivity\tspecificity\tprecision\tfallout\tmiss_rate\tFDR\tNPV\tmean_coverage")

    if args.iter is not None:
        i = 0
        block_size = args.N

        # collect alignments from the PCR directory
        pcr_alignments = [x for x in glob.glob(args.pcr)]
        shuffle(pcr_alignments)

        # collect alignments from genomic directory
        gen_alignments = [x for x in glob.glob(args.genomic)]
        shuffle(gen_alignments)

        if len(gen_alignments) >= len(pcr_alignments):
            limit = len(pcr_alignments)
        else:
            limit = len(gen_alignments)

        while i < args.iter:
            block_start = i * block_size
            block_end = block_start + block_size
            if block_end >= limit:
                break
            pcr_batch = pcr_alignments[block_start:block_end]
            gen_batch = gen_alignments[block_start:block_end]
            pcr_probs = parse_reads_probs(pcr_batch)
            gen_probs = parse_reads_probs(gen_batch)

            false_positive, true_negative, pcr_coverage = variant_call_stats(probs=pcr_probs,
                                                                             strand=args.strand,
                                                                             threshold=args.threshold,
                                                                             read_score=args.read_score)
            true_positive, false_negative, gen_coverage = variant_call_stats(probs=gen_probs,
                                                                             strand=args.strand,
                                                                             threshold=args.threshold,
                                                                             read_score=args.read_score)

            accuracy = (true_positive + true_negative) / (false_positive + true_negative +
                                                          true_positive + false_negative)

            sensitivity = true_positive / (true_positive + false_negative)  # hit rate, recall

            specificity = true_negative / (true_negative + false_positive)  # true negative rate

            precision = true_positive / (true_positive + false_positive)    # positive predictive value

            fall_out = false_positive / (false_positive + true_negative)    # false positive rate

            miss_rate = false_negative / (true_positive + false_negative)   # false negative rate

            fdr = false_positive / (true_positive + false_positive)         # false discovery rate

            npv = true_negative / (true_negative + false_negative)          # negative predictive value

            mean_cov = np.mean([pcr_coverage, gen_coverage])

            print(accuracy, sensitivity, specificity, precision, fall_out, miss_rate, fdr, npv, mean_cov, sep="\t", end="\n")
            i += 1
    else:
        pcr_dat = parse_reads_probs([args.pcr])
        genomic_dat = parse_reads_probs([args.genomic])
        false_positive, true_negative, pcr_coverage = variant_call_stats(probs=pcr_dat,
                                                                         strand=args.strand,
                                                                         threshold=args.threshold,
                                                                         read_score=args.read_score)

        true_positive, false_negative, gen_coverage = variant_call_stats(probs=genomic_dat,
                                                                         strand=args.strand,
                                                                         threshold=args.threshold,
                                                                         read_score=args.read_score)

        accuracy = (true_positive + true_negative) / (false_positive + true_negative + true_positive + false_negative)

        sensitivity = true_positive / (true_positive + false_negative)  # hit rate, recall

        specificity = true_negative / (true_negative + false_positive)  # true negative rate

        precision = true_positive / (true_positive + false_positive)  # positive predictive value

        fall_out = false_positive / (false_positive + true_negative)  # false positive rate

        miss_rate = false_negative / (true_positive + false_negative)  # false negative rate

        fdr = false_positive / (true_positive + false_positive)  # false discovery rate

        npv = true_negative / (true_negative + false_negative)          # negative predictive value

        mean_cov = np.mean([pcr_coverage, gen_coverage])

        print(accuracy, sensitivity, specificity, precision, fall_out, miss_rate, fdr, npv, mean_cov, sep="\t", end="\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv))

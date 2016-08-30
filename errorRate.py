#!/usr/bin/env python
"""errorRate.py randomly selects alignments from two experiments and tests their accuracy at a given coverage
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


def variant_call_stats(probs, read_score=40):
    assert len(probs) > 0, "Probs length 0"
    vclr = sp.Popen(["VClr", "-tool=methyl", "-s={}".format(read_score)], stdin=sp.PIPE, stdout=sp.PIPE)
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
        parser.add_argument("-N", action="store", dest="N", type=int, required=True)
        parser.add_argument("-i", action="store", dest="iter", type=int, required=True)
        args = parser.parse_args()
        return args
    args = parse_args()

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

    i = 0
    block_size = args.N
    while i <= args.iter:
        block_start = i * block_size
        block_end = block_start + block_size
        if block_end >= limit:
            break
        pcr_batch = pcr_alignments[block_start:block_end]
        gen_batch = gen_alignments[block_start:block_end]

        pcr_false_positive, pcr_true_negative, pcr_coverage = variant_call_stats(parse_reads_probs(pcr_batch))
        gen_true_positive, gen_false_negative, gen_coverage = variant_call_stats(parse_reads_probs(gen_batch))

        accuracy = (gen_true_positive + pcr_true_negative) / (pcr_false_positive + pcr_true_negative +
                                                              gen_true_positive + gen_false_negative)

        print(accuracy)
        i += 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))

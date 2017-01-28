#!/usr/bin/env python
"""This script takes a directory of assignments, which are tab-seperated files containing k-mer to ionic
current pairs with a posterior probability it consolidates them, and reformats the whole thing for 
HdpBuildUtil and runs it"""
from __future__ import print_function
import sys
import os
from argparse import ArgumentParser
from BuildModels import make_master_assignment_table, make_bulk_build_alignment


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-pcr", action="store", dest="pcr", required=True)
        parser.add_argument("-gen", action="store", dest="gen", required=True)
        parser.add_argument("-o", action="store", dest="out_dir", required=True)

        parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=30)
        parser.add_argument("-c", action="store", dest="methyl_assignments", required=False,
                            type=int, default=200)
        parser.add_argument("-t", action="store", dest="assignment_threshold", required=False, type=float, default=0.8)
        return parser.parse_args()

    args = parse_args()
    print("Commandline: %s " % " ".join(sys.argv[:]), file=sys.stderr)
    
    assignments = make_master_assignment_table([args.pcr, args.gen])
    
    out_aln_filename    = "build_alignment_CG_{mC}_{C}_{t}".format(mC=args.methyl_assignments,
                                                                   C=args.assignments,
                                                                   t=args.assignment_threshold)
    out_build_alignment = os.path.join(args.out_dir, out_aln_filename)

    make_bulk_build_alignment(assignments, "cytosine2",
                              n_canonical_assignments=args.assignments,
                              n_methyl_assignments=args.methyl_assignments,
                              threshold=args.assignment_threshold,
                              outfile=out_build_alignment,
                              strands=["t"])
    sys.exit(0)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

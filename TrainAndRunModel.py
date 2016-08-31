#!/usr/bin/env python
"""TrainAndRunModel.py trains the transitions for an HMM-HDP then runs a classification experiment
"""
from __future__ import print_function
import sys
import os
from argparse import ArgumentParser
from BuildModels import train_model_transitions, run_variant_calling_experiment

PATH_TO_SIGNALALIGN = os.path.abspath("../signalAlign/")
PATH_TO_BINS = PATH_TO_SIGNALALIGN + "/bin/"


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-r", action="store", dest="reference", required=True)
        parser.add_argument("--pcr_train", action="store", dest="pcr_fofn_train", required=True)
        parser.add_argument("--gen_train", action="store", dest="genomic_fofn_train", required=True)
        parser.add_argument("--pcr_test", action="store", dest="pcr_fofn_test", required=True)
        parser.add_argument("--gen_test", action="store", dest="genomic_fofn_test", required=True)
        parser.add_argument('--positions', action='append', dest='positions_file', required=True)
        parser.add_argument('--motif', action='append', dest='motif_file', required=True)
        parser.add_argument("-tH", action="store", dest="t_hdp", required=True)
        parser.add_argument("-cH", action="store", dest="c_hdp", required=True)
        parser.add_argument("-T", action="store", dest="t_hmm", required=True)
        parser.add_argument("-C", action="store", dest="c_hmm", required=True)
        parser.add_argument("-a", action="store", dest="batch", type=int, required=True)
        parser.add_argument("-i", action="store", dest="iterations", type=int, required=True)
        parser.add_argument("-n", action="store", dest="n_test_alns", type=int, required=True)
        parser.add_argument("-x", action="store", dest="degenerate", required=True)
        parser.add_argument("-j", action="store", dest="jobs", required=True)
        parser.add_argument("-o", action="store", dest="outpath", required=True)
        args = parser.parse_args()
        return args

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    args = parse_args()
    working_path = os.path.abspath(args.outpath)

    assert len(args.positions_file) == 2 and len(args.motif_file) == 2, "need to give training and testing " \
                                                                        "positions/motif files"
    for i in range(2):
        assert os.path.exists(args.positions_file[i]), "Didn't find positions file, looked " \
                                                       "{}".format(args.positions_file)
        assert os.path.exists(args.motif_file[i]), "Didn't find motif file, looked {}".format(args.motif_file)
    positions_file = args.positions_file[0]
    #motif_file = args.motif_file[0]
    test_positions = args.positions_file[1]
    test_motifs = args.motif_file[1]

    models = train_model_transitions(fasta=os.path.abspath(args.reference),
                                     pcr_fofn=args.pcr_fofn_train,
                                     genomic_fofn=args.genomic_fofn_train,
                                     degenerate=args.degenerate,
                                     jobs=args.jobs,
                                     positions_file=positions_file,
                                     iterations=args.iterations,
                                     batch_size=args.batch,
                                     outpath=working_path,
                                     stateMachine="threeStateHdp",
                                     t_hdp=args.t_hdp,
                                     c_hdp=args.c_hdp,
                                     t_model=os.path.abspath(args.t_hmm),
                                     c_model=os.path.abspath(args.c_hmm))

    run_variant_calling_experiment(fasta=os.path.abspath(args.reference),
                                   pcr_fofn=args.pcr_fofn_test,
                                   genomic_fofn=args.genomic_fofn_test,
                                   jobs=args.jobs,
                                   positions_file=test_positions,
                                   motif_file=test_motifs,
                                   t_model=models[0],
                                   c_model=models[1],
                                   t_hdp=models[2],
                                   c_hdp=models[3],
                                   outpath=working_path,
                                   n=args.n_test_alns,
                                   degenerate=args.degenerate)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

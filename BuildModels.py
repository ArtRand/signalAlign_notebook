#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import glob
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from subprocess import Popen
from itertools import product
from commonFunctions import get_first_seq, make_motif_file, get_all_sequence_kmers

PATH_TO_SIGNALALIGN = os.path.abspath("../signalAlign/")
PATH_TO_BINS = PATH_TO_SIGNALALIGN + "/bin/"


def gatc_motifs(core="GITC"):
    motifs = []
    for _ in product("ACGT", repeat=1):
        motifs.append(_[0] + core)
        motifs.append(core + _[0])
    return motifs


def find_gatc_motifs(sequence):
    motif = "GATC"
    motif_length = len(motif)
    for i, _ in enumerate(sequence):
        if sequence[i:i+motif_length] == motif:
            yield i + 1


def make_gatc_position_file(fasta, outfile):
    outfile = os.path.abspath(outfile)
    fH = open(outfile, 'w')
    fH.write("X\t")

    seq = get_first_seq(fasta)
    for i in find_gatc_motifs(seq):
        assert seq[i] == "A"
        fH.write("{}\t".format(i))
    fH.write("\n")
    fH.write("X\t")
    for i in find_gatc_motifs(seq):
        t_pos = i + 1
        assert seq[t_pos] == "T"
        fH.write("{}\t".format(t_pos))  # + 1 because this is for the reverse complement
    fH.write("\n")
    fH.close()
    return outfile


def make_gatc_motif_file(fasta, outfile):
    outfile = os.path.abspath(outfile)
    seq = get_first_seq(fasta)
    gatcs = [x for x in find_gatc_motifs(seq)]
    make_motif_file(gatcs, seq, outfile)
    return outfile


def run_guide_alignment(fasta, pcr_reads, genomic_reads, jobs, positions_file, motif_file, t_model, c_model,
                        outpath, n):
    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -d={reads} -r={fasta} -T={tModel} -C={cModel} -f=assignments " \
                              "-o={outpath} -p={positions} -q={targetFile} -X={sub} -n={n} -j={jobs}"

    read_sets = [pcr_reads, genomic_reads]
    labels = ["A", "I"]
    working_directories = [working_path + "/pcr_", working_path + "/genomic_"]

    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)

    for reads, label, working_directory in zip(read_sets, labels, working_directories):
        # assemble the command
        command = c.format(reads=reads, fasta=fasta, tModel=t_model, cModel=c_model, outpath=working_directory,
                           positions=positions_file, targetFile=motif_file, sub=label, n=n, jobs=int(jobs/2))
        commands.append(command)
    os.chdir(PATH_TO_BINS)
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in commands]
    status = [p.wait() for p in procs]
    os.chdir(working_path)
    working_directories = [d + "tempFiles_alignment/*.assignments" for d in working_directories]

    return working_directories


def run_variant_calling_experiment(fasta, pcr_reads, genomic_reads, jobs, positions_file, motif_file, t_model, c_model,
                                  outpath, n, degenerate, t_hdp, c_hdp):
    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -d={reads} -r={fasta} -T={tModel} -C={cModel} -f=variantCaller " \
                       "-o={outpath} -p={positions} -q={targetFile} -n={n} -j={jobs} -x={degenerate} " \
                       "-tH={tHdp} -cH={cHdp} -smt=threeStateHdp"
    read_sets = [pcr_reads, genomic_reads]
    working_directories = [working_path + "/vc_pcr_", working_path + "/vc_genomic_"]

    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)
    assert os.path.exists(t_hdp), "Didn't find template model, looked {}".format(t_hdp)
    assert os.path.exists(c_hdp), "Didn't find complement model, looked {}".format(c_hdp)

    for reads, working_directory in zip(read_sets, working_directories):
        # assemble the command
        command = c.format(reads=reads, fasta=fasta, tModel=t_model, cModel=c_model, outpath=working_directory,
                           positions=positions_file, targetFile=motif_file, n=n, jobs=int(jobs/2), tHdp=t_hdp,
                           cHdp=c_hdp)
        commands.append(command)

    os.chdir(PATH_TO_BINS)
    procs = [Popen(x.split(), stdout=sys.stdout, stderr=sys.stderr) for x in commands]
    status = [p.wait() for p in procs]
    os.chdir(working_path)
    return

def make_master_assignment_table(assignment_directories):
    def parse_assignment_file(file):
        data = pd.read_table(file,
                             usecols=(0, 1, 2, 3),
                             names=["kmer", "strand", "level_mean", "prob"],
                             dtype={"kmer": np.str, "strand": np.str, "level_mean": np.float64, "prob": np.float64},
                             header=None
                             )
        return data

    assignments = []
    for d in assignment_directories:
        assignments += [x for x in glob.glob(d) if os.stat(x).st_size != 0]
    assignment_dfs = []
    for f in assignments:
        assignment_dfs.append(parse_assignment_file(f))
    return pd.concat(assignment_dfs)


def train_model_transitions(fasta, pcr_reads, genomic_reads, jobs, positions_file, iterations, batch_size, outpath,
                            stateMachine="threeState", t_hdp=None, c_hdp=None):
    working_path = os.path.abspath(outpath) + "/"
    model_directory = working_path + "{}_".format(stateMachine)
    t_model = PATH_TO_SIGNALALIGN + "/models/testModelR9_5mer_acegit_template.model"
    c_model = PATH_TO_SIGNALALIGN + "/models/testModelR9_5mer_acegit_complement.model"
    assert os.path.exists(t_model), "Didn't find template model"
    assert os.path.exists(c_model), "Didn't find complement model"

    os.chdir(PATH_TO_BINS)
    c = "trainModels -d={pcr} -d={genomic} -X=C -X=I -r={fasta} -i={iter} -a={batch} --transitions -smt={smt} " \
        "-T={tModel} -C={cModel} -j={jobs} -x=adenosine -p={positions} -o={out} " \
        "".format(pcr=pcr_reads, genomic=genomic_reads, fasta=fasta, iter=iterations, batch=batch_size,
                  smt=stateMachine, tModel=t_model, cModel=c_model, jobs=jobs,
                  positions=positions_file, out=model_directory)
    if t_hdp is not None and c_hdp is not None:
        c += "-tH={tHdp} -cH={cHdp} ".format(tHdp=os.path.abspath(t_hdp), cHdp=os.path.abspath(c_hdp))
    c = PATH_TO_BINS + c
    os.system(c)
    models = [model_directory + "tempFiles_expectations/template_trained.hmm",
              model_directory + "tempFiles_expectations/complement_trained.hmm"]
    if t_hdp is not None and c_hdp is not None:
        models.append(model_directory + "tempFiles_expectations/template.singleLevelPriorEcoli.nhdp")
        models.append(model_directory + "tempFiles_expectations/complement.singleLevelPriorEcoli.nhdp")
    os.chdir(working_path)
    return models


def make_build_alignment(assignments, ref_fasta, num_assignments, outfile, threshold=0.1):
    seq = get_first_seq(ref_fasta)
    kmers = get_all_sequence_kmers(seq, 5).keys()
    kmers += gatc_motifs("GITC")
    fH = open(outfile, "w")
    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"
    for strand in ["t", "c"]:
        by_stand = assignments.ix[(assignments['strand'] == strand) & (assignments['prob'] >= threshold)]
        for k in kmers:
            kmer_assignments = by_stand.ix[by_stand['kmer'] == k]
            kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)
            n = 0
            for _, r in kmer_assignments.iterrows():
                fH.write(entry_line.format(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob']))
                n += 1
                if n >= num_assignments:
                    break
    fH.close()
    return outfile


def build_hdp(build_alignment_path, template_model, complement_model, outpath, samples=15000):
    build_alignment = os.path.abspath(build_alignment_path)
    t_model = os.path.abspath(template_model)
    c_model = os.path.abspath(complement_model)
    outpath = os.path.abspath(outpath) + "/"
    hdp_pipeline_dir = outpath + "hdpPipeline/"
    os.makedirs(hdp_pipeline_dir)
    os.chdir(PATH_TO_BINS)
    c = "hdp_pipeline --build_alignment={build} -tM={tModel} -cM={cModel} -Ba=1 -Bb=1 -Ma=1 -Mb=1 -La=1 -Lb=1 " \
        "-s={samples} --verbose --grid_start=60 --grid_end=130 --grid_length=1400 --verbose -o={out} " \
        "--hdp_type=ecoli".format(build=build_alignment, tModel=t_model, cModel=c_model, samples=samples,
                                  out=hdp_pipeline_dir)
    c = PATH_TO_BINS + c
    os.system(c)
    os.chdir(outpath)
    return [hdp_pipeline_dir + "template.singleLevelPriorEcoli.nhdp",
            hdp_pipeline_dir + "complement.singleLevelPriorEcoli.nhdp"]


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-r", action="store", dest="reference", required=True)
        parser.add_argument("-pcr", action="store", dest="pcr_reads", required=True)
        parser.add_argument("-gen", action="store", dest="genomic_reads", required=True)
        parser.add_argument("-o", action="store", dest="outpath", required=True)
        parser.add_argument("-j", action="store", dest="jobs", required=False, default=4, type=int)
        parser.add_argument("-i", action="store", dest="iterations", required=False, type=int, default=20)
        parser.add_argument("-a", action="store", dest="batch", required=False, type=int, default=15000)
        parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=100)
        parser.add_argument("-n", action="store", dest="n_aligns", required=False, type=int, default=400)
        parser.add_argument("-e", action="store", dest="n_test_alns", required=False, type=int, default=400)
        parser.add_argument("-t", action="store", dest="assignment_threshold", required=False, type=float, default=0.3)
        parser.add_argument("-g", action="store", dest="samples", required=False, type=int, default=15000)
        args = parser.parse_args()
        return args
    args = parse_args()

    # make the positions and motif file
    working_path = os.path.abspath(args.outpath)
    positions_file = make_gatc_position_file(args.reference, working_path + "/gatc_positions.positions")

    # make the motif file
    motif_file = make_gatc_motif_file(args.reference, working_path + "/gatc_target.target")

    # train the transitions
    models = train_model_transitions(fasta=os.path.abspath(args.reference),
                                     pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                     genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                     jobs=args.jobs,
                                     positions_file=positions_file,
                                     iterations=args.iterations,
                                     batch_size=args.batch,
                                     outpath=working_path,
                                     )
    # do the initial alignments
    assignment_dirs = run_guide_alignment(fasta=os.path.abspath(args.reference),
                                          pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                          genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                          jobs=args.jobs,
                                          positions_file=positions_file,
                                          motif_file=motif_file,
                                          n=args.n_aligns,
                                          t_model=models[0],
                                          c_model=models[1],
                                          outpath=working_path)
    # concatenate the assignments into table
    master = make_master_assignment_table(assignment_dirs)
    build_alignment = make_build_alignment(assignments=master,
                                           ref_fasta=os.path.abspath(args.reference),
                                           num_assignments=args.assignments,
                                           outfile=working_path + "/buildAlignment.tsv",
                                           threshold=args.assignment_threshold)
    # build hdp
    hdps = build_hdp(build_alignment_path=build_alignment,
                     template_model=models[0],
                     complement_model=models[1],
                     outpath=working_path,
                     samples=args.samples)
    # train HMM/HDP
    hdp_models = train_model_transitions(fasta=os.path.abspath(args.reference),
                                         pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                         genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                         jobs=args.jobs,
                                         positions_file=positions_file,
                                         iterations=args.iterations,
                                         batch_size=args.batch,
                                         outpath=working_path,
                                         stateMachine="threeStateHdp",
                                         t_hdp=hdps[0],
                                         c_hdp=hdps[1],
                                         )
    # run methylation variant calling experiment
    run_variant_calling_experiment(fasta=os.path.abspath(args.reference),
                                   pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                   genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                   jobs=args.jobs,
                                   positions_file=positions_file,
                                   motif_file=motif_file,
                                   t_model=hdp_models[0],
                                   c_model=hdp_models[1],
                                   outpath=working_path,
                                   n=args.n_test_alns,
                                   degenerate="adenosine",
                                   t_hdp=hdp_models[2],
                                   c_hdp=hdp_models[3])

if __name__ == "__main__":
    sys.exit(main(sys.argv))







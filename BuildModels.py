#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import glob
import pandas as pd
import numpy as np
from subprocess import Popen, check_call
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
    return


def make_gatc_motif_file(fasta, outfile):
    seq = get_first_seq(fasta)
    gatcs = [x for x in find_gatc_motifs(seq)]
    make_motif_file(gatcs, seq, outfile)


def run_guide_alignment(fasta, pcr_reads, genomic_reads, jobs, positions_file, motif_file, n, t_model, c_model,
                        outpath):
    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -d={reads} -r={fasta} -T={tModel} -C={cModel} -f=assignments " \
                              "-o={outpath} -p={positions} -q={targetFile} -X={sub} -n={n} -j={jobs}"

    read_sets = [pcr_reads, genomic_reads]
    labels = ["A", "I"]
    working_directories = [working_path + "/pcr_", working_path + "/genomic_"]

    assert os.path.exists(t_model), "Didn't find template model"
    assert os.path.exists(c_model), "Didn't find complement model"

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

    t_model = PATH_TO_SIGNALALIGN + "/models/testModelR9_5mer_acegit_template.model"
    c_model = PATH_TO_SIGNALALIGN + "/models/testModelR9_5mer_acegit_complement.model"

    os.chdir(PATH_TO_BINS)
    c = "trainModels -d={pcr} -d={genomic} -X=C -X=I -r={fasta} -i={iter} -a={batch} --transitions -smt={smt} " \
        "-T={tModel} -C={cModel} -j={jobs} -x=adenosine -p={positions} -o={out} " \
        "".format(pcr=pcr_reads, genomic=genomic_reads, fasta=fasta, iter=iterations, batch=batch_size,
                  smt=stateMachine, tModel=t_model, cModel=c_model, jobs=jobs,
                  positions=os.path.abspath(positions_file), out=working_path)
    if t_hdp is not None and c_hdp is not None:
        c += "-tH={tHdp} -cH={cHdp} ".format(tHdp=os.path.abspath(t_hdp), cHdp=os.path.abspath(c_hdp))
    c = PATH_TO_BINS + c
    os.system(c)
    models = [working_path + "tempFiles_expectations/template_trained.hmm",
              working_path + "tempFiles_expectations/complement_trained.hmm"]
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
    return


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
    return


def main(args):
    # train the transitions
    models = train_model_transitions("/Users/Rand/projects/marginAlign/notebooks/pUC19_SspI.fa",
                                     "/Users/Rand/projects/marginAlign/sandbox/R9_pUC/NB07_pcr/",
                                     "/Users/Rand/projects/marginAlign/sandbox/R9_pUC/NB08_genomic/", 4,
                                     "/Users/Rand/projects/marginAlign/notebooks/pUC_gatc.positions",
                                     1, 2000, "./")
    # do the initial alignments
    assignment_dirs = run_guide_alignment(fasta="/Users/Rand/projects/marginAlign/notebooks/pUC19_SspI.fa",
                                          pcr_reads="/Users/Rand/projects/marginAlign/sandbox/R9_pUC/NB07_pcr/",
                                          genomic_reads="/Users/Rand/projects/marginAlign/sandbox/R9_pUC/NB08_genomic/",
                                          jobs=4,
                                          positions_file="/Users/Rand/projects/marginAlign/notebooks/pUC_gatc.positions",
                                          motif_file="/Users/Rand/projects/marginAlign/notebooks/pUC_gatc.target",
                                          n=2,
                                          t_model=models[0],
                                          c_model=models[1],
                                          outpath="/Users/Rand/projects/marginAlign/notebooks/")

    master = make_master_assignment_table(assignment_dirs)
    make_build_alignment(master, "/Users/Rand/projects/marginAlign/notebooks/pUC19_SspI.fa", 3, "./buildAlignment.tsv")
    build_hdp("./buildAlignment.tsv", models[0], models[1], "./", 200)

if __name__ == "__main__":
    sys.exit(main(sys.argv))







#!/usr/bin/env python
from __future__ import print_function, division
import os
import sys
import glob
import pandas as pd
import numpy as np
import string
from argparse import ArgumentParser
from subprocess import Popen
from itertools import product
from commonFunctions import get_first_seq, make_motif_file, get_all_sequence_kmers, make_CCWGG_positions_file, \
    find_ccwgg_motifs

PATH_TO_SIGNALALIGN = os.path.abspath("../signalAlign/")
PATH_TO_BINS = PATH_TO_SIGNALALIGN + "/bin/"


def kmer_length_from_model(model_file):
    with open(model_file, "r") as model:
        line = model.readline().split()
        kmer_length = int(line[-1])
        model.close()
    assert kmer_length == 5 or kmer_length == 6
    return kmer_length


def make_positions_file(fasta, degenerate, outfile):
    if degenerate == "adenosine":
        return make_gatc_position_file(fasta, outfile)
    else:
        return make_CCWGG_positions_file(fasta, outfile)


def gatc_kmers(sequence_kmers, kmerlength):
    assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
    # NNNNGATCNNN
    methyl_core = "GITC"
    normal_core = "GATC"
    nucleotides = "ACGT"

    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    labeled_kmers = []

    # add NNNNGA*
    if kmerlength == 6:
        for fourmer in fourmers:
            labeled_kmer = (fourmer + methyl_core)[:kmerlength]
            normal_kmer = (fourmer + normal_core)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNNGA*T and NNNGA*
    for threemer in threemers:
        labeled_kmer = (threemer + methyl_core)[:kmerlength]
        normal_kmer = (threemer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # A*TCNNN
        if kmerlength == 6:
            labeled_kmer = (methyl_core + threemer)[1:]
            normal_kmer = (normal_core + threemer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNGA*TC and NNGA*T
    for twomer in twomers:
        labeled_kmer = (twomer + methyl_core)[:kmerlength]
        normal_kmer = (twomer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # A*TCNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + twomer)[1:]
            normal_kmer = (normal_core + twomer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # NGA*TCN
        if kmerlength == 6:
            labeled_kmer = (twomer[0] + methyl_core + twomer[1])
            normal_kmer = (twomer[0] + normal_core + twomer[1])
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    if kmerlength == 5:
        for onemer in "ACTG":
            labeled_kmer = onemer + methyl_core
            normal_kmer = onemer + normal_core
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = methyl_core + onemer
            normal_kmer = normal_core + onemer
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    return set(labeled_kmers)


def ctag_kmers(sequence_kmers, kmerlength):
    assert kmerlength == 5 or kmerlength == 6, "only works with kmer lengths 5 and 6"
    # NNNCTAGNNNN
    methyl_core = "CTIG"
    normal_core = "CTAG"
    nucleotides = "ACGT"

    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    labeled_kmers = []

    # add A*GNNNN
    if kmerlength == 6:
        for fourmer in fourmers:
            labeled_kmer = (methyl_core + fourmer)[2:]
            normal_kmer = (normal_core + fourmer)[2:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    # add NNNCTA*
    for threemer in threemers:
        if kmerlength == 6:
            labeled_kmer = (threemer + methyl_core)[:kmerlength]
            normal_kmer = (threemer + normal_core)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = (methyl_core + threemer)[1:]
            normal_kmer = (normal_core + threemer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # A*GNNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + threemer)[2:]
            normal_kmer = (normal_core + threemer)[2:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)

    # add NNCTA*G and NNCTA*
    for twomer in twomers:
        labeled_kmer = (twomer + methyl_core)[:kmerlength]
        normal_kmer = (twomer + normal_core)[:kmerlength]
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(labeled_kmer)
        # CTA*GNN
        if kmerlength == 6:
            labeled_kmer = (methyl_core + twomer)[:kmerlength]
            normal_kmer = (normal_core + twomer)[:kmerlength]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
        # TA*GNN
        if kmerlength == 5:
            labeled_kmer = (methyl_core + twomer)[1:]
            normal_kmer = (normal_core + twomer)[1:]
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)

    if kmerlength == 5:
        for onemer in nucleotides:
            labeled_kmer = onemer + methyl_core
            normal_kmer = onemer + normal_core
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
            labeled_kmer = methyl_core + onemer
            normal_kmer = normal_core + onemer
            if normal_kmer in sequence_kmers:
                labeled_kmers.append(labeled_kmer)
    return set(labeled_kmers)


def ccwgg_kmers(sequence_kmers, kmer_length):
    def check_and_add(methyl_kmer):
        normal_kmer = string.translate(methyl_kmer, demethylate)
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(methyl_kmer)

    labeled_kmers = []

    methyl_core1 = "CEAGG"
    methyl_core2 = "CETGG"
    demethylate = string.maketrans("E", "C")

    nucleotides = "ACGT"
    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]
    # NNNNCC*WGGNN

    # NNNNCC*
    if kmer_length == 6:
        for fourmer in fourmers:
            labeled_kmer1 = (fourmer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (fourmer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    # NNNCC*W and NNNCC*
    for threemer in threemers:
        labeled_kmer1 = (threemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (threemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # NNCC*WG and NNCC*W
    for twomer in twomers:
        labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # C*WGGNN
        if kmer_length == 6:
            labeled_kmer1 = (methyl_core1 + twomer)[1:]
            labeled_kmer2 = (methyl_core2 + twomer)[1:]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    for onemer in nucleotides:
        # CC*WGGN and C*WGGN
        labeled_kmer1 = methyl_core1 + onemer
        labeled_kmer2 = methyl_core2 + onemer
        if kmer_length == 6:
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)
        if kmer_length == 5:
            check_and_add(labeled_kmer1[1:])
            check_and_add(labeled_kmer2[1:])
        labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    if kmer_length == 5:
        check_and_add(methyl_core1)
        check_and_add(methyl_core2)

    return set(labeled_kmers)


def ggwcc_kmers(sequence_kmers, kmer_length):
    def check_and_add(methyl_kmer):
        normal_kmer = string.translate(methyl_kmer, demethylate)
        if normal_kmer in sequence_kmers:
            labeled_kmers.append(methyl_kmer)

    labeled_kmers = []

    methyl_core1 = "GGAEC"
    methyl_core2 = "GGTEC"
    demethylate = string.maketrans("E", "C")

    nucleotides = "ACGT"
    fourmers = [''.join(x) for x in product(nucleotides, repeat=4)]
    threemers = [''.join(x) for x in product(nucleotides, repeat=3)]
    twomers = [''.join(x) for x in product(nucleotides, repeat=2)]

    # NNGGWC*CNNN

    # C*CNNNN
    for fourmer in fourmers:
        labeled_kmer1 = (methyl_core1 + fourmer)[3:]
        labeled_kmer2 = (methyl_core2 + fourmer)[3:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # WC*CNNN and C*CNNN
    for threemer in threemers:
        labeled_kmer1 = (methyl_core1 + threemer)[2:] if kmer_length == 6 else (methyl_core1 + threemer)[3:]
        labeled_kmer2 = (methyl_core2 + threemer)[2:] if kmer_length == 6 else (methyl_core2 + threemer)[3:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    # GWC*CNN and WC*CNN
    for twomer in twomers:
        labeled_kmer1 = (methyl_core1 + twomer)[1:] if kmer_length == 6 else (methyl_core1 + twomer)[2:]
        labeled_kmer2 = (methyl_core2 + twomer)[1:] if kmer_length == 6 else (methyl_core2 + twomer)[2:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # NNGGWC*
        if kmer_length == 6:
            labeled_kmer1 = (twomer + methyl_core1)[:kmer_length]
            labeled_kmer2 = (twomer + methyl_core2)[:kmer_length]
            check_and_add(labeled_kmer1)
            check_and_add(labeled_kmer2)

    for onemer in nucleotides:
        # NGGWC* and NGGWC*C
        labeled_kmer1 = (onemer + methyl_core1)[:kmer_length]
        labeled_kmer2 = (onemer + methyl_core2)[:kmer_length]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)
        # GGWC*CN GWC*CN
        labeled_kmer1 = methyl_core1 + onemer if kmer_length == 6 else (methyl_core1 + onemer)[1:]
        labeled_kmer2 = methyl_core2 + onemer if kmer_length == 6 else (methyl_core2 + onemer)[1:]
        check_and_add(labeled_kmer1)
        check_and_add(labeled_kmer2)

    if kmer_length == 5:
        check_and_add(methyl_core1)
        check_and_add(methyl_core2)

    return set(labeled_kmers)


def motif_kmers(core, kmer_length=5, multiplier=5):
    motifs = []
    repeat = kmer_length - len(core)
    if repeat == 0:
        return [core] * multiplier
    else:
        for k in product("ACGT", repeat=repeat):
            fix = ''.join(k)
            motifs.append(fix + core)
            motifs.append(core + fix)
        return motifs * multiplier


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


def make_gatc_or_ccwgg_motif_file(fasta, degenerate, outfile):
    if degenerate == "adenosine":
        motif_finder = find_gatc_motifs
    else:
        motif_finder = find_ccwgg_motifs
    outfile = os.path.abspath(outfile)
    seq = get_first_seq(fasta)
    positions = [x for x in motif_finder(seq)]
    make_motif_file(positions, seq, outfile)
    return outfile


def run_guide_alignment(fasta, pcr_reads, genomic_reads, jobs, positions_file, motif_file, t_model, c_model,
                        outpath, n, degenerate):
    def get_labels():
        if degenerate == "adenosine":
            return ["A", "I"]
        else:
            return ["C", "E"]

    working_path = os.path.abspath(outpath)
    commands = []
    c = PATH_TO_BINS + "runSignalAlign -d={reads} -r={fasta} -T={tModel} -C={cModel} -f=assignments " \
                       "-o={outpath} -p={positions} -q={targetFile} -X={sub} -n={n} -j={jobs}"

    read_sets = [pcr_reads, genomic_reads]
    labels = get_labels()
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
    working_directories = [working_path + "/{}_pcr_".format(degenerate),
                           working_path + "/{}_genomic_".format(degenerate)]

    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)
    assert os.path.exists(t_hdp), "Didn't find template model, looked {}".format(t_hdp)
    assert os.path.exists(c_hdp), "Didn't find complement model, looked {}".format(c_hdp)

    for reads, working_directory in zip(read_sets, working_directories):
        # assemble the command
        command = c.format(reads=reads, fasta=fasta, tModel=t_model, cModel=c_model, outpath=working_directory,
                           positions=positions_file, targetFile=motif_file, n=n, jobs=int(jobs/2), tHdp=t_hdp,
                           cHdp=c_hdp, degenerate=degenerate)
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


def train_model_transitions(fasta, pcr_reads, genomic_reads, degenerate, jobs, positions_file, iterations, batch_size,
                            outpath, t_model, c_model, stateMachine="threeState", t_hdp=None, c_hdp=None):
    working_path = os.path.abspath(outpath) + "/"
    model_directory = working_path + "{}_".format(stateMachine)
    assert os.path.exists(t_model), "Didn't find template model, looked {}".format(t_model)
    assert os.path.exists(c_model), "Didn't find complement model, looked {}".format(c_model)

    if degenerate == "adenosine":
        methyl_char = "I"
    else:
        methyl_char = "E"

    os.chdir(PATH_TO_BINS)
    c = "trainModels -d={pcr} -d={genomic} -X=C -X={methylChar} -r={fasta} -i={iter} -a={batch} --transitions " \
        "-smt={smt} -T={tModel} -C={cModel} -j={jobs} -x=adenosine -p={positions} -o={out} " \
        "".format(pcr=pcr_reads, genomic=genomic_reads, fasta=fasta, iter=iterations, batch=batch_size,
                  smt=stateMachine, tModel=t_model, cModel=c_model, jobs=jobs, methylChar=methyl_char,
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


def make_build_alignment(assignments, degenerate, kmer_length, ref_fasta, n_canonical_assignments,
                         n_methyl_assignments, outfile, threshold):
    def write_kmers(max_assignments, kmer_list):
        for strand in ["t", "c"]:
            by_stand = assignments.ix[(assignments['strand'] == strand) & (assignments['prob'] >= threshold)]
            for k in kmer_list:
                kmer_assignments = by_stand.ix[by_stand['kmer'] == k]
                if kmer_assignments.empty:
                    print("missing kmer {}, continuing".format(k))
                    continue
                kmer_assignments = kmer_assignments.sort_values(['prob'], ascending=0)
                n = 0
                for _, r in kmer_assignments.iterrows():
                    fH.write(
                        entry_line.format(strand=r['strand'], kmer=r['kmer'], event=r['level_mean'], prob=r['prob']))
                    n += 1
                    if n >= max_assignments:
                        break
                if n < max_assignments:
                    print("WARNING didn't find {max} requested assignments for {kmer} only found {found}"
                          "".format(max=max_assignments, kmer=k, found=n))
    seq = get_first_seq(ref_fasta)
    sequence_kmers = get_all_sequence_kmers(seq, kmer_length).keys()
    if degenerate == "adenosine":
        methyl_kmers = list(gatc_kmers(sequence_kmers=sequence_kmers, kmerlength=kmer_length))
        methyl_kmers += list(ctag_kmers(sequence_kmers=sequence_kmers, kmerlength=kmer_length))
    else:
        #print("C-methylation not fixed yet", file=sys.stderr)
        #sys.exit(1)
        methyl_kmers = list(ccwgg_kmers(sequence_kmers=sequence_kmers, kmer_length=kmer_length))
        methyl_kmers += list(ggwcc_kmers(sequence_kmers=sequence_kmers, kmer_length=kmer_length))
    fH = open(outfile, "w")
    entry_line = "blank\t0\tblank\tblank\t{strand}\t0\t0.0\t0.0\t0.0\t{kmer}\t0.0\t0.0\t{prob}\t{event}\t0.0\n"
    write_kmers(n_canonical_assignments, sequence_kmers)
    write_kmers(n_methyl_assignments, methyl_kmers)
    fH.close()
    return outfile


def build_hdp(build_alignment_path, template_model, complement_model, outpath, samples=15000):
    working_path = os.path.abspath(outpath) + "/"
    build_alignment = os.path.abspath(build_alignment_path)
    t_model = os.path.abspath(template_model)
    c_model = os.path.abspath(complement_model)
    outpath = os.path.abspath(outpath) + "/"
    hdp_pipeline_dir = outpath + "hdpPipeline/"
    os.makedirs(hdp_pipeline_dir)
    os.chdir(PATH_TO_BINS)
    c = "hdp_pipeline --build_alignment={build} -tM={tModel} -cM={cModel} -Ba=1 -Bb=1 -Ma=1 -Mb=1 -La=1 -Lb=1 " \
        "-s={samples} --verbose --grid_start=50 --grid_end=140 --grid_length=1800 --verbose -o={out} " \
        "--hdp_type=ecoli".format(build=build_alignment, tModel=t_model, cModel=c_model, samples=samples,
                                  out=hdp_pipeline_dir)
    c = PATH_TO_BINS + c
    os.system(c)
    os.chdir(working_path)
    return [hdp_pipeline_dir + "template.singleLevelPriorEcoli.nhdp",
            hdp_pipeline_dir + "complement.singleLevelPriorEcoli.nhdp"]


def main(args):
    def parse_args():
        parser = ArgumentParser(description=__doc__)
        parser.add_argument("-r", action="store", dest="reference", required=True)
        parser.add_argument("-pcr", action="store", dest="pcr_reads", required=True)
        parser.add_argument("-gen", action="store", dest="genomic_reads", required=True)
        parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm', required=True, type=str)
        parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm', required=True, type=str)
        parser.add_argument("-o", action="store", dest="outpath", required=True)
        parser.add_argument("-x", action="store", dest="degenerate", required=True)
        parser.add_argument("-j", action="store", dest="jobs", required=False, default=4, type=int)
        parser.add_argument("-i", action="store", dest="iterations", required=False, type=int, default=20)
        parser.add_argument("-a", action="store", dest="batch", required=False, type=int, default=15000)
        parser.add_argument("-s", action="store", dest="assignments", required=False, type=int, default=30)
        parser.add_argument("-c", action="store", dest="methyl_assignments", required=False, type=int, default=200)
        parser.add_argument("-n", action="store", dest="n_aligns", required=False, type=int, default=1000)
        parser.add_argument("-e", action="store", dest="n_test_alns", required=False, type=int, default=1000)
        parser.add_argument("-t", action="store", dest="assignment_threshold", required=False, type=float, default=0.8)
        parser.add_argument("-g", action="store", dest="samples", required=False, type=int, default=15000)
        args = parser.parse_args()
        return args

    command_line = " ".join(sys.argv[:])
    print("Command Line: {cmdLine}\n".format(cmdLine=command_line), file=sys.stderr)

    args = parse_args()

    # make the positions and motif file
    working_path = os.path.abspath(args.outpath)
    positions_file = make_positions_file(fasta=args.reference,
                                         degenerate=args.degenerate,
                                         outfile=working_path + "/{}_positions.positions".format(args.degenerate))

    # make the motif file
    motif_file = make_gatc_or_ccwgg_motif_file(fasta=args.reference,
                                               degenerate=args.degenerate,
                                               outfile=working_path + "/{}_target.target".format(args.degenerate))

    # train the transitions
    models = train_model_transitions(fasta=os.path.abspath(args.reference),
                                     pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                     genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                     degenerate=args.degenerate,
                                     jobs=args.jobs,
                                     positions_file=positions_file,
                                     iterations=args.iterations,
                                     batch_size=args.batch,
                                     outpath=working_path,
                                     t_model=os.path.abspath(args.in_T_Hmm),
                                     c_model=os.path.abspath(args.in_C_Hmm))
    # do the initial alignments
    assignment_dirs = run_guide_alignment(fasta=os.path.abspath(args.reference),
                                          pcr_reads=os.path.abspath(args.pcr_reads) + "/",
                                          genomic_reads=os.path.abspath(args.genomic_reads) + "/",
                                          jobs=args.jobs,
                                          positions_file=positions_file,
                                          motif_file=motif_file,
                                          n=args.n_aligns,
                                          degenerate=args.degenerate,
                                          t_model=models[0],
                                          c_model=models[1],
                                          outpath=working_path)
    assert kmer_length_from_model(models[0]) == kmer_length_from_model(models[1]), "Models had different kmer lengths"
    # concatenate the assignments into table
    master = make_master_assignment_table(assignment_dirs)
    build_alignment = make_build_alignment(assignments=master,
                                           degenerate=args.degenerate,
                                           kmer_length=kmer_length_from_model(models[0]),
                                           ref_fasta=os.path.abspath(args.reference),
                                           n_canonical_assignments=args.assignments,
                                           n_methyl_assignments=args.methyl_assignments,
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
                                         degenerate=args.degenerate,
                                         jobs=args.jobs,
                                         positions_file=positions_file,
                                         iterations=args.iterations,
                                         batch_size=args.batch,
                                         outpath=working_path,
                                         stateMachine="threeStateHdp",
                                         t_hdp=hdps[0],
                                         c_hdp=hdps[1],
                                         t_model=os.path.abspath(args.in_T_Hmm),
                                         c_model=os.path.abspath(args.in_C_Hmm))
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
                                   degenerate=args.degenerate,
                                   t_hdp=hdp_models[2],
                                   c_hdp=hdp_models[3])
    # run the control experiment
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
                                   degenerate="variant",
                                   t_hdp=hdp_models[2],
                                   c_hdp=hdp_models[3])

if __name__ == "__main__":
    sys.exit(main(sys.argv))







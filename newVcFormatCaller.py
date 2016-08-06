#!/usr/bin/env python
from __future__ import division, print_function
import sys
import pandas as pd
import numpy as np
from collections import Counter

def parse_vc_alignment_file(filepath):
    table = pd.read_table(filepath,
                          usecols=(0, 1, 2, 3, 4, 5, 6),
                          names=['event_idx', 'ref_pos', 'base', 'prob', 'strand', 'forward', 'read_label'],
                          dtype={
                              'event_idx': np.int64,
                              'ref_pos': np.int64,
                              'base': np.str,
                              'prob': np.float64,
                              'strand': np.str,
                              'forward': np.str,
                              'read_label': np.str},
                          header=None)
    return table


def normalize_probs(probs):
    total = sum(probs.values())
    for base in probs:
        probs[base] /= total


def degenerate_probs_map(degenerate):
    if degenerate == "twoWay":
        return {'C': 0.0, 'E': 0.0}
    elif degenerate == "threeWay":
        return {'C': 0.0, 'E': 0.0, "O": 0.0}
    else:
        return {"A": 0.0, "C": 0.0, "G": 0.0, "T": 0.0}


def zymo_analysis(alignment_filepath, correct_label, incorrect_labels, degenerate, threshold, read_score_threshold,
                  read_accuracy_outfile, site_accuracy_outfile):
    def read_score(read_df):
        total_prob = sum(read_df['prob'])
        return 100 * total_prob / len(read_df['prob'])

    def call_reads(calls_df):
        # group the calls DataFrame by read
        for read, read_df in calls_df.groupby('read_label'):
            # call each site for that read
            score = read_score(read_df=read_df)
            if score < read_score_threshold:
                continue
            site_calls = {}
            for site, probs in call_sites(read_df):
                site_calls[site] = max(probs, key=probs.get)
            # read: read_label, site_calls {site: call}, score
            yield read, site_calls, score

    def call_sites(calls_df):
        # group the aligned pairs by the reference position they are reporting on
        for site, site_df in calls_df.groupby('ref_pos'):
            # probs is a map base : probabilty
            probs = degenerate_probs_map(degenerate=degenerate)
            # only take pairs above a threshold posterior probability
            thresholded_df = site_df.ix[site_df['prob'] >= threshold]
            if thresholded_df.empty:
                print("skipping", site, file=sys.stderr)
                continue
            for _, row in thresholded_df.iterrows():
                posterior_prob = row['prob']
                called_base = row['base']
                if posterior_prob < threshold:
                    print("Error sorting by posterior probs")
                    sys.exit(1)
                probs[called_base] += posterior_prob
            normalize_probs(probs=probs)
            yield site, probs

    def read_accuracy():
        fH = open(read_accuracy_outfile, "w")
        for strand in strands:
            by_strand = aln.ix[aln['strand'] == strand]
            for read, site_calls, score in call_reads(by_strand):
                correct = site_calls.values().count(correct_label)
                incorrect = 0
                for incorrect_label in incorrect_labels:
                    incorrect += site_calls.values().count(incorrect_label)
                total = len(site_calls)
                accuracy = correct / total
                print("{accuracy}\t{score}\t{strand}\t{read_label}" \
                      "".format(accuracy=accuracy, read_label=read, score=score, strand=strand), file=fH)
        fH.close()

    def site_accuracy():
        fH = open(site_accuracy_outfile, "w")
        for strand in strands:
            by_strand = aln.ix[aln['strand'] == strand]
            correct_site_counts = Counter()
            incorrect_site_counts = Counter()
            for _, site_calls, _ in call_reads(by_strand):
                for site in site_calls:
                    called_base = site_calls[site]
                    if called_base == correct_label:
                        correct_site_counts[site] += 1
                    else:
                        incorrect_site_counts[site] += 1
            all_sites = set(correct_site_counts.keys() + incorrect_site_counts.keys())
            for site in all_sites:
                try:
                    accuracy = correct_site_counts[site] / (correct_site_counts[site] + incorrect_site_counts[site])
                except KeyError:
                    accuracy = 0.0
                print("{accuracy}\t{site}\t{strand}\n".format(accuracy=accuracy, site=site, strand=strand),
                      file=fH)
        fH.close()
        return

    strands = ["t", "c"]
    aln = parse_vc_alignment_file(alignment_filepath)
    site_accuracy()
    read_accuracy()
    return

f = sys.argv[1]
zymo_analysis(f, "C", ["E"], "twoWay", 0.0, 0.0, sys.argv[2], sys.argv[3])



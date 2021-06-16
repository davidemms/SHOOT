#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import random
import subprocess

import fasta_writer


def run(din):
    # 1. Create diamond profile for orthogroup assignment
    create_og_lookup(din)


def get_orthogroups(clustersFilename):
    """
    Get OGs as a list of sets of gene id-pair strings
    Args:
        clustersFilename - filename
    """
    predictedOGs = []
    nOGsString = ""
    qContainsProfiles = False
    with open(clustersFilename, 'rb') as clusterFile:
        header = True
        og = set()
        for line in clusterFile:
            if header:
                if line.count("begin"):
                    header = False
            else:
                if line.find(")") != -1:
                    break
                if line[-2] == "$":
                    line = line[:-3]
                if line[0] == " ":
                    # continuation of group
                    x = line.split()
                    y = [x_ for x_ in x if not x_.startswith('Prof')]
                    og = og.union(y)
                else:
                    # new OG
                    if len(og) != 0:
                        predictedOGs.append(og)
                    nOGsString, line = line.split(" ", 1)
                    x = line.split()
                    y = [x_ for x_ in x if not x_.startswith('Prof')]
                    if len(x) != len(y):
                        qContainsProfiles = True
                    og = set(y)
        if len(og) > 0:
            predictedOGs.append(og)
    if not qContainsProfiles:
        assert(len(predictedOGs) == int(nOGsString) + 1)
    return predictedOGs


def create_og_lookup(din, n_for_profile=5):
    wd = din + "WorkingDirectory/"
    fn_fasta = din + "diamond_profile_sequences.fa"
    fn_diamond_db = fn_fasta + ".db"
    ogs = get_orthogroups(wd + "clusters_OrthoFinder_I1.5.txt_id_pairs.txt")
    fw = fasta_writer.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    for iog, og in enumerate(ogs):
        n = len(og)
        if n < 4:
            og_id = "x%07d_" % iog      # indicates no tree
        else:
            og_id = "%07d_" % iog
        s = random.sample(og, min(n_for_profile, n))
        seq_write.extend(s)
        for ss in s:
            seq_convert[ss] = og_id + ss
    fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    subprocess.call(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db])
    

if __name__ == "__main__":
    din = sys.argv[1]
    if not din.endswith("/"):
        din += "/"
    run(din)

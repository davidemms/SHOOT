#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob
import os
import sys
import random
import subprocess

import fasta_writer
import sample_genes


def run(din):
    # 1. Create diamond profile for orthogroup assignment
    create_profiles_database(din)


def get_orthogroups(clustersFilename):
    """
    Read orthogroups from file
    Args:
        clustersFilename - filename
    Returns:
        ogs - list of sets of strings (isp_iseq)
    """
    predictedOGs = []
    nOGsString = ""
    qContainsProfiles = False
    with open(clustersFilename, 'r') as clusterFile:
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


def sample_random(og, n_max):
    """
    Sample min(n_max, |og|) genes randomly from clade
    Args:
        og - set of strings
        n_max - max number of genes to sample
    Returns:
        genes - list of genes
    """
    return random.sample(og, min(n_max, len(og)))


def create_profiles_database(din, q_kmeans = True, n_for_profile=5, q_ids=True):
    """
    Create a fasta file with profile genes from each orthogroup
    Args:
        din - Input OrthoFinder results directory
        n_for_profile - The number of genes to use from each orthogroup, when available
        q_kmeans - Should kmeans be used instead of random sampling
    Notes:
    If the trees have been split into subtres then profiles will be created for 
    the subtrees instead.
    A file "diamond_profile_sequences.fa.db.dmnd" is created in the OrthoFinder
    directory suitable for use by SHOOT. Each sequence is named <OG>_isp_iseq.
    """
    wd = din + "WorkingDirectory/"
    pat_super = din + "Gene_Trees/subtrees/super/OG%07d.super.tre"
    pat_sub_msa_glob = din + "Gene_Trees/subtrees/msa_sub/OG%07d.*.fa"
    fn_fasta = din + "diamond_profile_sequences.new.fa"
    fn_diamond_db = fn_fasta + ".db"
    ogs = get_orthogroups(wd + "clusters_OrthoFinder_I1.5.txt_id_pairs.txt")
    fw = fasta_writer.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    print("WARNING: Check all gene names, can't start with '__'")
    # If there are subtrees then we need to convert their IDs in the profile file
    # back to internal IDs
    if q_ids:
        import ofids
        ids = ofids.OrthoFinderIDs(wd).SequenceDict()
        ids_rev = {v:k for k,v in ids.items()}
    for iog, og in enumerate(ogs):
        if iog % 10 == 0:
            print(iog)
        n = len(og)
        og_id = "%07d" % iog
        # if n < 4:
        #     og_id = "x%07d" % iog      # indicates no tree
        # else:
        #     og_id = "%07d" % iog
        q_subtrees = os.path.exists(pat_super % iog)
        if q_subtrees:
            fns_msa = list(glob.glob(pat_sub_msa_glob % iog))
        else:
            fns_msa = [wd + "Alignments_ids/OG%07d.fa" % iog, ]
        for fn in fns_msa:
            i_part = os.path.basename(fn).rsplit(".", 2)[1]
            if q_kmeans:
                if q_subtrees:
                    # MSA needs to be modified
                    fn_temp = "/tmp/shoot_db_create" + os.path.basename(fn)
                    fw_temp = fasta_writer.FastaWriter(fn)
                    fw_temp.WriteSeqsToFasta([g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")],
                                            fn_temp)
                    fn = fn_temp
                # Don't trim as OrthoFinder has already trimmed by default
                s = sample_genes.select_from_aligned(fn, n_for_profile, q_trim=False)
                if q_subtrees:
                    os.remove(fn_temp)
            else:
                if q_subtrees:
                    fw_temp = fasta_writer.FastaWriter(fn)
                    og = [g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")]
                s = sample_random(og, n_for_profile)
            if (q_ids and q_subtrees):
                s = [ids_rev[ss] for ss in s]
            if q_subtrees:
                og_id_full = og_id + "." + i_part
            else:
                og_id_full = og_id
            seq_write.extend(s)
            for ss in s:
                seq_convert[ss] = og_id_full + "_" + ss
        # re-write it for each OG so I can check it's going ok
        # print(seq_write)
        # print(list(fw.SeqLists.keys())[:10])
        fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    subprocess.call(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db])


if __name__ == "__main__":
    din = sys.argv[1]
    if not din.endswith("/"):
        din += "/"
    run(din)

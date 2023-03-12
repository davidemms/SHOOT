#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Steps to fix up datasets:
1. Gene IDs, detailed elsewhere, but should put it here
2. EPA requires bifurcating trees. See bifurcating_trees.py
3. IQTREE: must eliminate sequences that are all gaps.
"""

import glob
import os
import sys
import random
import subprocess
import string

from shoot import fasta_writer, ofids, sample_genes


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


def create_mmseqs_database(din, subtrees_dir="Gene_Trees/subtrees_2000"):
    """
    Create a stockholm file of all the MSA
    Args:
        din - Input OrthoFinder results directory
        subtrees_dir - directory containing the split trees to use
    Notes:
    If the trees have been split into subtrees then profiles will be created for
    the subtrees instead.
    """
    wd = din + "WorkingDirectory/"
    subtrees_label = os.path.split(subtrees_dir)[1]
    pat_super = din + subtrees_dir + "/super/OG%07d.super.tre"
    pat_sub_msa_glob = din + subtrees_dir + "/msa_sub/OG%07d.*.fa"
    fn_stockholm = din + "all_msa.%s.sto" % subtrees_label
    fn_mmseqs_db = fn_stockholm + ".db"
    ogs = get_orthogroups(wd + "clusters_OrthoFinder_I1.5.txt_id_pairs.txt")
    seq_write = []
    seq_convert = dict()
    # delete the existing file
    with open(fn_stockholm, 'w') as outfile:
        pass
    for iog, og in enumerate(ogs):
        # if iog < 20 or iog > 29:
        #     continue
        if iog % 10000 == 0:
            print(iog)
        n = len(og)
        og_id = "%07d" % iog
        q_subtrees = os.path.exists(pat_super % iog)
        fn_msa_unsplit = din + "MultipleSequenceAlignments/OG%07d.fa" % iog
        fn_seq_unsplit = din + "Orthogroup_Sequences/OG%07d.fa" % iog
        if q_subtrees:
            fns_msa = list(glob.glob(pat_sub_msa_glob % iog))
        elif os.path.exists(fn_msa_unsplit):
            fns_msa = [fn_msa_unsplit, ]
        else:
            fns_msa = [fn_seq_unsplit, ]
        for fn in fns_msa:
            i_part = os.path.basename(fn).rsplit(".", 2)[1] if q_subtrees else None
            fw_temp = fasta_writer.FastaWriter(fn)
            if q_subtrees:
                og = [g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")]
                og_id_full = og_id + "." + i_part
            else:
                og = None
                og_id_full = og_id
            fw_temp.AppendToStockholm(og_id_full, fn_stockholm, og)
    # subprocess.call(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db])


def create_profiles_database(din, selection="kmeans", min_for_profile=20, q_ids=True, divide=10,
                            subtrees_dir=""):
    """
    Create a fasta file with profile genes from each orthogroup
    Args:
        din - Input OrthoFinder results directory
        selection - "kmeans|random|all"  - method to use for sampling genes from orthogroups
        n_for_profile - The number of genes to use from each orthogroup, when available
        q_ids - Convert subtrees (with user gene accessions) back to IDs for profiles database (should always be true I think)
    Notes:
    If the trees have been split into subtrees then profiles will be created for
    the subtrees instead.
    A file "diamond_profile_sequences.fa.db.dmnd" is created in the OrthoFinder
    directory suitable for use by SHOOT. Each sequence is named <OG>_isp_iseq.
    """
    if selection not in ("kmeans", "random", "all"):
        raise RuntimeError("selection method '%s' is not defined" % selection)
    wd = din + "WorkingDirectory/"
    if subtrees_dir:
        subtrees_label = "." + os.path.split(subtrees_dir)[1]
        pat_super = din + subtrees_dir + "/super/OG%07d.super.tre"
        pat_sub_msa_glob = din + subtrees_dir + "/msa_sub/OG%07d.*.fa"
    else:
        subtrees_label = ""
    ogs = get_orthogroups(wd + "clusters_OrthoFinder_I1.5.txt_id_pairs.txt")
    fw = fasta_writer.FastaWriter(wd + "Species*fa", qGlob=True)
    seq_write = []
    seq_convert = dict()
    # print("WARNING: Check all gene names, can't start with '__'")
    # If there are subtrees then we need to convert their IDs in the profile file
    # back to internal IDs
    if q_ids:
        ids = ofids.OrthoFinderIDs(wd).Spec_SeqDict()
        ids_rev = {v: k for k, v in ids.items()}
    nToDo = len(ogs)
    for iog, og in enumerate(ogs):
        if iog >= 0 and divmod(iog, 10 if nToDo <= 200 else 100 if nToDo <= 2000 else 1000)[1] == 0:
            print("Done %d of %d" % (iog, nToDo))
        og_id = "%07d" % iog
        q_subtrees = subtrees_dir and os.path.exists(pat_super % iog)
        fn_msa = wd + "Alignments_ids/OG%07d.fa" % iog
        if q_subtrees: 
            print("Subtrees: %d" % iog)
            fns_msa = list(glob.glob(pat_sub_msa_glob % iog))
        elif os.path.exists(fn_msa):
            fns_msa = [fn_msa, ]
        else:
            fns_msa = [wd + "Sequences_ids/OG%07d.fa" % iog, ]
        for fn in fns_msa:
            if not os.path.exists(fn):
                print("File does not exist, skipping: %s" % fn)
                continue
            i_part = os.path.basename(fn).rsplit(".", 2)[1] if q_subtrees else None
            fw_temp = fasta_writer.FastaWriter(fn)
            n_in_og = len(fw_temp.SeqLists)
            n_for_profile = n_in_og // divide
            if n_for_profile < min_for_profile:
                n_for_profile = min_for_profile
            if selection == "kmeans" and len(fw_temp.SeqLists) > n_for_profile:
                if q_subtrees:
                    # MSA needs to be modified
                    letters = string.ascii_lowercase
                    fn_temp = "/tmp/shoot_db_create" + "".join(random.choice(letters) for i in range(6)) + os.path.basename(fn)
                    fw_temp.WriteSeqsToFasta([g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")],
                                            fn_temp)
                    fn = fn_temp
                # Don't trim as OrthoFinder has already trimmed by default
                s = sample_genes.select_from_aligned(fn, n_for_profile, q_trim=False)
                if q_subtrees:
                    os.remove(fn_temp)
            elif selection == "kmeans" or selection == "random":
                if q_subtrees:
                    fw_temp = fasta_writer.FastaWriter(fn)
                    og = [g for g in fw_temp.SeqLists if not g.startswith("SHOOTOUTGROUP_")]
                s = sample_random(og, n_for_profile)
            else:
                s = [g for g in fasta_writer.FastaWriter(fn).SeqLists.keys() if not g.startswith("SHOOTOUTGROUP_")]
            if q_ids and q_subtrees:
                s = [ids_rev[ss] for ss in s]
            if q_subtrees:
                og_id_full = og_id + "." + i_part
            else:
                og_id_full = og_id
            seq_write.extend(s)
            for ss in s:
                seq_convert[ss] = og_id_full + "_" + ss
    if selection == "all":
        fn_fasta = din + "profile_sequences%s.all.fa" % subtrees_label
    else:
        fn_fasta = din + "profile_sequences%s.d%d_min%d_%s.fa" % (subtrees_label, divide, min_for_profile, selection)
    fw.WriteSeqsToFasta_withNewAccessions(seq_write, fn_fasta, seq_convert)
    fn_diamond_db = fn_fasta + ".db"
    subprocess.call(["diamond", "makedb", "--in", fn_fasta, "-d", fn_diamond_db])
    return fn_diamond_db


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Args: input_directory [full,profiles]")
        sys.exit() 
    din = sys.argv[1]
    full_or_profiles = sys.argv[2]
    if full_or_profiles not in ["full", "profiles", "mmseqs"]:
        print("Args: input_directory [full,profiles,mmseqs]")
        sys.exit() 
    if not din.endswith("/"):
        din += "/"
    if full_or_profiles == "full":
        create_profiles_database(din, selection="all", subtrees_dir="")
    elif full_or_profiles == "profiles":
        create_profiles_database(din, 
                                selection="kmeans",
                                min_for_profile=20, 
                                q_ids=True, 
                                divide= 10, 
                                subtrees_dir="")
    elif full_or_profiles == "mmseqs":
        create_mmseqs_database(din)



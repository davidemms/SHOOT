#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import sys
import argparse

import og_assigner
import tree_grafter
import msa_grafter
import quartets_pairwise_align


def main(d_db, infn, q_msa, q_print=False):
    query_name, query_renamed, infn = rename_gene(infn)
    og_assign = og_assigner.OGAssignDIAMOND(d_db)
    og_part = og_assign.assign(infn)
    if og_part is not None:
        print("Gene assigned to: OG%s" % og_part)
    else:
        print("No homologs found for gene in this database")
        return ""

    warn_str = ""
    if q_msa:
        # do a tree using an MSA
        graft = msa_grafter.MSAGrafter(d_db)
        fn_tree, warn_str = graft.add_gene(og_part, infn, query_name, query_renamed)
    else:
        quart = quartets_pairwise_align.PairwiseAlignQuartets(d_db, og_part, infn)
        search = tree_grafter.TreeGrafter(quart, d_db)
        # This doesn't feel right, iog is fixed in the constructor of quart and hence
        # in the constructor of search, and yet is passed as a variable here.
        # search.place_gene(iog)              
        search.place_gene() 

    print("Tree: %s" % fn_tree)  
    if warn_str != "":
        print("WARNING: " + warn_str) 
    if q_print:
        with open(fn_tree, 'r') as infile:
            print(next(infile).rstrip())   # remove any trailing newline characters
    return fn_tree        


def rename_gene(infn):
    preface_str = "YdN3Z" # this will not occur in the shoot database
    with open(infn, 'r') as infile:
        acc = next(infile)
        if acc.startswith(">"):
            acc = acc[1:].rstrip()
            acc_new = preface_str + acc
        else:
            acc = "Query"
            acc_new = preface_str + acc
        lines = [l for l in infile] # python incorrectly complains about readlines
    fn_new = infn + ".rn.fa"
    with open(fn_new, 'w') as outfile:
        outfile.write(">" + acc_new + "\n")
        for l in lines:
            outfile.write(l)
    return acc, acc_new, fn_new


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    parser.add_argument("-m", "--msa", action="store_true", help= "Use an MSA tree")
    parser.add_argument("-p", "--print_tree", action="store_true", help= "Print tree as final line")
    args = parser.parse_args()
    main(args.db, args.infile, args.msa, args.print_tree)
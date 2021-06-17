#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import sys
import argparse

import og_assigner
import tree_grafter
import msa_grafter
import quartets_pairwise_align


def main(d_db, infn, q_msa, q_print=False):
    og_assign = og_assigner.OGAssignDIAMOND(d_db)
    iog, q_tree = og_assign.assign(infn)
    if iog != -1:
        print("Gene assigned to: OG%07d" % iog)
    # if not q_tree:
    #     print("No tree, analysis complete")
    #     return
    # else:

    warn_str = ""
    if q_msa:
        # do a tree using an MSA
        graft = msa_grafter.MSAGrafter(d_db)
        fn_tree, warn_str = graft.add_gene(iog, infn)
    else:
        quart = quartets_pairwise_align.PairwiseAlignQuartets(d_db, iog, infn)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    parser.add_argument("-m", "--msa", action="store_true", help= "Use an MSA tree")
    parser.add_argument("-p", "--print_tree", action="store_true", help= "Print tree as final line")
    args = parser.parse_args()
    main(args.db, args.infile, args.msa, args.print_tree)
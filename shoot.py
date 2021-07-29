#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import sys
import re
import argparse

import og_assigner
import tree_grafter
import msa_grafter
import msa_grafter_epa
import quartets_pairwise_align

import ete3

def main(d_db, infn, q_msa, nU, nL, tree_method, q_print=False):
    """
    Run SHOOT
    Args:
        d_db - Input directory containing SHOOT database
        infn - Query FASTA filename
        q_msa - Use MSA method for tree inference
        nU - Upper limit for tree, unless nL is not None
        nL - Exceed nU if alternative is a tree < nL
        tree_method - method to use for tree inference
    Returns:
        fn_tree - Filename for tree or None
    """
    # rapid check for FASTA format
    ok = True
    with open(infn, 'r') as infile:
        ok = bool(re.match("^>.", next(infile))) 
        ok = ok and bool(re.match("^[a-zA-Z]+$", next(infile))) 
    if not ok:
        print("ERROR: Input file should be FASTA format")
        return

    og_assign = og_assigner.OGAssignDIAMOND(d_db)
    og_part = og_assign.assign(infn)
    if og_part is not None:
        print("Gene assigned to: OG%s" % og_part)
    else:
        print("No homologs found for gene in this database")
        return 

    warn_str = ""
    if q_msa:
        # do a tree using an MSA
        if "iqtree" == tree_method:
            graft = msa_grafter.MSAGrafter(d_db)
        elif "epa" == tree_method:
            graft = msa_grafter_epa.MSAGrafter_EPA(d_db)
        else:
            print("ERROR: %s method has not been implemented" % tree_method)
            return
        fn_tree, query_gene, warn_str = graft.add_gene(og_part, infn)
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
        
    q_need_to_print = q_print
    if nU is not None:
        # if only nL is specified alone that has no effect
        t = ete3.Tree(fn_tree)
        if len(t) > nU:
            node = t & query_gene
            while len(node) < nU:
                node_prev = node
                n_taxa_prev = len(node)
                node = node.up
            # now there are more than nU genes in this tree, step down one
            # unless it is fewer than nL
            node = node_prev if (nL is None or n_taxa_prev >= nL) else node
            nwk_str = node.write()
            with open(fn_tree, 'w') as outfile:
                outfile.write(nwk_str)
            if q_need_to_print:
                print(nwk_str)
                q_need_to_print = False

    if q_need_to_print:
        with open(fn_tree, 'r') as infile:
            print(next(infile).rstrip())   # remove any trailing newline characters
    return fn_tree        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    parser.add_argument("-m", "--msa", action="store_true", help= "Use an MSA tree")
    parser.add_argument("-u", "--upper", type=int, help= "Upper limit for tree, unless -l")
    parser.add_argument("-l", "--lower", type=int, help= "Exceed -u if alternative is < -l")
    parser.add_argument("-p", "--print_tree", action="store_true", help= "Print tree as final line")
    parser.add_argument("-t", "--tree_method", default="epa", choices={"epa", "iqtree"})
    args = parser.parse_args()
    main(args.db, args.infile, args.msa, args.upper, args.lower, args.tree_method, args.print_tree)
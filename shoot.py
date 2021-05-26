import sys
import argparse

import og_assigner
import tree_searcher
import quartets_pairwise_align


def main(d_db, infn):
    og_assign = og_assigner.OGAssignDIAMOND(d_db)
    iog, q_tree = og_assign.assign(infn)
    if not q_tree:
        print("No tree, analysis complete")
        return

    quart = quartets_pairwise_align.PairwiseAlignQuartets(d_db, iog, infn)
    search = tree_searcher.TreeSearcher(quart, d_db)

    # This doesn't feel right, iog is fixed in the constructor of quart and hence
    # in the constructor of search, and yet is passed as a variable here.
    search.place_gene(iog)              

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    args = parser.parse_args()
    main(args.db, args.infile)
import sys
import argparse

import og_assigner


def main(d_db, infn):
    og_assign = og_assigner.OG_Assign_D(d_db)
    iog, q_tree = og_assign.assign(infn)
    if not q_tree:
        print("No tree, analysis complete")
        return
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    args = parser.parse_args()
    main(args.db, args.infile)
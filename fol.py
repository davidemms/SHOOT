"""
Places a gene in the Forest of Life -- prototype:
"""
import gzip
import csv
import argparse
import subprocess
from collections import Counter


def run(fn_query, fn_out_base, d_ox_db):
    if not d_ox_db.endswith("/"):
        d_ox_db += "/"
    iog, q_tree = assign_to_orthogroup(fn_query, fn_out_base, d_ox_db)
    if not q_tree:
        print("No tree, exiting")
        return

def blast(fn_query, fn_db, w=3, T=11, X=2):
    """
    Args:
        w - word length
        T - threshold for high-scoring words
        X - max score drop-off
    """

    # get all query words
    
    # get all high scoring words for this list (APM/BLOSUM)


    # find locations of these in database


    # extend in both directions until score drops


    # sort all extended hits by alignment score



def assign_to_orthogroup(fn_query, fn_out_base, d_ox_db):
    fn_results = run_diamond(fn_query, fn_out_base, d_ox_db)
    iog, q_tree = og_from_diamond_results(fn_results)
    print("Belongs of OG%07d" % iog)
    return iog, q_tree


def run_diamond(fn_query, fn_out_base, d_ox_db):
    fn_db = d_ox_db + "diamond_profile_sequences.fa.db.dmnd"
    fn_og_results_out = fn_out_base + ".ogs.txt"
    subprocess.call(["diamond", "blastp", "-d", fn_db, "-q", fn_query, "-o", fn_og_results_out, "--quiet", "-e", "0.001", "--compress", "1"])
    return fn_og_results_out + ".gz"


def og_from_diamond_results(fn_og_results_out):
    """
    Returns:
        iog        : int
        q_has_tree : bool
    """
    ogs = []
    scores = []
    with gzip.open(fn_og_results_out, 'rt') as infile:
        reader = csv.reader(infile, delimiter="\t")
        for l in reader:
            ogs.append(l[1].split("_")[0])
            scores.append(float(l[-2]))
    sortedTuples = sorted(zip(scores, ogs))
    scores = [i for i, j in sortedTuples]
    ogs = [j for i, j in sortedTuples]
    c = Counter(ogs[:5])
    a, b = c.most_common(2)
    if a[1] > 2:
        og = a[0]
    elif b[1] == 2:
        # we have a tie, get the highest scoring of the two
        og = next(r for r in ogs if (r == a[0] or r == b[0]))
    elif a[1] == 2:
        # a is still the winner
        og = a[0]
    else:
        # all scored 1, get the highest scoring
        og = next(r for r in ogs if (r == a[0] or r == b[0]))
    if og.startswith("x"):
        return int(og[:1]), False
    else:
        return int(og), True


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("query", help="input query FASTA file")
    parser.add_argument("db", help="OX database directory")
    parser.add_argument("out", help="Output filename")
    # OX database directory should contain:
    # 1. DIAMOND profiles for mapping to orthogroups: diamond_profile_sequences.fa.db.dmnd
    # 2. OX database: trees with various lookup data required

    args = parser.parse_args()
    run(args.query, args.out, args.db)
    



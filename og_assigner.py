import sys
import gzip
import csv
import argparse
import subprocess
from collections import Counter


class OGAssigner(object):
    """
    Assign a sequence to an orthogroup
    """
    def assign(self, infn):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        raise NotImplemented()


class OGAssignDIAMOND(OGAssigner):
    """
    Assign using DIAMOND
    """
    def __init__(self, d_db):
        self.d_db = d_db

    def assign(self, infn):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        fn_results = self.run_diamond(infn, infn)
        iog, q_tree = self.og_from_diamond_results(fn_results)
        return iog, q_tree
    
    def run_diamond(self, fn_query, fn_out_base):
        """
        Run DIAMOND against database of orthogroup profiles
        Args:
            fn_query - input query FASTA filename
            fn_out_base - base filename on which to create internal filenames
        Returns:
            fn_og_results_out - DIAMOND results filename
        """
        fn_db = self.d_db + "diamond_profile_sequences.fa.db.dmnd"
        fn_og_results_out = fn_out_base + ".ogs.txt"
        subprocess.call(["diamond", "blastp", "-d", fn_db, "-q", fn_query, "-o", fn_og_results_out, "--quiet", "-e", "0.001", "--compress", "1"])
        return fn_og_results_out + ".gz"

    def og_from_diamond_results(self, fn_og_results_out):
        """
        Get the OG based on the DIAMOND results
        Args:
            fn_og_results_out - fn of compressed DIAMOND results
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
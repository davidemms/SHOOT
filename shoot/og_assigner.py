import os 
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
    def assign(self):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            og - index of OG & subpart "int.int" or "int"
        """
        raise NotImplemented()


class OGAssignDIAMOND(OGAssigner):
    """
    Assign using DIAMOND
    """
    def __init__(self, d_db, nthreads, profiles_db_name="diamond_profile_sequences.fa.db.dmnd"):
        self.d_db = d_db
        self.nthreads = nthreads
        self.profiles_db_name = profiles_db_name

    def assign(self, infn, q_ultra_sens=False):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        ext_cleaned = ".sh.cleaned"
        if infn.endswith(ext_cleaned):
            fn_base = infn.replace(ext_cleaned, "")
        else:
            fn_base = infn
        fn_results = self.run_diamond(infn, fn_base, q_ultra_sens=q_ultra_sens)
        iog = self.og_from_diamond_results(fn_results)
        if not q_ultra_sens and iog is None:
            # Try again with a more sensitive search
            fn_results = self.run_diamond(infn, fn_base, q_ultra_sens=True)
            iog = self.og_from_diamond_results(fn_results)
        return iog
    
    def run_diamond(self, fn_query, fn_out_base, q_ultra_sens=False):
        """
        Run DIAMOND against database of orthogroup profiles
        Args:
            fn_query - input query FASTA filename
            fn_out_base - base filename on which to create internal filenames
        Returns:
            fn_og_results_out - DIAMOND results filename
        """
        fn_db = self.d_db + self.profiles_db_name
        fn_og_results_out = fn_out_base + ".sh.ogs.txt"
        with open(os.devnull, 'w') as FNULL:
            cmd_list = ["diamond", "blastp", "-d", fn_db, "-q", fn_query, "-o", fn_og_results_out, "--quiet", "-e", "0.001", "--compress", "1", "-p", str(self.nthreads)]
            if q_ultra_sens:
                cmd_list += ["--ultra-sensitive"]
            subprocess.call(cmd_list,stdout=FNULL, stderr=FNULL)
        return fn_og_results_out + ".gz"

    def og_from_diamond_results(self, fn_og_results_out, q_ignore_sub = False):
        """
        Get the OG based on the DIAMOND results
        Args:
            fn_og_results_out - fn of compressed DIAMOND results
            q_ignore_sub - ignore any subtrees and just look at overall OGs
        Returns:
            iog        : str "int.int" or "int" otherwise None
        """
        ogs = []
        scores = []
        with gzip.open(fn_og_results_out, 'rt') as infile:
            reader = csv.reader(infile, delimiter="\t")
            for l in reader:
                og = l[1].split("_")[0]
                ogs.append(og.split(".",1)[0] if q_ignore_sub else og)
                scores.append(float(l[-2]))
        sortedTuples = sorted(zip(scores, ogs))
        scores = [i for i, j in sortedTuples]
        ogs = [j for i, j in sortedTuples]
        if len(ogs) == 0:
            # no match
            return None
        else:
            return ogs[0]
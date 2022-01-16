import os 
import sys
import gzip
import csv
import argparse
import subprocess
from collections import Counter

import og_assigner


class OGAssigner_mmseqs(og_assigner.OGAssigner):
    """
    Assign using MMseqs
    """
    def __init__(self, d_db, profiles_db_name="all_msa.subtrees_2000.profile"):
        self.d_db = d_db
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
        fn_results = self.run_mmseqs(infn, fn_base, q_ultra_sens=q_ultra_sens)
        ogs, scores = self.og_from_mmseqs_results(fn_results)
        if not q_ultra_sens and len(ogs) == 0:
            # Try again with a more sensitive search
            fn_results = self.run_mmseqs(infn, fn_base, q_ultra_sens=True)
            ogs, scores = self.og_from_mmseqs_results(fn_results)
        return ogs, scores
    
    def run_mmseqs(self, fn_query, fn_out_base, q_ultra_sens=False):
        """
        Run MMseqs against database of orthogroup profiles
        Args:
            fn_query - input query FASTA filename
            fn_out_base - base filename on which to create internal filenames
        Returns:
            fn_og_results_out - DIAMOND results filename
        """
        fn_db = self.d_db + self.profiles_db_name
        fn_og_results_out = fn_out_base + ".sh.ogs.m8"
        cmd_list = ["mmseqs", "easy-search", fn_query, fn_db, fn_og_results_out, "/tmp", "-v", "0", "--db-load-mode", "2"]
        # print(" ".join(cmd_list))
        if q_ultra_sens:
            cmd_list += ["-s", "7.5"]
        else:
            cmd_list += ["-s", "2"]
        subprocess.call(cmd_list,)
        return fn_og_results_out

    def og_from_mmseqs_results(self, fn_og_results_out, q_ignore_sub = False):
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
        bits_minus = []
        with open(fn_og_results_out, 'rt') as infile:
            reader = csv.reader(infile, delimiter="\t")
            for l in reader:
                og = l[1]
                ogs.append(og.split(".",1)[0] if q_ignore_sub else og)
                scores.append(float(l[-2]))
                bits_minus.append(-float(l[-1]))
        if len(ogs) == 0:
            return None
        sortedTuples = sorted(zip(scores, ogs))
        scores = [i for i, j in sortedTuples]
        ogs = [j for i, j in sortedTuples]
        if len(ogs) > 0:
            # filter out all ogs other than those which are potentially worth considering
            sliver = np.nextafter(0, 1)
            scores_ml10 = [-np.log10(s+sliver) for s in scores]
                s0 = scores_ml10[0]
            scores_ml10 = [s for s in scores_ml10 if s0-s<10]
            ogs = ogs[:len(scores_ml10)]
            scores = scores[:len(scores_ml10)]
        else:
            return ogs, scores
        return ogs[0]

import os 
import sys
import gzip
import csv
import argparse
import subprocess
from collections import Counter

import numpy as np


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
    def __init__(self, d_db, nthreads, q_all_seqs=False):
        self.d_db = d_db
        self.nthreads = nthreads
        self.q_all_seqs_default = q_all_seqs
        self.profiles_db_name = "diamond_profile_sequences.fa.db.dmnd"
        self.all_seqs_db_name = "diamond_all_sequences.fa.db.dmnd"

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
        ogs, scores = self.og_from_diamond_results(fn_results)
        if ((not q_ultra_sens) or (not self.q_all_seqs_default)) and len(ogs) == 0:
            # Try again with a more sensitive search & all sequences
            fn_results = self.run_diamond(infn, fn_base, q_ultra_sens=True, q_all_seqs=True)
            ogs, scores = self.og_from_diamond_results(fn_results)
        return ogs, scores
    
    def run_diamond(self, fn_query, fn_out_base, q_ultra_sens=False, q_all_seqs=False):
        """
        Run DIAMOND against database of orthogroup profiles
        Args:
            fn_query - input query FASTA filename
            fn_out_base - base filename on which to create internal filenames
        Returns:
            fn_og_results_out - DIAMOND results filename
        """
        fn_db = self.d_db + (self.all_seqs_db_name if q_all_seqs else self.profiles_db_name) 
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
            iog        : List[str] "int.int" or "int" of ordered, ambiguous assignments
        Info:
            Hits to other groups are returned if np.log10 difference is less than 10 
            and -np.log10 score > difference.
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
        if len(ogs) > 0:
            # filter out all ogs other than those which are potentially worth considering
            sliver = np.nextafter(0, 1)
            scores_ml10 = [-np.log10(s+sliver) for s in scores]
                s0 = scores_ml10[0]
            # in a test of 15k sequences only 12 passed the first test but failed s>(s0-s)
            # it is not worth arguing over whether it's a good second criteris 
            # scores = [s for s in scores if s0-s<10 and s>(s0-s)]
            scores_ml10 = [s for s in scores_ml10 if s0-s<10]
            ogs = ogs[:len(scores_ml10)]
            scores = scores[:len(scores_ml10)]
        else:
            return ogs, scores

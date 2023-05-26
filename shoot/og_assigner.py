import os
import sys
import re
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
        Returns:
            assignations: a dictionay of gene name with their:
                iogs - List[index of OG]
                scores - List[e-values]
        """
        ext_cleaned = ".sh.cleaned"
        if infn.endswith(ext_cleaned):
            fn_base = infn.replace(ext_cleaned, "")
        else:
            fn_base = infn
        fn_results = self.run_diamond(infn, fn_base, q_ultra_sens=q_ultra_sens)
        assignations = self.og_from_diamond_results(fn_results)

        # Get genes.
        infile = open(infn, 'r')
        fasta_data = infile.read()
        genes = dict(re.findall(r"^>([^\n]+)\n([^>]+)", fasta_data, re.MULTILINE))
        infile.close()

        reassign = []
        for gene in genes:
            if gene not in assignations:
                assignations[gene] = {
                    'ogs': [],
                    'scores': [],
                }
            if ((not q_ultra_sens) or (not self.q_all_seqs_default)) and (len(assignations[gene]['ogs']) == 0):
                reassign.append(gene)
        if 0 < len(reassign):
            # Generate a new FASTA with the gene to reassign
            new_fasta_data = ''
            for gene in reassign:
                new_fasta_data = new_fasta_data + '>' + gene + genes[gene] + '\n'

            sensfile_name = infn + '.sens.fa'
            with open(sensfile_name, 'w') as sensfile:
                sensfile.write(new_fasta_data)

            # Try again with a more sensitive search & all sequences
            fn_results = self.run_diamond(sensfile_name, fn_base, q_ultra_sens=True, q_all_seqs=True)
            assignations = self.og_from_diamond_results(fn_results)
        return assignations
    
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
        with open(fn_og_results_out + '.log', 'w') as log:
            cmd_list = ["diamond", "blastp", "-d", fn_db, "-q", fn_query, "-o", fn_og_results_out, "--quiet", "-e", "0.001", "--compress", "1", "-p", str(self.nthreads)]
            if q_ultra_sens:
                cmd_list += ["--ultra-sensitive"]
            subprocess.call(cmd_list,stdout=log, stderr=log)
        return fn_og_results_out + ".gz"

    @staticmethod
    def og_from_diamond_results(fn_og_results_out, q_ignore_sub=False):
        """
        Get the OG based on the DIAMOND results
        Args:
            fn_og_results_out - fn of compressed DIAMOND results
            q_ignore_sub - ignore any subtrees and just look at overall OGs
        Returns:
            iog: a dictionay which keys are gene names and values are dicts:
              ogs: the ordered list of assigned OGs
              scores: the list of corresponding e-values scores
        Info:
            Hits to other groups are returned if np.log10 difference is less than 10 
            and -np.log10 score > difference.
        """
        assignations = {}
        ogs = []
        scores = []
        with gzip.open(fn_og_results_out, 'rt') as infile:
            reader = csv.reader(infile, delimiter="\t")
            for l in reader:
                gene = l[0]
                if not gene in assignations:
                    assignations[gene] = {
                        'ogs': [],
                        'scores': [],
                    }
                og = l[1].split("_")[0]
                assignations[gene]['ogs'].append(og.split(".",1)[0] if q_ignore_sub else og)
                assignations[gene]['scores'].append(float(l[-2]))
        for gene, values in assignations.items():
            sortedTuples = sorted(zip(values['scores'], values['ogs']))
            scores = [i for i, j in sortedTuples]
            ogs = [j for i, j in sortedTuples]
            unique = [i for i in range(len(ogs)) if ogs[i] not in ogs[:i]]
            ogs = [ogs[i] for i in unique]
            scores = [scores[i] for i in unique]
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
            assignations[gene]['ogs'] = ogs
            assignations[gene]['scores'] = scores
        return assignations

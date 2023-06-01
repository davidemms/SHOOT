"""
Use RAxML's EPA to placed the gene in the tree. This is faster than IQTREE and therefore
supports placement of a gene in a larger tree (perhaps 2500 taxon tree in 10s rather
than a 500 taxon tree with IQTREE). USeing larger trees allows higher accuracy as
we don't have to place a gene so specifically (in such a small subtree) using DIAMOND
and instead can rely more on tree inference methods to put a gene in the correct 
place.
"""
import json
import os
import re
import subprocess

import ete3

import msa_grafter

class MSAGrafter_EPA(msa_grafter.MSAGrafter):
    def __init__(self, *args, **kwargs):
        super(MSAGrafter_EPA, self).__init__(*args, **kwargs)
        self.re_disallowed_seq_chars = re.compile("[^-ABCDEFGHIKLMNPQRSTVWYZ\n]")

    def run_tree_inference(self, og_part, n_seqs_123many, fn_msa, q_subtree, query_names):
        """
        Run EPA
        Args:
            og_part - OG ID string, either iog or iog.ipart
            n_seqs_123many - Lower bound on number of seqs in original tree, one of {1,2,3,4}
            fn_msa - the MSA
            q_subtree - is the gene assigned to a subtree
            query_names - list of query genes
        Implementation:
        - >=4 genes: unroot the starting tree & run iqtree
        - 3 genes: run iqtree from scratch
        """
        if n_seqs_123many < 3:
            return
        if n_seqs_123many >= 4:
            # have an original tree with 4 or more taxa + new seq
            # then will have a complete tree
            fn_tree_orig = self.db.fn_tree(og_part)
            # fn_unrooted = fn_tree_orig + ".unroot.tre"
            # if not os.path.exists(fn_unrooted):
            #     try:
            #         t = ete3.Tree(fn_tree_orig)
            #     except ete3.parser.newick.NewickError:
            #         t = ete3.Tree(fn_tree_orig, format=1)
            #     t.unroot()
            #     fn_unrooted = fn_tree_orig + ".un.tre"
            #     t.write(outfile=fn_unrooted)
            # Split the MSA into reference and query files
            fn_ref, fn_query = self.split_reference_query(fn_msa, query_names)
            # Run EPA
            out_dir = fn_msa + "_epa/"
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            fn_tree_bifur = fn_tree_orig + ".bifur.tre"
            fn_tree_input = fn_tree_bifur if os.path.exists(fn_tree_bifur) else fn_tree_orig
            # print("epa-ng -T 32 -m LG --redo --tree %s --ref-msa %s --query %s --preserve-rooting on --outdir %s" % (fn_tree_input, fn_ref, fn_query, out_dir))
            with open(os.devnull, 'w') as FNULL:
                subprocess.call("epa-ng -T %d -m LG --redo --tree %s --ref-msa %s --query %s --preserve-rooting on --outdir %s" % (self.nthreads, fn_tree_input, fn_ref, fn_query, out_dir), shell=True, stdout=FNULL, stderr=FNULL)
            results_fn = out_dir + "epa_result.jplace"
            fn_tree_new = self.place_gene_gappa(results_fn)
            self.add_support(results_fn, fn_tree_new)
        elif n_seqs_123many == 3:
            # have minimum number of sequences for a tree, by no original tree
            subprocess.call("iqtree -nt %d -quiet -m LG -fast -redo -s %s" % (self.nthreads, fn_msa), shell=True)
            fn_tree_new = fn_msa + ".treefile"
        return fn_tree_new


    # def extract_epa_tree(self, results_fn):
    #     """
    #     Extract the Newick tree string from the EPA results and write a newick tree file
    #     Args:
    #         results_fn - EPA results file, typically epa_result.jplace
    #     Returns:
    #         fn_tree - Newick tree file
    #     """
    #     fn_tree = results_fn + ".tre"
    #     padding = '  "tree": "'
    #     with open(results_fn, 'r') as infile:
    #         next(infile)
    #         t_str = next(infile)[len(padding):]
    #         # tree string is surrounded by quotes
    #         while not t_str.endswith('"'):
    #             t_str = t_str[:-1]
    #         t_str = t_str[:-1]
    #     with open(fn_tree, 'w') as outfile:
    #         outfile.write(re.sub('{\d+}', '', t_str))
    #     # ete can't be used as it can't read EPA's stupid format
    #     return fn_tree

    @staticmethod
    def add_support(results_fn, fn_tree_new):
        t = ete3.Tree(fn_tree_new)
        max_support = max(n.support for n in t.traverse())
        gene, support = MSAGrafter_EPA.get_gene_support(results_fn)
        if max_support > 1.0:
            support *= 100.
        n = (t & gene).up
        n.support = support
        t.write(outfile=fn_tree_new)

    @staticmethod
    def get_gene_support(results_fn):
        i_support = 2
        with open(results_fn, 'r') as infile:
            x = json.load(infile)
        gene = x['placements'][0]['n'][0]
        support = x['placements'][0]['p'][0][i_support]
        return gene, support

    def place_gene_gappa(self, jplace_fn):
        d, fn = os.path.split(jplace_fn)
        with open(os.devnull, 'w') as FNULL:
            subprocess.call(["gappa", "examine", "graft", "--jplace-path", jplace_fn, "--out-dir", d, "--file-prefix", fn + ".", "--allow-file-overwriting", "--threads", "1"],stdout=FNULL, stderr=FNULL)
        return jplace_fn + ".epa_result.newick"

    def split_reference_query(self, fn_msa, query_names):
        """
        Split an MSA into separate files for the reference & query sequences, as 
        required by EPA
        Args:
            fn_msa - input MSA (last sequence is the query)
        Returns:
            fn_ref - the reference MSA filename
            fn_query - the query MSA filename
        """
        fn_ref = fn_msa + ".ref.fa"
        fn_query = fn_msa + ".query.fa"
        infile = open(fn_msa, 'r')
        fasta_data = infile.read()
        genes = dict(re.findall(r"^>([^\n]+)\n([^>]+)", fasta_data, re.MULTILINE))
        infile.close()
        with open(fn_ref, 'w') as out_ref, open(fn_query, 'w') as out_query:
            for gene, seq in genes.items():
                if gene in query_names:
                    out_query.write('>' + gene + '\n' + seq + '\n')
                else:
                    out_ref.write('>' + gene + '\n' + seq + '\n')

        return fn_ref, fn_query

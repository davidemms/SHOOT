"""
Add the new sequence to an existing MSA and then infer the tree from that. If possible
use the original tree as a start point - but not sure if I've got anything that can do this.
"""
import os
import ete3
import subprocess

import database

class MSAGrafter(object):
    def __init__(self, d_db):
        """
        Add a gene to a tree by extending the MSA and re-inferring the tree.
        Args:
            d_db - Database directory
        """
        self.db = database.Database(d_db)

    def add_gene(self, iog, infn, method = "iqtree"):
        """
        Args:
            iog - the OG to search in
            infn - FASTA filename containing the gene sequence
        """
        warn_string = ""
        fn_msa_orig = self.db.fn_msa(iog)
        if not os.path.exists(fn_msa_orig):
            fn_msa_orig = self.db.fn_seqs(iog)
        n_seqs_orig_lower_bound = self.n_seqs_1234(fn_msa_orig)
        fn_msa_new = infn + ".msa.fa"

        subprocess.call("mafft --retree 1 --maxiterate 0 --nofft --thread 16 --quiet --add %s %s > %s" % (infn, fn_msa_orig, fn_msa_new), shell=True)
        # print(fn_msa_new)
        if n_seqs_orig_lower_bound > 3:
            # have an original tree with 4 or more taxa + new seq
            # then won't have a complete tree
            fn_tree_orig = self.db.fn_tree(iog)
            try:
                t = ete3.Tree(fn_tree_orig)
            except ete3.parser.newick.NewickError:
                t = ete3.Tree(fn_tree_orig, format=1)
            t.unroot()
            fn_unrooted = fn_tree_orig + ".un.tre"
            t.write(outfile=fn_unrooted)
            subprocess.call("iqtree -nt 32 -quiet -m LG -fast -redo -g %s -s %s" % (fn_unrooted, fn_msa_new), shell=True)
            fn_tree_new = fn_msa_new + ".treefile"
        elif n_seqs_orig_lower_bound == 3:
            # have minimum number of sequences for a tree, by no original tree
            subprocess.call("iqtree -nt 32 -quiet -m LG -fast -redo -s %s" % fn_msa_new, shell=True)
            fn_tree_new = fn_msa_new + ".treefile"
        else:
            # not enough taxa for a tree
            genes = []
            with open(fn_msa_new, 'r') as infile:
                for l in infile:
                    if l.startswith(">"):
                        genes.append(l[1:].rstrip())
            newick_str = "(" + ",".join(genes) + ");"
            fn_final_tree = fn_msa_new + ".tre"
            with open(fn_final_tree, 'w') as outfile:
                outfile.write(newick_str)
            warn_string = "Tree could not be rooted"
            return fn_final_tree, warn_string

        if n_seqs_orig_lower_bound < 4:
            # no original tree to root with
            warn_string = "Tree could not be rooted"
            return fn_tree_new, warn_string

        # Root the tree as it was previously rooted
        try:
            t_orig = ete3.Tree(fn_tree_orig)
        except ete3.parser.newick.NewickError:
            t_orig = ete3.Tree(fn_tree_orig, format=1)
        chs = t_orig.get_children()
        n = [len(ch) for ch in chs]
        i = n.index(min(n))
        outgroup_names = chs[i].get_leaf_names()
        if method == "iqtree":
            outgroup_names = self.iqtree_names_adjust(outgroup_names)
        #print(outgroup_names)
        t_new = ete3.Tree(fn_tree_new, format=1)
        for n in t_new.traverse():
            if not n.is_leaf() and n.name != "":
                try:
                    n.support = n.name.split("/")[1]
                except:
                    pass
        if not self.check_monophyly(t_new, outgroup_names):
            # find the largest set that is monophyletic
            # inexact method, add genes until fails
            n_genes = 0
            while n_genes < len(outgroup_names):
                n_genes += 1
                test_mono_outgroup_names = outgroup_names[:1]
                if not self.check_monophyly(t_new, test_mono_outgroup_names):
                    outgroup_names = test_mono_outgroup_names[:-1]
        # root here
        if len(outgroup_names) == 1:
            t_new.set_outgroup(outgroup_names[0])
        else:
            root = t_new.get_common_ancestor(outgroup_names)
            if root != t_new:
                print("Rerooting")
                t_new.set_outgroup(root)
        # transfer_support_values(t_orig, t_new)
        fn_final_tree = infn + ".grafted.msa.tre"
        t_new.write(outfile=fn_final_tree)
        return fn_final_tree, warn_string

    @staticmethod
    def iqtree_names_adjust(genes):
        return [g.replace("|", "_") for g in genes]

    @staticmethod
    def n_seqs_1234(infn):
        """
        Returns the number of sequences if it is less than 4, otherwise returns 4
        """
        n = 0
        with open(infn, 'r') as infile:
            for l in infile:
                if l.startswith(">"):
                    n+=1
                    if n >= 4:
                        break
        return n

    @staticmethod
    def check_monophyly(node, taxa):
        """ This should be used as a wrapper of the ete method since that method returns false for a single node. It should return true
        Args:
            node - the node under which to search
            taxa - the list of taxa
        Returns
            bool - are the taxa monophyletic
        """
        if len(taxa) == 1:
            return (list(taxa)[0] in node)
        else: 
            return node.check_monophyly(target_attr='name', values=taxa)[0]



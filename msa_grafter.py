"""
Add the new sequence to an existing MSA and then infer the tree from that. If possible
use the original tree as a start point - but not sure if I've got anything that can do this.
"""
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

    def add_gene(self, iog, infn, method="iqtree"):
        """
        Args:
            iog - the OG to search in
            infn - FASTA filename containing the gene sequence
        """
        fn_msa_orig = self.db.fn_msa(iog)
        fn_msa_new = infn + ".msa.fa"

        # print("mafft --quiet --add %s %s > %s" % (infn, fn_msa_orig, fn_msa_new))
        subprocess.call("mafft --quiet --add %s %s > %s" % (infn, fn_msa_orig, fn_msa_new), shell=True)
        print(fn_msa_new)
        if method == "iqtree":
            fn_tree_orig = self.db.fn_tree(iog)
            t = ete3.Tree(fn_tree_orig)
            t.unroot()
            fn_unrooted = fn_tree_orig + ".un.tre"
            t.write(outfile=fn_unrooted)
            subprocess.call("iqtree -quiet -g %s -s %s" % (fn_unrooted, fn_msa_new), shell=True)
            fn_tree_new = fn_msa_new + ".treefile"
        elif method == "fasttree":
            fn_tree_new = infn + ".msa.tre"
            # print("FastTree -quiet %s > %s" % (fn_msa_new, fn_tree_new))
            subprocess.call("FastTree -quiet %s > %s" % (fn_msa_new, fn_tree_new), shell=True)
            print(fn_tree_new)
        else:
            raise NotImplementedError("Unknown method: %s" % method)
        # Root the tree as it was previously rooted
        fn_tree = self.db.fn_tree(iog)
        t_orig = ete3.Tree(fn_tree)
        chs = t_orig.get_children()
        n = [len(ch) for ch in chs]
        i = n.index(min(n))
        outgroup_names = chs[i].get_leaf_names()
        print(outgroup_names)
        t_new = ete3.Tree(fn_tree_new)
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
        fn_final_tree = infn + ".grafted.msa.tre"
        t_new.write(outfile=fn_final_tree)
        return fn_final_tree


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



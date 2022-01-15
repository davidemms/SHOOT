# -*- coding: utf-8 -*-
import os
import sys
import glob 

import ete3

from . import og_assigner
from . import database
from . import deeptree

class OGAssignQuartets(og_assigner.OGAssignDIAMOND):
    """
    Use diamond to assign the gene to the overall homolog tree and then use the 
    quartets method to assign it to the correct subtree
    """
    def __init__(self, d_db):
        self.d_db = d_db
        self.db = database.Database(d_db)

    def assign(self, infn):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        fn_results = self.run_diamond(infn, infn)
        iog = self.og_from_diamond_results(fn_results, q_ignore_sub=True)
        if iog is None:
            # Try again with a more sensitive search
            fn_results = self.run_diamond(infn, infn, q_ultra_sens=True)
            iog = self.og_from_diamond_results(fn_results)
        if iog is None:
            return iog
        print("Got: " + iog)
        # if there is no supertree then there is nothing more to do
        if not os.path.exists(self.db.fn_tree_super(int(iog.split(".")[0]))):
            return iog

        # otherwise, use the quartets method to place it in the right subtree
        with open(infn, 'r') as infile:
            query_gene = next(infile)[1:].rstrip()
        fn_msa_db = self.db.fn_msa(iog)
        genes_in_og = []
        with open(fn_msa_db, 'r') as infile:
            for l in infile:
                if not l.startswith(">"):
                    continue
                genes_in_og.append(l[1:].rstrip())
        if not query_gene in genes_in_og:
            # wrong overall tree, so it's going to be wrong whatever
            return iog
        # Otherwise, open the gene tree, remove our gene and use shoot to put it 
        # back in. Root the tree and see where it ends up
        t = ete3.Tree(self.db.fn_tree(iog))
        # get outgroup for later rooting
        N = len(t)
        assert(len(t.children) == 2)
        for ch in t.children:
            if 2*len(ch) <= N:
                break
        outgroup = ch.get_leaf_names()
        n = t & query_gene
        n.detach()
        t.unroot()
        fn_tree_removed = infn + ".removed.tre"
        t.write(outfile=fn_tree_removed)
        msa_fn = self.db.fn_msa(iog)
        fn_tree_out = infn + ".regrafted.tre"
        deeptree.main(msa_fn, fn_tree_out, guide=fn_tree_removed)

        t = ete3.Tree(fn_tree_out)
        try:
            t.set_outgroup(t.get_common_ancestor(outgroup))
        except:
            pass
        n = t & query_gene
        sister_clade = n.up.get_leaf_names()
        sister_clade.remove(query_gene)

        # see what subtree we've put it in
        d_gene_to_subtree = dict()
        for fn in glob.glob(self.d_db + "Gene_Trees/subtrees/msa_sub/OG%s.*.fa" % iog):
            t = os.path.basename(fn)[2:].split(".")
            og_part = t[0] + "." + t[1]
            with open(fn, 'r ') as infile:
                for l in infile:
                    if not l.startswith(">"):
                        continue
                    g = l[1:].rstrip()
                    d_gene_to_subtree[g] = og_part
        # print(d_gene_to_subtree.keys()[:5])
        og_part = d_gene_to_subtree[next(g for g in sister_clade)]
        # print("Got: " + og_part)
        return og_part

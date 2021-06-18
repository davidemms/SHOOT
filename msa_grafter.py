"""
Add the new sequence to an existing MSA and then infer the tree from that. If possible
use the original tree as a start point - but not sure if I've got anything that can do this.
"""
import os
import subprocess

import ete3

import database

class MSAGrafter(object):
    def __init__(self, d_db):
        """
        Add a gene to a tree by extending the MSA and re-inferring the tree.
        Args:
            d_db - Database directory
        """
        self.db = database.Database(d_db)

    def add_gene(self, og_part, infn, name_orig = None, name_temp = None, method = "iqtree"):
        """
        Args:
            og_part - the OG or OG.PART to search in
            infn - FASTA filename containing the gene sequence
            name_orig - Gene name from user
            name_temp - Working gene name to prevent clash
        Notes:
            - If name_orig and name_temp are both not None, changes query name back
            to name_orig in the final tree
        """
        q_subtree = "." in og_part
        warn_string = ""
        fn_msa_orig = self.db.fn_msa(og_part)
        # if not os.path.exists(fn_msa_orig):
        #     fn_msa_orig = self.db.fn_seqs(iog)
        n_seqs_orig_lower_bound = self.n_seqs_1234(fn_msa_orig)
        fn_msa_new = infn + ".msa.fa"

        subprocess.call("mafft --retree 1 --maxiterate 0 --nofft --thread 16 --quiet --add %s %s > %s" % (infn, fn_msa_orig, fn_msa_new), shell=True)
        # print(fn_msa_new)
        if n_seqs_orig_lower_bound >= 4:
            # have an original tree with 4 or more taxa + new seq
            # then will have a complete tree
            fn_tree_orig = self.db.fn_tree(og_part)
            fn_unrooted = fn_tree_orig + ".unroot.tre"
            if not os.path.exists(fn_unrooted):
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
            if not q_subtree:
                # write the tree, otherwise we will graft if in to the full tree below
                newick_str = "(" + ",".join(genes) + ");"
                fn_final_tree = fn_msa_new + ".tre"
                with open(fn_final_tree, 'w') as outfile:
                    outfile.write(newick_str)
                warn_string = "Tree could not be rooted"
                return fn_final_tree, warn_string

        # Rooting
        if not q_subtree:
            if n_seqs_orig_lower_bound < 4:
                # no original tree to root with
                warn_string = "Tree could not be rooted"
                if (name_orig is not None and name_temp is not None):
                    t = ete3.Tree(fn_tree_new)
                    n = t & name_temp
                    n.name = name_orig
                    t.write(outfile = fn_tree_new)
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
            if (name_orig is not None and name_temp is not None):
                n = t_new & name_temp
                n.name = name_orig
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
        else:
            # SUBTREE 
            n_taxa_lower_bound = n_seqs_orig_lower_bound + 1
            # graft the subtree into the already rooted supertree
            # new subtree size. 2 can't occur (query¸ single subtree gene, outgroup gene)
            # Options
            # 3  split branch from supertree in half, have the node this dist and each of the genes this dist from node
            # >=4: infer tree, root on outgroup, get newick subtree
            iog, i_part_hit = list(map(int, og_part.split(".")))
            try:
                t_sup = ete3.Tree(self.db.fn_tree_super(iog))
            except ete3.parser.newick.NewickError:
                t_sup = ete3.Tree(self.db.fn_tree_super(iog), format=1)
            parent_node_name = "PART." + og_part.split(".")[1]
            n_parts = len(t_sup)
            # now graft in the inferred gene tree
            if n_taxa_lower_bound == 3:
                p = t_sup & parent_node_name
                d = p.dist / 2.0
                genes = [g for g in genes if not g.startswith("SHOOTOUTGROUP")]
                if (name_orig is not None and name_temp is not None):
                    genes[genes.index(name_temp)] = name_orig
                nwk_sub = "(%s:%f,%s:%f)1:" % (genes[0], d, genes[1], d)
            elif n_taxa_lower_bound >= 4:
                t_new = ete3.Tree(fn_tree_new, format=1)
                if ((name_orig is not None) and (name_temp is not None)):
                    n = t_new & name_temp
                    n.name = name_orig
                g_out = next(g for g in t_new.get_leaf_names() if g.startswith("SHOOTOUTGROUP"))
                t_new.set_outgroup(g_out)
                n = next(ch for ch in t_new.children if ch.name != g_out)
                nwk_sub = n.write()
                # need to remove the final distance, should end ")"
                while nwk_sub[-1] != ")":
                    nwk_sub = nwk_sub[:-1]
                nwk_sub += "1:"
            else:
                print("ERROR, unreachable state")
            # now put all the subtrees into the super tree and return
            d_sub_nwks = []
            for i_part in range(n_parts):
                if i_part == i_part_hit:
                    d_sub_nwks.append(nwk_sub)
                else:
                    with open(self.db.fn_trees_sub_without(iog, i_part), 'r') as infile:
                        this_nwk = next(infile).rstrip()
                        assert(this_nwk.endswith(";"))
                        this_nwk = this_nwk[:-1]
                        # each will have a distance on its root that is duplicated 
                        # in the super tree
                        # Either: "Monodelphis_domestica_A0A5F8H6S2:3.96746;""
                        # Or: ...295)1:0.04477;
                        # It needs to end ":"
                        if this_nwk.startswith("("):
                            while this_nwk[-1] != ")":
                                this_nwk = this_nwk[:-1]
                            this_nwk += "1:"
                        else:
                            while this_nwk[-1] != ":":
                                this_nwk = this_nwk[:-1]
                    d_sub_nwks.append(this_nwk)
            nwk_super = t_sup.write()
            nwk_sup_splits = nwk_super.split("PART.")
            nwk = nwk_sup_splits[0]
            for c in nwk_sup_splits[1:-1]:
                i, remainder = c.split(":", 1)
                nwk += d_sub_nwks[int(i)] + remainder
            nwk += nwk_sup_splits[-1]
            print(nwk.count("("))
            print(nwk.count(")"))
            fn_final_tree = infn + ".grafted.msa.tre"
            with open(fn_final_tree, 'w') as outfile:
                outfile.write(nwk + "\n")
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



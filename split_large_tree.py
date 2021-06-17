import os 
import argparse

import ete3

import fasta_writer

def split_tree(fn_tree, fn_msa, n_taxa, q_outgroup):
    """
    Split a tree into smaller subtrees
    Args:
        fn_tree - input tree filename
        fn_msa - input MSA filename
        n_taxa - number of taxa to aim for in each sub tree
        q_outgroup - include an outgroup taxon in each subtree to allow downstream root placement
    Returns:
        n - number of profile sequences required
    Implementation:
        Each tree with more than 5 sequences requires 5 profile sequences. Below 
        that, it is the number of sequences in the tree.
        Clades left in the tree are treated as requiring 5 profile sequences too.
    """
    n_profile = 0
    fw = fasta_writer.FastaWriter(fn_msa)
    d, fn = os.path.split(fn_tree)
    if d == "":
        d = "./"
    dout_main = d + "/subtrees/"
    d_out_sup = dout_main + "super/"
    d_out_sub = dout_main + "sub/"
    d_out_msa_sub = dout_main + "msa_sub/"
    for d_to_make in [dout_main, d_out_sup, d_out_sub, d_out_msa_sub]:
        if not os.path.exists(d_to_make):
            os.mkdir(d_to_make)
    t = ete3.Tree(fn_tree)
    fn_out_pat = d_out_sub + fn + ".%d.tre"
    fn_out_msa_pat = d_out_msa_sub + fn + ".%d.fa"
    fn_out_mega_pat = d_out_sup + fn + ".super.tre"
    i_part = 0
    if len(t) <= n_taxa:
        t.write(outfile = fn_out_pat % i_part)
        # print("No need to split tree")
        return min(5, len(t))
    # traverse tree
    sizes = []
    stop_fn = lambda node : len(node) <= n_taxa
    for n in t.traverse("preorder", is_leaf_fn = stop_fn):
        l = len(n)
        if l <= n_taxa:
            # split here 
            n_profile += min(5, l)
            n.write(outfile = fn_out_pat % i_part)
            # get a representative from each split
            fw.WriteSeqsToFasta(n.get_leaf_names(), fn_out_msa_pat % i_part)
            sizes.append(l)
            p = n.up
            n = n.detach()
            if len(n) > 2:
                n.unroot()
            n.write(outfile = (fn_out_pat % i_part) + ".unroot.tre")
            ch = p.add_child(name = "PART.%d-%d_genes" % (i_part, l))
            i_part += 1
    # print(sorted(sizes))
    t.write(outfile=fn_out_mega_pat)
    t.unroot()
    t.write(outfile=fn_out_mega_pat + ".unroot.tre")
    # print(fn_out_mega_pat)
    return n_profile

def count_profiles():
    n_target_taxa = 500
    d = "/lv01/data/emms/SHOOT/DATA/UniProt_RefProteomes_homologs/"
    fn_tree_pat = d + "Gene_Trees/OG%07d_tree.txt"
    fn_msa_pat = d + "MultipleSequenceAlignments/OG%07d.fa"
    n_trees = 17125
    n_profiles = 0
    for i in range(n_trees):
        if i % 100 == 0:
            print(i)
        n_profiles += split_tree(fn_tree_pat % i, fn_msa_pat % i, n_target_taxa)
    print("%d profile sequences required" % n_profiles)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("tree", help="Input tree file")
    parser.add_argument("msa", help="Input MSA file")
    parser.add_argument("-n", "--ntaxa", help="Number of taxa to aim for for each subtree", type=int, default=500)
    parser.add_argument("-o", "--outgroup", action="store_true",
                        help="Include an outgroup gene in each subtree. Required for SHOOT", )
    args = parser.parse_args()
    split_tree(args.tree, args.msa, args.ntaxa, q_outgroup=args.outgroup)

    # count_profiles()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# import sys
import re
import csv
import argparse
from collections import defaultdict

import og_assigner
import tree_grafter
import msa_grafter
import msa_grafter_epa
# import quartets_pairwise_align

import ete3


gene_name_disallowed_chars_re = '[^A-Za-z0-9_\\-.]'


class Options(object):
    def __init__(
            self, 
            nthreads,
            q_profiles_all=False,
            search_high_sens=False,    # Use high sens. for homolog group search
            q_msa=True,                # Use MSA method for tree inference
            q_mafft_accelerated=True,
            tree_method="epa",         # method to use for tree inference
            nU=None,                   # Upper limit for tree, unless nL is not None 
            nL=None,                   # Exceed nU if alternative is a tree < nL
            q_print=False,             # print tree to stdout 
            q_orthologs=False,         # write orthologs to INFN.sh.orthologs.tsv
            ):
        self.nthreads=nthreads
        self.q_profiles_all=q_profiles_all
        self.search_high_sens=search_high_sens
        self.q_msa=q_msa,
        self.q_mafft_accelerated=q_mafft_accelerated
        self.tree_method=tree_method
        self.nU = nU
        self.nL = nL
        self.q_print = q_print
        self.q_orthologs = q_orthologs


def main(d_db, infn, opts):
    """
    Run SHOOT
    Args:
        d_db - Input directory containing SHOOT database
        infn - Query FASTA filename
        opts - Options object

    Returns:
        fn_tree - Filename for tree or None
    """
    d_db += "/"
    # rapid check for FASTA format
    ok = True
    with open(infn, 'r') as infile:
        acc = next(infile)
        ok = bool(re.match("^>.", acc)) 
        ok = ok and bool(re.match("^[a-zA-Z]+$", next(infile))) 
    if not ok:
        print("ERROR: Input file should be FASTA format")
        return
    # now fix up accession if required
    fn_for_use = clean_fasta(infn)

    og_assign = og_assigner.OGAssignDIAMOND(d_db, opts.nthreads, opts.q_profiles_all)
    og_part = og_assign.assign(fn_for_use, q_ultra_sens=opts.search_high_sens)
    if og_part is not None:
        print("Gene assigned to: OG%s" % og_part)
    else:
        print("No homologs found for gene in this database")
        return 

    warn_str = ""
    if opts.q_msa:
        # do a tree using an MSA
        if "iqtree" == opts.tree_method:
            graft = msa_grafter.MSAGrafter(d_db, opts.nthreads)
        elif "epa" == opts.tree_method:
            graft = msa_grafter_epa.MSAGrafter_EPA(d_db, opts.nthreads)
        else:
            print("ERROR: %s method has not been implemented" % opts.tree_method)
            return
        fn_tree, query_gene, warn_str = graft.add_gene(
                og_part, 
                fn_for_use, 
                infn, 
                q_mafft_acc=opts.q_mafft_accelerated,
                )
    else:
        raise Exception("Option not available")
        # quart = quartets_pairwise_align.PairwiseAlignQuartets(d_db, og_part, fn_for_use)
        # search = tree_grafter.TreeGrafter(quart, d_db)
        # This doesn't feel right, iog is fixed in the constructor of quart and hence
        # in the constructor of search, and yet is passed as a variable here.
        # search.place_gene(iog)              
        # search.place_gene() 
    
    print("Tree: %s" % fn_tree)  
    if warn_str != "":
        print("WARNING: " + warn_str) 
        
    q_need_to_print = opts.q_print
    if opts.nU is not None:
        # if only nL is specified alone that has no effect
        t = ete3.Tree(fn_tree)
        if len(t) > opts.nU:
            node = t & query_gene
            while len(node) < opts.nU:
                node_prev = node
                n_taxa_prev = len(node)
                node = node.up
            # now there are more than nU genes in this tree, step down one
            # unless it is fewer than nL
            node = node_prev if (opts.nL is None or n_taxa_prev >= opts.nL) else node
            nwk_str = node.write()
            with open(fn_tree, 'w') as outfile:
                outfile.write(nwk_str)
            if q_need_to_print:
                print(nwk_str)
                q_need_to_print = False

    if q_need_to_print:
        with open(fn_tree, 'r') as infile:
            print(next(infile).rstrip())   # remove any trailing newline characters

    if opts.q_orthologs:
        fn_ologs = fn_for_use + ".sh.orthologs.tsv"
        write_orthologs(fn_tree, fn_ologs, query_gene)
    return fn_tree        


def clean_fasta(infn):
    """
    If the accession contains problematic characters write a new file and return
    its filename, otherwise return the original filename
    Args:
        infn - User input filename
    Returns:
        fn_for_use - file to use going forward for SHOOT
    """
    with open(infn, 'r') as infile:
        acc = next(infile)
        name = acc.rstrip()[1:]
        name_cleaned = re.sub(gene_name_disallowed_chars_re, '_', name)
        if name == name_cleaned:
            return infn
        # otherwise, need to create a new file
        fn_for_use = infn + ".sh.cleaned"
        with open(fn_for_use, 'w') as outfile:
            outfile.write(">%s\n" % name_cleaned)
            for l in infile:
                outfile.write(l)
        return fn_for_use


def gene_to_species(name):
    """
    Convert a gene name to species name
    Args:
        name - gene name
    Returns 
        species_name - species name
    """
    return "_".join(name.split("_")[:2])


def write_orthologs(fn_tree, fn_ologs, gene_name):
    """
    Write the orthologs of the query gene to file
    Args:
        fn_tree - filename for the SHOOT tree
        fn_ologs - filename to write orthologs
        gene_name - query gene name
    """
    t = ete3.Tree(fn_tree)
    with open(fn_ologs, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter="\t")
        writer.writerow(["Species", "Orthologs"])
        n = t & gene_name
        n_prev = n
        species_query_branch = set()
        while n.up is not None:
            n = n.up
            n_other = [node for node in n.children if node != n_prev]
            genes_other = [gene for node in n_other for gene in node.get_leaf_names()]
            species_other = set([gene_to_species(gene) for gene in genes_other])
            overlap = species_query_branch.intersection(species_other)
            o = len(overlap)
            m = min(len(species_other), len(species_query_branch))
            if o == 0 or (m >= 4 and o == 1) or (m >= 9 and o ==2):
                add_orthologs_to_file(writer, genes_other, overlap)
            n_prev = n
            species_query_branch.update(species_other)


def add_orthologs_to_file(writer, genes, overlap):
    """
    Add the genes below the nodes in n_other to the the csv writer as orthologs
    Args:
        writer - CSV writer
        genes - list of gene names
        overlap - set of species in the overlap (these will not be written)
    """
    species_to_genes = defaultdict(list)
    for g in genes:
        species_to_genes[gene_to_species(g)].append(g)
    for sp in sorted(species_to_genes.keys()):
        if sp in overlap:
            continue
        writer.writerow([sp, ", ".join(species_to_genes[sp])])


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by fol_create_dp.py")
    # parser.add_argument("-m", "--msa", action="store_true", help= "Use an MSA tree")
    # Output options
    parser.add_argument("-u", "--upper", type=int, help= "Upper limit for tree, unless -l")
    parser.add_argument("-l", "--lower", type=int, help= "Exceed -u if alternative is < -l")
    parser.add_argument("-p", "--print_tree", action="store_true", help= "Print tree as final line")
    parser.add_argument("-o", "--orthologs", action="store_true", help="Requires database gene names to be 'genus_species_geneID")
    # search options
    parser.add_argument("--high_sens", action="store_true", 
                        help= "High sensitivity for homology group search")
    parser.add_argument("--profiles_all", action="store_true", 
                        help= "Use all genes for the homology group profiles")
    parser.add_argument("--mafft_defaults", action="store_true", 
                        help= "Use mafft defaults rather than accelerated options")
    parser.add_argument("-t", "--tree_method", default="epa", choices={"epa", "iqtree"})
    # threads
    parser.add_argument("-n", "--nthreads", type=int, default=16)
    args, unknown = parser.parse_known_args()
    opts = Options(
            nthreads=args.nthreads,
            q_profiles_all=args.profiles_all,
            search_high_sens=args.high_sens,
            q_msa=True,
            q_mafft_accelerated=not args.mafft_defaults,
            tree_method=args.tree_method,
            nU=args.upper,
            nL=args.lower, 
            q_print=args.print_tree, 
            q_orthologs=args.orthologs,
            )
    main(args.db, args.infile, opts)


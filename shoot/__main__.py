#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import csv
import argparse
from collections import defaultdict

import og_assigner
import tree_grafter
import msa_grafter
import msa_grafter_epa
import utils
# import quartets_pairwise_align

import ete3


gene_name_disallowed_chars_re = '[^A-Za-z0-9_\\-.]'


class Options(object):
    def __init__(
            self, 
            nthreads,
            q_profiles_all=False,
            search_high_sens=False,    # Use high sens. for homolog group search
            q_mafft_accelerated=True,
            tree_method="epa",         # method to use for tree inference
            nU=None,                   # Upper limit for tree, unless nL is not None 
            nL=None,                   # Exceed nU if alternative is a tree < nL
            q_print=False,             # print tree to stdout 
            q_orthologs=False,         # write orthologs to INFN.sh.orthologs.tsv
            q_ambiguous=False,         # If ambiguity, do each potential tree
            q_full_tree=False,         # Use phylogenetic methods on full tree regardless
            ):
        self.nthreads=nthreads
        self.q_profiles_all=q_profiles_all
        self.search_high_sens=search_high_sens
        self.q_mafft_accelerated=q_mafft_accelerated
        self.tree_method=tree_method
        self.nU = nU
        self.nL = nL
        self.q_print = q_print
        self.q_orthologs = q_orthologs
        self.q_ambiguous = q_ambiguous
        self.q_full_tree= q_full_tree


def main(d_db, infn, opts):
    """
    Run SHOOT
    Args:
        d_db - Input directory name (standard format) containing SHOOT database 
        infn - Query FASTA filename
        opts - Options object

    Returns:
        fn_tree - Filename for tree or None
    """
    if not is_fasta(infn):
        print("ERROR: Input file should be FASTA format")
        return
    # Fix up accession if required
    fn_for_use, _ = clean_fasta(infn)

    # Get genes.
    infile = open(fn_for_use, 'r')
    fasta_data = infile.read()
    genes = dict(re.findall(r"^>([^\n]+)\n([^>]+)", fasta_data, re.MULTILINE))
    infile.close()

    # Assign to tree 
    og_assign = og_assigner.OGAssignDIAMOND(d_db, 
                                            opts.nthreads, 
                                            opts.q_profiles_all)
    # assignations is a dict which keys are gene names and values are dicts
    # with 2 keys:
    #   'ogs': stores the ordered list of assigned OGs
    #   'scores': stores the ordered list of the corresponding e-values
    assignations = og_assign.assign(fn_for_use, q_ultra_sens=opts.search_high_sens)
    
    # orphans is a list of orphan gene names.
    orphans = []
    # aogs are lists of gene names (as values) assigned to OGs identifiers (as keys).
    aogs = {}
    for gene, ga in assignations.items():
        ogs = ga['ogs']
        scores = ga['scores']
        if len(ogs) == 0:
            print("No homologs found for sequence %s in this database" % gene)
            orphans.append(gene)
            continue
        if opts.q_full_tree:
            ogs, scores = utils.og_parts_to_ogs(ogs, scores)
        if len(ogs) > 1:
            print("Assignment of %s to homolog group is ambiguous, %d possibilities:" % (gene, len(ogs)))
            print("Group       E-value")
            for og, score in zip(ogs, scores):
                print("%-11s %g" % (og, score))
            if opts.q_ambiguous:
                print("\nWill attempt placement in each of these trees.")
            else:
                ogs = ogs[:1]
                print("\nWill attempt placement n best hit tree: %s" % ogs[0])
        else:
            # Unambiguous placement
            print("Sequence %s assigned to homolog group OG%s with e-value %g" % (gene, ogs[0], scores[0])) 
        assignations[gene]['ogs'] = ogs
        assignations[gene]['scores'] = scores
        for og in ogs:
            if og not in aogs:
                aogs[og] = []
            aogs[og].append(gene)


    db_name = os.path.basename(d_db[:-1]) 
    with open(infn + ".assign.txt", 'w') as outfile:
        outfile.write("%s\n" % (db_name))

    #for i_tree, (og_part, score) in enumerate(zip(ogs, scores)):
    for i_tree, (og_part, og_genes) in enumerate(aogs.items()):
        og_infn = infn + '.OG' + og_part + '.fa'
        # Place in tree
        warn_str = ""
        if "iqtree" == opts.tree_method:
            graft = msa_grafter.MSAGrafter(d_db, opts.nthreads)
        elif "epa" == opts.tree_method:
            graft = msa_grafter_epa.MSAGrafter_EPA(d_db, opts.nthreads)
        else:
            print("ERROR: %s method has not been implemented" % opts.tree_method)
            return
        # Prepare OG fasta.
        og_fasta = ''
        for og_gene in og_genes:
            og_fasta = og_fasta + '>' + og_gene + '\n' + genes[og_gene] + '\n'
        with open(og_infn, 'w') as outfile:
            outfile.write(og_fasta)
        fn_tree, query_gene_names_final, warn_str = graft.add_gene(
                                        og_part, 
                                        og_infn,
                                        og_infn,
                                        q_mafft_acc=opts.q_mafft_accelerated,
                                        tree_method=opts.tree_method,
                                        )
        scores = []
        for qgene in query_gene_names_final:
            # get the gene score for the selected OG (og_part).
            qogs = assignations[qgene]['ogs']
            qscores = assignations[qgene]['scores']
            qscores = dict(zip(qogs, qscores))
            scores.append(qscores[og_part])

        with open(infn + ".assign.txt", 'a') as outfile:
            outfile.write("%d\t%s\t%s\t%s\n" % (i_tree, og_part, ','.join(map("{:g}".format, scores)), ','.join(query_gene_names_final)))
        if len(ogs) > 1:
            fn_old = fn_tree
            fn_tree = os.path.splitext(fn_tree)[0] + ".%d.tre" % i_tree
            os.rename(fn_old, fn_tree) 
        
        print("Tree: %s" % fn_tree)  
        if warn_str != "":
            print("WARNING: " + warn_str) 

        q_need_to_print = opts.q_print
        if opts.nU is not None:
            # if only nL is specified alone that has no effect
            for og_gene in og_genes:
                t = ete3.Tree(fn_tree)
                if len(t) > opts.nU:
                    node = t & og_gene
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
            fn_ologs = og_infn + ".sh.orthologs.tsv"
            write_orthologs(fn_tree, fn_ologs, og_genes)
 

def is_fasta(infn):
    # rapid check for FASTA format
    try:
        if not os.path.exists(infn):
            return False
        with open(infn, 'r') as infile:
            acc = next(infile)
            ok = bool(re.match("^>.", acc)) 
            ok = ok and bool(re.match("^[a-zA-Z]+$", next(infile))) 
    except:
        return False
    return ok


def clean_fasta(infn):
    """
    If accessions contain problematic characters, write a new file and return
    its filename, otherwise return the original filename
    Args:
        infn - User input filename
    Returns:
        fn_for_use - file to use going forward for SHOOT
    """
    with open(infn, 'r') as infile:
        name_mapping = {}
        fasta_data = infile.read()
        gene_name_matches = re.findall(r"^>([^\n]+)", fasta_data, re.MULTILINE)
        already_clean = True
        for gene_name in gene_name_matches:
            name_cleaned = re.sub(gene_name_disallowed_chars_re, '_', gene_name)
            if gene_name != name_cleaned:
                fasta_data = fasta_data.replace('>' + gene_name, '>' + name_cleaned)
                already_clean = False
            name_mapping[gene_name] = name_cleaned

        if already_clean:
            fn_for_use = infn
        else:
            # otherwise, need to create a new file
            fn_for_use = infn + ".sh.cleaned"
            with open(fn_for_use, 'w') as outfile:
                outfile.write(fasta_data)
            map_file = infn + ".sh.map"
            with open(map_file, 'w') as outfile:
                for gene_name, name_cleaned in name_mapping.items():
                    outfile.write(gene_name + '\t' + name_cleaned + '\n')
        return fn_for_use, name_cleaned


def gene_to_species(name):
    """
    Convert a gene name to species name
    Args:
        name - gene name
    Returns 
        species_name - species name
    """
    return "_".join(name.split("_")[:2])


def write_orthologs(fn_tree, fn_ologs, gene_names):
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
        writer.writerow(["Query", "Species", "Orthologs"])
        for gene_name in gene_names:
            try:
                n = t & gene_name
            except ete3.coretype.tree.TreeError as ex:
                print('Gene %s not found in tree!' % gene_name)
                print(t)
                raise
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
                    add_orthologs_to_file(writer, gene_name, genes_other, overlap)
                n_prev = n
                species_query_branch.update(species_other)


def add_orthologs_to_file(writer, gene_name, genes, overlap):
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
        writer.writerow([gene_name, sp, ", ".join(species_to_genes[sp])])


def clean_dir_name(dname):
    while dname.endswith("/"):
        dname = dname[:-1]
    dname += "/"
    return dname
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help= "Input FASTA filename of the query sequence")
    parser.add_argument("db", help= "Database directory, prepared by create_shoot_db.py")
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
    parser.add_argument("-a", "--ambiquity", action='store_true', 
                        help="Add gene to each tree in case of ambiguous assignment")
    parser.add_argument("-f", "--full", action='store_true', 
                        help="Always add phylogenetically using full tree not super-tree method")
    # threads
    parser.add_argument("-n", "--nthreads", type=int, default=16)
    args, unknown = parser.parse_known_args()
    db = clean_dir_name(args.db)
    opts = Options(
            nthreads=args.nthreads,
            q_profiles_all=args.profiles_all,
            search_high_sens=args.high_sens,
            q_mafft_accelerated=not args.mafft_defaults,
            tree_method=args.tree_method,
            nU=args.upper,
            nL=args.lower, 
            q_print=args.print_tree, 
            q_orthologs=args.orthologs,
            q_ambiguous=args.ambiquity,
            q_full_tree=args.full,
            )
    try:
        main(db, args.infile, opts)
    finally:
        with open(args.infile + ".sh.log", 'w') as outfile:
            outfile.write("Finished\n")


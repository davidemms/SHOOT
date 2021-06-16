#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse

import ete3

import fasta_writer

def write_fasta(args):
    """
    Write the fasta for a target clade
    Args:
        args:
    """
    fn_intree = args.intree
    fn_infasta = args.infasta
    seq = args.seq
    i = args.level
    m = args.max
    # fn_out = args.outfile
    
    try:
        try:
            t = ete3.Tree(fn_intree)
        except ete3.parser.newick.NewickError:
            t = ete3.Tree(fn_intree, format=1)
    except:
        sys.stderr.write("Could not open tree file")
        return
    
    if not os.path.exists(fn_infasta):
        sys.stderr.write("Could not open FASTA file")
        return
    fw = fasta_writer.FastaWriter(fn_infasta)

    if not seq in t:
        sys.stderr.write("Sequence is not in gene tree")
        return

    n = t & seq
    while (i > 0 and n.up is not None):
        if len(n.up) > m:
            break
        n = n.up
        i-=1

    req_genes = n.get_leaf_names()
    # address one common problem. I'll actually solve this higher up the workflow instead though
    if any(r not in fw.SeqLists for r in req_genes):
        fw.SeqLists = {k.replace("|", "_"):v for k, v in fw.SeqLists.items()}
    fw.Print(req_genes)


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Helper utility for SHOOT")
    subparsers = parser.add_subparsers()

    parser_write_fasta = subparsers.add_parser('write_fasta', help="Write fasta file for a clade at specified level above a target gene")
    parser_write_fasta.add_argument("infasta")
    parser_write_fasta.add_argument("intree")
    parser_write_fasta.add_argument("seq", help="target sequence name")
    parser_write_fasta.add_argument("level", type=int, help="level above target sequence for required clade")
    # parser_write_fasta.add_argument("outfile")
    parser_write_fasta.add_argument("-m", "--max", type=int, help="max number of sequences", default=1000)
    parser_write_fasta.set_defaults(func=write_fasta)

    args = parser.parse_args()
    args.func(args)

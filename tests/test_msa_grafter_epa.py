import os
import sys

import pytest
import ete3

d_here = os.path.dirname(__file__) + "/"
sys.path.append(d_here + "..")
d_data = d_here + "data/"

from shoot import msa_grafter_epa as mg

tree167 = "((((Anopheles_gambiae_Q7PY88:1.00341,Schizosaccharomyces_pombe_P87061_TEA1:0.647048)1.0:0.103407,Nematostella_vectensis_A7RUR7:1.03727)1.0:0.259703,(((Oryza_sativa_A0A0N7KFI1:0.023099,Oryza_sativa_A0A0P0VB25:0.21072)1.0:0.122292,Zea_mays_A0A1D6Q3E6:0.187593)1.0:0.253092,(Oryza_sativa_Q53NA8:0.064834,(Zea_mays_A0A1D6F561:0.020896,Zea_mays_C0PDT7:0.031719)1.0:0.057139)1.0:0.377634)1.0:0.518882)0.85:0.441375,(SHOOTOUTGROUP_Leishmania_major_Q4Q9Y2:1.184299,Arabidopsis_thaliana_Q0WW40_FBK5:0.453692):1.253581);\n"



def test_load_single():
    gene, support = mg.MSAGrafter_EPA.get_gene_support(d_data + "epa_result.167.jplace")
    assert gene == "Arabidopsis_thaliana_Q0WW40_FBK5"
    assert support == 1.0
    
    
def test_load_multi():
    gene, support = mg.MSAGrafter_EPA.get_gene_support(d_data + "epa_result.323.jplace")
    assert gene == "Caenorhabditis_elegans_O62337"
    assert support == pytest.approx(0.3490507267)


def test_add_support_out_of_100(tmpdir):
    d = tmpdir.mkdir("subdir")
    fh = d.join("tree323.tre")
    # read-only tree file containing input data, don't edit/over-write
    with open(d_data + "epa_result.323.jplace.epa_result.newick", 'r') as infile:
        fh.write(infile.read())
    fn_tree = os.path.join(fh.dirname, fh.basename)
    mg.MSAGrafter_EPA.add_support(d_data + "epa_result.323.jplace", fn_tree)
        
    t = ete3.Tree(fn_tree)
    n = (t & "Caenorhabditis_elegans_O62337").up
    assert n.support == pytest.approx(34.90507267)


def test_add_support_fractions(tmpdir):
    d = tmpdir.mkdir("subdir")
    fh = d.join("tree167.tre")
    fh.write(tree167)
    fn_tree = os.path.join(fh.dirname, fh.basename)
    mg.MSAGrafter_EPA.add_support(d_data + "epa_result.167.jplace", fn_tree)
        
    t = ete3.Tree(fn_tree)
    n = (t & "Arabidopsis_thaliana_Q0WW40_FBK5").up
    assert n.support == 1.0

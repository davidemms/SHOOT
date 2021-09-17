# SHOOT.bio - the phylogenetic search engine

SHOOT is a phylogenetic alternative to BLAST. Instead of returning a list of similar sequences to a query sequence it returns a maximum likelihood phylogenetic tree with your query sequence embedded in it.

Try it out: https://shoot.bio/

Preprint: https://www.biorxiv.org/content/10.1101/2021.09.01.458564

## Using the SHOOT command line tool
SHOOT allows you to search a protein sequence against a database of gene trees. It returns your gene grafted into the correct position within its corresponding gene tree.

### Preparing a SHOOT phylogenetic database
0. Install dependencies:
    - Python ete3 library
    - DIAMOND
    - MAFFT
    - EPA-ng & gappa (https://github.com/lczech/gappa)
    - Alternatively, IQ-TREE can be used instead of the combination EPA-ng + gappa
2. Run an OrthoFinder analysis on your chosen species, using the multiple sequence alginment option for tree inference, "-M msa".
    - Paper: Emms, D.M., Kelly, S. OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biol 20, 238 (2019). https://doi.org/10.1186/s13059-019-1832-y
    - GitHub: https://github.com/davidemms/OrthoFinder
    - Tutorials: https://davidemms.github.io/
3. Run `python create_shoot_db.py RESULTS_DIRECTORY`, replacing "RESULTS_DIRECTORY" with the path to the OrthoFinder results directory from step 1. 
4. Resolve polytomies (only necessary if using EPA-ng): `python bifurcating_trees.py RESULTS_DIRECTORY`

The OrthoFinder RESULTS_DIRECTORY is now a SHOOT database.

### Running SHOOT
```
python shoot.py INPUT_FASTA SHOOT_DB
```
where INPUT_FASTA is a fasta file containing the amino acid sequence for the search and SHOOT_DB is the SHOOT database directory created using the steps above.

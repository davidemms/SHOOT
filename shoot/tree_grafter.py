import sys

sys.path.append("/home/emms/workspace/p4/trunk/DeepTreePrototype/")


class TreeGrafter(object):
    """
    Places a gene within a gene tree
    """
    def __init__(self, quart, d_db):
        """
        Args:
            quart - Quartets class
            d_db - Database directory
        """
        self.quart = quart
        self.d_db = d_db

    def place_gene(self, iog):
        """
        Main public method for placing a gene in a tree
        Args:
            infn - FASTA filename containing the gene sequence
            iog - the OG containing the gene
        """
        # load tree

        # ask the quartets code what the mapping is from taxa to their indices, 
        # change the names on the trees to the indices.

        # create the tree builder class

        # overwrite its self.t with the current tree

        # call the AddTaxa method

        # DONE
        # after this is working I can make improvements to it. Below are my intial 
        # thoughts on how this code should work:


        # could potentially use the first quartet from the initial diamond assignment, 
        # but that's a potential optimisation not a necessity.


        # decide on first quartet

        # get answer from oracle

        # update location in tree

        # interatively suggest further quartets


    def suggest_first_quartet(self):
        pass

    def suggest__subsequent_quartet(self):
        pass
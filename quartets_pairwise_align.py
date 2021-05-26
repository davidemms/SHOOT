import quartets


class PairwiseAlignQuartets(quartets.Quartets):
    def __init__(self, d_db, iog, infn):
        """
        Args:
            d_db - Database directory
            iog - the OG to search in
            infn - FASTA filename containing the gene sequence
        """
        self.d_db = d_db
        self.iog = iog
        self.infn = infn


    def quartet(self, xyz):
        """
        Return placement of query gene, q, in quartet qxyz
        Args:
            xyz - list of 3 gene IDs for genes x,y & z
        Returns:
            i_sister - sister gene: 0,1 or 2
        """
        raise NotImplemented()


    def ResolveQuartet(self, taxa0):
        """
        Method required by TreeBuilder
        Args:
            taxa - list of taxa indices. First one is the one to be placed
        Returns:
            i    - index of the taxa it should be placed with (1,2 or 3 in base 0)
        """
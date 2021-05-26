import sys

sys.path.append("/home/emms/workspace/p4/trunk/DeepTreePrototype/")

import deeptree


class Quartets(object):
    def __init__(self):
        """
        Abstract base class for quartet oracle classes
        """


class DeepQuartets(Quartets):
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
        self.msa = deeptree.ReadMSA(d_db)
        raise Exception("DeepTree uses an MSA, I don't have this yet")

    def quartet(self, xyz):
        """
        Return placement of query gene, q, in quartet qxyz
        Args:
            xyz - list of 3 gene IDs for genes x,y & z
        Returns:
            i_sister - sister gene: 0,1 or 2
        """
        raise NotImplemented()
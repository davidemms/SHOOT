import sys

# for the time being, import the code I wrote to assign a gene to its orthogroup
sys.path.append("/home/emms/workspace/p4/trunk/OrthoFinder-Accelerate/")
import fol


class OG_Assigner(object):
    """
    Assign a sequence to an orthogroup
    """
    def assign(self, infn):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        raise NotImplemented()


class OG_Assign_D(object):
    """
    Assign using DIAMOND
    """
    def __init__(self, d_db):
        self.d_db = d_db

    def assign(self, infn):
        """
        Assign a sequence to an orthogroup
        Args:
            infn - Input FASTA filename
        Returns
            iog - index of OG
        """
        iog, q_tree = fol.assign_to_orthogroup(infn, infn, self.d_db)

    return iog, q_tree
"""
The database of trees
"""

class Database(object):
    def __init__(self, d_db):
        self.d_db = d_db
        if not self.d_db.endswith("/"):
            self.d_db += "/"
        # Do some checking that everything is present
        # If different options are possible (with/without MSA) then can provide 
        # this info too

    def fn_msa(self, og_part):
        """
        Returns MSA filename
        Args:
            og_part: str, the OG or OG.PART
        """
        if "." in og_part:
            return self.d_db + "Gene_Trees/subtrees/msa_sub/OG%s.fa" % og_part
        else:
            return self.d_db + "MultipleSequenceAlignments/OG%s.fa" % og_part

    def fn_tree(self, og_part):
        """
        Returns tree filename
        Args:
            og_part: str, the OG or OG.PART
        """
        if "." in og_part:
            return self.d_db + "Gene_Trees/subtrees/sub/OG%s.tre" % og_part
        else:
            return self.d_db + "Gene_Trees/OG%s_tree.txt" % og_part
        
        
    def fn_seqs(self, iog):
        return self.d_db + "Orthogroup_Sequences/OG%07d.fa" % iog
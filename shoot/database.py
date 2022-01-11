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

    def fn_tree_super(self, iog):
        """
        Returns tree filename
        Args:
            iog: int
        """
        return self.d_db + "Gene_Trees/subtrees/super/OG%07d.super.tre" % iog

    def fn_trees_sub_without(self, iog, ipart):
        """
        Returns glob pattern for sub trees for this OG
        Args:
            iog: int
        """
        return self.d_db + "Gene_Trees/subtrees/sub/OG%07d.%d.without.tre" % (iog, ipart)

        
    def fn_seqs(self, og_part):
        og = og_part.split(".")[0]
        return self.d_db + "Orthogroup_Sequences/OG%s.fa" % og
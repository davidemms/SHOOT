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

    def fn_msa(self, iog):
        return self.d_db + "MultipleSequenceAlignments/OG%07d.fa" % iog

    def fn_tree(self, iog):
        return self.d_db + "Gene_Trees/OG%07d_tree.txt" % iog
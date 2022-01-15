import sys
import Bio.SeqIO

import database


class Quartets(object):
    def __init__(self, d_db, iog, fn_query):
        """
        Args:
            d_db - Database directory
            iog - the OG to search in
            fn_query - FASTA filename containing the gene sequence
        """
        self.db = database.Database(d_db)
        self.iog = iog
        self.fn_query = fn_query
        fn_og_fasta = self.db.fn_seqs(iog)
        self.seqs, self.taxa_list, self.taxa_lookup = self.read_seqs(fn_query, fn_og_fasta)


    @staticmethod
    def read_seqs(infn):
        """
        Read the sequences in a fasta file
        Args:
            infn - input fasta filename
        Returns:
            seqs - list of biopython records
            taxa_list - ordered list of taxa names
            taxa_lookup - dict from name to index
        """
        seqs = list(Bio.SeqIO.parse(infn, "fasta"))
        taxa_list = [r.name for r in seqs]
        taxa_lookup = {k:i for i,k in enumerate(taxa_list)}
        return seqs, taxa_list, taxa_lookup


# class DeepQuartets(Quartets):
#     def __init__(self, d_db, iog, infn):
#         """
#         Args:
#             d_db - Database directory
#             iog - the OG to search in
#             infn - FASTA filename containing the gene sequence
#         """
#         self.d_db = d_db
#         self.iog = iog
#         self.infn = infn
#         self.msa = deeptree.ReadMSA(d_db)
#         raise Exception("DeepTree uses an MSA, I don't have this yet")

#     def quartet(self, xyz):
#         """
#         Return placement of query gene, q, in quartet qxyz
#         Args:
#             xyz - list of 3 gene IDs for genes x,y & z
#         Returns:
#             i_sister - sister gene: 0,1 or 2
#         """
#         raise NotImplemented()

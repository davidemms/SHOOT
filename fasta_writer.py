# -*- coding: utf-8 -*-

import sys
import gzip
import string
import os.path
import glob
count = 0 

class FastaWriter(object):
    def __init__(self, sourceFastaFilename, qUseOnlySecondPart=False, qGlob=False, qFirstWord=False, qWholeLine=False):
        self.SeqLists = dict()
        qFirst = True
        accession = ""
        sequence = ""
        files = glob.glob(sourceFastaFilename) if qGlob else [sourceFastaFilename]
        for fn in files:
            with (gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn, 'r')) as fastaFile:
                for line in fastaFile:
                    if line[0] == ">":
                        # deal with old sequence
                        if not qFirst:
                            self.SeqLists[accession] = sequence
                            sequence = ""
                        qFirst = False
                        # get id for new sequence
                        if qWholeLine:
                            accession = line[1:].rstrip()
                        else:
                            accession = line[1:].rstrip().split()[0]
                        if qUseOnlySecondPart:
                            try:
                                accession = accession.split("|")[1]
                            except:
                                sys.stderr(line + "\n")
                        elif qFirstWord:
                            accession = accession.split()[0]
                    else:
                        sequence += line
            try:
                self.SeqLists[accession] = sequence
            except:
                sys.stderr(accession + "\n")
    
    def WriteSeqsToFasta(self, seqs, outFilename):
        with open(outFilename, 'w') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % seq)
                    outFile.write(self.SeqLists[seq])
                else:
                    sys.stderr("ERROR: %s not found" % seq)
            
    def Print(self, seqs):
        if type(seqs[0]) is not str:
            seqs = [s.ToString() for s in seqs]
        for seq in seqs:
            if seq in self.SeqLists:
                sys.stdout.write(">%s\n" % seq)
                sys.stdout.write(self.SeqLists[seq])
            else:
                sys.stderr.write("ERROR: %s not found\n" % seq)

                    
    def WriteSeqsToFasta_withNewAccessions(self, seqs, outFilename, idDict):
        with open(outFilename, 'w') as outFile:
            for seq in seqs:
                if seq in self.SeqLists:
                    outFile.write(">%s\n" % idDict[seq])
                    outFile.write(self.SeqLists[seq])
                else:
                    sys.stderr.write(seq + "\n")
    

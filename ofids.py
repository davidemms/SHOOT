# -*- coding: utf-8 -*-

class FullAccession(object):
    def __init__(self, idsFilename):
        # only want the first part and nothing else (easy!)
        self.idToNameDict = dict()
        with open(idsFilename, 'r') as idsFile:
            for line in idsFile:
                line = line.rstrip()
                if not line: continue
#                if line.startswith("#"): continue
                id, accession = line.split(": ", 1)
                id = id.replace("#", "")
                id = id.strip()
                # Replace problematic characters
                accession = accession.replace(":", "_").replace(",", "_").replace("(", "_").replace(")", "_") #.replace(".", "_")
                if id in self.idToNameDict:
                    raise RuntimeError("ERROR: A duplicate id was found in the fasta files: % s" % id)
                self.idToNameDict[id] = accession 
                
    def GetIDToNameDict(self):
        return self.idToNameDict
        
        
class OrthoFinderIDs(object):
    def __init__(self, wd, idExtractor = FullAccession):
        self.wd = wd
        self.speciesIDsEx = FullAccession(wd + "/SpeciesIDs.txt")
        self.species_dict = self.SpeciesDict()
        self.nsp = len(self.species_dict)
        self._Spec_SeqIDs = None
        self._extractor = idExtractor
        self.seqIDsEx = None
        self.qAddSpeciesToIDs = True

    def SequenceDict(self):
        if self.seqIDsEx == None:
            try:
                self.seqIDsEx = self._extractor(self.wd + "/SequenceIDs.txt")
            except RuntimeError as error:
                print(str(error))
                print("Tried to use only the first part of the accession in order to list the sequences in each orthogroup\nmore concisely but these were not unique. The full accession line will be used instead.\n")     
                self.seqIDsEx = FullAccession(self.wd + "/SequenceIDs.txt")
        return self.seqIDsEx.GetIDToNameDict()
        
    def SpeciesDict(self):
        d = self.speciesIDsEx.GetIDToNameDict()
        return {k:v.rsplit(".",1)[0] for k,v in d.items()}
        
    def Spec_SeqDict(self):
        if self._Spec_SeqIDs != None:
            return self._Spec_SeqIDs
        seqs = self.SequenceDict()
        seqs = {k:v for k,v in seqs.items() if int(k.split("_")[0]) in range(self.nsp)}
        if not self.qAddSpeciesToIDs:
            self._Spec_SeqIDs = seqs
            return seqs
        specs = self.species_dict
        specs_ed = {k:v.replace(".", "_").replace(" ", "_") for k,v in specs.items()}
        self._Spec_SeqIDs = {seqID:specs_ed[seqID.split("_")[0]] + "_" + name for seqID, name in seqs.items()}
        return self._Spec_SeqIDs
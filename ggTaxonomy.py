''' read in gg_13_5.fasta (downloaded from http://greengenes.secondgenome.com/downloads)
and gg_13_5_taxonomy.txt and output Fasta of greengenes reads with phyla assigned'''
import re
taxonomyDict = {} # ID and then taxonomy to genus level (I just need genus to test correctness with BLAST)

def readInTaxa(file):
    '''creates a dictionary of each sequence ID as keys and the 7-level taxonomy as values'''
    with open(file, 'r') as f:
        for line in f:
            split_line = line.split('\t')
            ID = split_line[0]
            taxa = split_line[1].strip('\n')
            taxonomyDict[ID] = taxa

def assignTaxonomy(file):
    '''creates a new fasta file with each ID and taxonomy assigned'''
    assignedOutput = open("gg_13_5_taxonomy.fasta", 'w')
    with open(file, 'r') as fastaFile:
        missed = False
        for line in fastaFile:
            if line[0] == ">":
                ID = line[1:].strip('\n')
                taxonomy = taxonomyDict[ID]
                assignedOutput.write('>' + ID + '\t' + taxonomy + '\n') 
            else:
                assignedOutput.write(line)


readInTaxa("gg_13_5_taxonomy.txt")
assignTaxonomy("gg_13_5.fasta")

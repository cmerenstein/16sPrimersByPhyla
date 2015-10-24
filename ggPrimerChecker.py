#Carter Merenstein
#Middlebury College
''' assigns sequences to primers. This method doesn't rely on alignment or 
prior knowledge of where primers are in the 16S sequence. This is an advantage,
as sequences in GreenGenes may not all align exactly to the E. coli 16S in length'''

import csv
import time

primers_by_name = {} #final output will be by name, not sequence
primers_as_keys = {} #primer sequences in keys. Values will be dict with phylum:number format
phyla = {} #total number of each phyla

min_k = 999        #smallest primer size
max_k = 0       #largest

def readInPrimers(file):
    ''' reads file of primers and builds dict with name as key and primer options as values in a list'''
    temp_pbn = {} #will be returned to set primers_by_name
    with open(file, 'r') as primer_file:
        for line in primer_file:
            line = line.strip("\n")
            line_list = line.split(',')
            name = line_list[0]
            primers = line_list[1:]
            temp_pbn[name] = primers
    return temp_pbn

def buildPAK(dictIn):
    ''' takes primers by name and builds dict of each primer, not associated with a name.
    separate functions it makes it easier to quickly see the difference between the two dictionaries'''
    temp_pak = {}
    for name in primers_by_name.keys():         # keys are primer names
        for value in primers_by_name[name]:     # value is the primer sequence
            temp_pak[value] = {} #dict will be phylum:number
    return temp_pak

def getMinPrimer(primers):
    ''' later we'll need to know the minimum word size to look for primers'''
    temp_min = 999
    for primer in primers:
        if len(primer) < temp_min:
            temp_min = len(primer)
    return temp_min

def getMaxPrimer(primers):
    ''' later we'll need to know the maximum word size to look for primers'''
    temp_max = 0
    for primer in primers:
        if len(primer) > temp_max:
            temp_max = len(primer)
    return temp_max

def findPrimers(file, pak):
    ''' takes in file of sequences and dictionary with the primers as keys.'''
    
    with open(file, 'r') as fastaFile:
        phylum = ''
        ID = ''
        for line in fastaFile:
            if line[0] == '>':
                ID = line.split('\t')[0].strip('>') #get the ID
                phylum = line.split(';')[1].strip(' ') #get the phylum from the whole taxonomy
                try:
                    phyla[phylum] +=1  #count the total times a phylum occurs
                except:
                    phyla[phylum] = 1
            else:
                '''iterate through sequence checking each word of size min_k to max_k. Check if it is in the primers_as_keys
    dictionary. The values are dictionaries of phylum[number], increment number'''
                i = 0 # start of word
                j = 0 #used to check different sized primers
                for j in range((max_k +1 )- min_k):   ##this lets us go through the whole sequence with a word size from min to max, to get every primer
                    k = min_k + j # end of word
                    for l in range(len(line)-k + 1): #this actually goes through the whole sequence
                        word = line[(i+l):(k+l)].lower()
                        try:
                            pak[word] #check if it matches a primer
                            try:
                                pak[word][phylum] += 1 #increment as an occurrence of that phylum
                            except:
                                 pak[word][phylum] = 1  #if the phylum hasn't been found yet, but the word matches a primer
                        except:
                            pass # if not found, it's not a primer. move on
    return pak

def testOutput(primerNames, pak):
	''' only use this when testing on a truncated fasta file. Will take forever with whole GreenGenes'''
    for primer in primerNames.keys():
        print(primer)
        for seq in primerNames[primer]:
            print(pak[seq])
            
def outputPrimersTaxonomy(primers_with_taxonomy, primer_varients, total_phyla):
    ''' primers_with_taxonomy is equivalent to primers_as_keys, primer_varients is primers_by_name
       Will output CSV with primers (named) in rows with phyla as columns and occurrence (%) as values'''
    with open("ggPrimerPhylaSkipped.csv", 'w', newline = '') as out:
        w = csv.writer(out)
        phyla_list = []
        header = ['Phyla:']
        for phylum in total_phyla.keys():
            phyla_list.append(phylum)
            header.append(phylum.strip('p_'))
        w.writerow(header) #Names of phyla
        n = ['n'] #row to count occurrence of phyla
        for phylum in phyla_list:
            n.append(total_phyla[phylum])
        w.writerow(n)
        for primer in primer_varients.keys():
            line = [primer]
            for phylum in phyla_list:
                count = 0
                for varient in primer_varients[primer]:
                    try:
                        count += primers_with_taxonomy[varient][phylum]
                    except KeyError:
                        pass
                coverage = count/total_phyla[phylum]
                line.append(coverage)
            w.writerow(line)
        
                
primers_by_name = readInPrimers("universal_primers_R.txt") #reverse primers were reverse complimented
primers_as_keys = buildPAK(primers_by_name)
min_k = getMinPrimer(primers_as_keys.keys())
max_k = getMaxPrimer(primers_as_keys.keys())
timer = time.time()
primers_as_keys = findPrimers('gg_13_5_taxonomy.fasta', primers_as_keys)
now = time.time()
elapsed  = now - timer

outputPrimersTaxonomy(primers_as_keys, primers_by_name, phyla)
##testOutput(primers_by_name, primers_as_keys)


print(elapsed)

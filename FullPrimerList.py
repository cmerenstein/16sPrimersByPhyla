''' takes a list of primers and names in IUPAC base form and outputs a list of all the primers
with each iteration written out in only standard bases.
Reverse primers have the reverse compliment taken, (and therefore end up
corresponding to the 5'->3' of the GreenGenes database)'''
from Bio.Seq import Seq #used to do reverse complement

primers = {}
primers_iupac = {}

iupac_bases = {'a':['a'], 'c':['c'], 't':['t'],'g':['g'],'r': ['a','g'], 'y':['c','t'],'s':['g','c'],'w':['a','t'],'k':['g','t'],'m':['a','c'],'b':['c','g','t'],'d':['a','g','t'],'h':['a','c','t'],'v':['a','c','g'],'n':['a','c','g','t']}
single_bases = ['a', 'c', 't', 'g']

def revComp(string):
    ''' returns reverse compliment, used for reverse primers '''
    my_seq = Seq(string)
    my_seq = my_seq.reverse_complement()
    return str(my_seq)

with open("universal_primers.txt", 'r') as primer_txt:
    ''' in order to later check denegrated primers using exact matching
    we turn the denegrated forms into a list of all possible options'''
    for line in primer_txt:
        line_arr = line.split(' ')
        primers[line_arr[0]] = []   #each primer name starts with a blank list
        primers_iupac[line_arr[1].lower()] = line_arr[0] #iupac bases for primer with name as value
    for primer in primers_iupac.keys():
  
        complete_list = ['']
        for base in primer:
            new_list = []
            for item in complete_list:
                for option in iupac_bases[base]:
                    new_item = item + option
                    new_list.append(new_item)
            complete_list = new_list
            primers[primers_iupac[primer]] = complete_list

with open("universal_primers_complete_R.txt", 'w') as output:
    ''' writes output as name first, then primers, all separated by commas'''
    for primer in primers.keys():
        output.write(primer)
        if (primer[-1] == 'R' or primer[-2] == 'R'): # primers labelled F or R for direction
            for seq in primers[primer]:
                output.write(',' + revComp(seq))
        else:
            for seq in primers[primer]:
                output.write(',' + seq)
        output.write("\n")


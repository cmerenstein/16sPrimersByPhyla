Carter Merenstein
10/10/2015

--- All primers are from Seorgel et al. ---
universal_primers.txt contains all primers and initial source.
universal_primers_complete_R.txt contains a translated list wherein each primer with
ambiguous bases is translated into multiple primers made up of only standard bases.
Reverse primers are given with their reverse compliment in order to match greengeens
direction.
This is done by FullPrimerList.py
accuracy of FullPrimerList is checked by u_p_c_check.py

--- GreenGenes data is from most recent (May 2013) release ---
http://greengenes.secondgenome.com/downloads
16s sequences with IDs is available as gg_13_5.fasta.gz
Taxonomy associated with each ID is gg_13_5_taxonomy.txt.gz

--- Fasta file with ID and taxonomy is created by ggTaxonomy.py ---
This produces the file gg_13_5_taxonomy.fasta, which is used for primer matching.
Proper matching of taxonomy will be tested with BLAST on some n sequences

--- ggPrimerChecker.py reads in primer file and fasta ---
It uses exact matching to check which primers are present in each sequence.
Takes about 3 hours to run on the full GreenGenes database
from bioio import fastqRead, fastaRead, system, reverseComplement
import os, sys, itertools
from collections import Counter
from math import log

"""Runs kmer analysis"""

def countKmers(kmerSize, referenceFastaFile, readFastqFile):
	refKmers, readKmers = Counter(), Counter()

	for name, seq in fastaRead(referenceFastaFile):
	    for i in xrange(kmerSize, len(seq)):
		s = seq[ i - kmerSize : i ]
		if "N" not in s:
		    refKmers[s] += 1
                    #refKmers[reverseComplement(s)] += 1


	for name, seq, qual in fastqRead(readFastqFile):
	    for i in xrange(kmerSize, len(seq)):
		s = seq[ i - kmerSize : i ]
		if "N" not in s:
		    readKmers[s] += 1
                    #readKmers[reverseComplement(s)] += 1

	return (refKmers, readKmers)

def analyzeCounts(refKmers, readKmers, name):
	refSize, readSize = sum(refKmers.values()), sum(readKmers.values())
	outf = open(name + "kmer_counts.txt", "w")
	outf.write("kmer\trefCount\trefFraction\treadCount\treadFraction\tlogFoldChange\n")

	for kmer in itertools.product("ATGC", repeat=5):
	    kmer = "".join(kmer)
	    refFraction, readFraction = 1.0 * refKmers[kmer] / refSize, 1.0 * readKmers[kmer] / readSize
	    if refFraction == 0:
		foldChange = "-Inf"
	    elif readFraction == 0:
		foldChange = "Inf"
	    else:
		foldChange = -log(readFraction / refFraction)
	    outf.write("\t".join(map(str,[kmer, refKmers[kmer], refFraction, readKmers[kmer], readFraction, foldChange]))+"\n")
	outf.close()

	#system("Rscript ./kmer_analysis.R {} {} {} {} {}".format(name + "kmer_counts.txt", name + "pval_kmer_counts.txt", name + "top_bot_sigkmer_counts.txt", name + "volcano_plot.pdf", "Kmer"))

# args
referenceFastaFile = sys.argv[1]
readFastqFile = sys.argv[2]
file_name = str(sys.argv[3])
kmerSize = 5 #sys.argv[4]

#analyze kmers across both files
refKmers, readKmers = countKmers(kmerSize, referenceFastaFile, readFastqFile)
if len(refKmers) > 0 and len(readKmers) > 0:
    analyzeCounts(refKmers, readKmers, file_name + "_all_bases_")



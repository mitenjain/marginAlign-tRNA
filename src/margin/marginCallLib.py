import pysam, sys, os, collections
from jobTree.src.bioio import fastaRead, system, fastaWrite, logger
import numpy as np
from margin.utils import *
try:
    import cPickle 
except ImportError:
    import pickle as cPickle

def marginCallTargetFn(target, samFile, referenceFastaFile, outputVcfFile,  options):
    """Calculates the posterior probabilities of all the matches in between
    each pairwise alignment of a read to the reference. The collates this 
    posterior probabilities and uses them to call SNVs.
    """
    target.setFollowOnTargetFn(paralleliseSamProcessingTargetFn, 
                               args=(samFile, 
                                     readFastqFile, referenceFastaFile, outputVcfFile, 
                                     posteriorProbabilityCalculationTargetFn, 
                                     variantCallSamFileTargetFn, options))
   
def posteriorProbabilityCalculationTargetFn(target, exonerateCigarStringFile, 
                referenceSequenceName, referenceSequence, querySequenceFile, 
                outputPosteriorProbsFile, options):
    """Calculates the posterior probabilities of matches in a set of pairwise
    alignments between a reference sequence and a set of reads. 
    """
    #Temporary files
    tempRefFile = os.path.join(target.getLocalTempDir(), "ref.fa")
    tempReadFile = os.path.join(target.getLocalTempDir(), "read.fa")
    
    #Write the temporary reference file.
    fastaWrite(tempRefFile, referenceSequenceName, referenceSequence) 
    
    #Hash to store posterior probabilities in
    expectationsOfBasesAtEachPosition = {}
    
    #For each cigar string
    for exonerateCigarString, (querySequenceName, querySequence) in \
    zip(open(exonerateCigarStringFile, "r"), fastaRead(querySequenceFile)):
        fastaWrite(tempReadFile, querySequenceName, querySequence)
        #Call to cactus_realign
        tempPosteriorProbsFile = os.path.join(target.getLocalTempDir(), "posteriorProbs.txt")
        system("echo %s | cactus_realign %s %s --diagonalExpansion=10 \
        --splitMatrixBiggerThanThis=100 --outputAllPosteriorProbs=%s --loadHmm=%s" % \
                   (exonerateCigarString[:-1], tempRefFile, tempReadFile, 
                    tempPosteriorProbsFile, options.hmmFile))
        
        #Now collate the reference position expectations
        for refPosition, readPosition, posteriorProb in \
        map(lambda x : map(float, x.split()), open(tempPosteriorProbsFile, 'r')):
            key = (referenceSequenceName, int(refPosition))
            if key not in expectationsOfBasesAtEachPosition:
                expectationsOfBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases)))
            readBase = readSeq[int(readPosition)].upper()
            if readBase in bases:
                expectationsOfBasesAtEachPosition[key][readBase] += posteriorProb
    
    #Pickle the posterior probs
    fileHandle = open(outputPosteriorProbsFile, 'w')
    cPickle.dump(expectationsOfBasesAtEachPosition, fileHandle, cPickle.HIGHEST_PROTOCOL)
    fileHandle.close() 

def variantCallSamFileTargetFn(target, samFile, referenceFastaFile, 
                            outputVcfFile, tempPosteriorProbFiles, options):
    """Collates the posterior probabilities and calls SNVs for each reference base.
    """
    #Hash to store posterior probabilities in
    expectationsOfBasesAtEachPosition = {}

    #Read in the posterior probs from the pickles from the tempPosteriorProbFiles files
    for tempPosteriorProbFile in tempPosteriorProbFiles:
        fileHandle = open(tempPosteriorProbFile, 'r')
        expectationsOfBasesAtEachPosition2 = cPickle.load(fileHandle)
        for key in expectationsOfBasesAtEachPosition2:
            if key not in expectationsOfBasesAtEachPosition2:
                expectationsOfBasesAtEachPosition[key] = dict(zip(bases, [0.0]*len(bases)))
            for base in bases:
                expectationsOfBasesAtEachPosition[key][base] += expectationsOfBasesAtEachPosition2[key][base]
        fileHandle.close()
    
    #Array to store the VCF calculations
    variantCall = [] #Each key is of the form (referenceName, referencePosition, baseCall)
    
    #Now do the SNV calculations
    
    #For each call write out a VCF line representing the output.
    pass

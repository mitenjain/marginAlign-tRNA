import pysam, sys, os, collections
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, cigarReadFromString, PairwiseAlignment, system, fastaWrite, fastqWrite, cigarRead, logger, nameValue, absSymPath
from cactus.bar import cactus_expectationMaximisation
from cactus.bar.cactus_expectationMaximisation import Hmm, SYMBOL_NUMBER
import numpy as np

def variantCallSamFileTargetFn(target, samFile, vcfFile, options):
    """Calls SNVs using the margin-call method for the given SAM file.
    """
    #Chain the sam file
    tempSamFile = os.path.join(target.getGlobalTempDir(), "temp.sam")
    chainSamFile(samFile, tempSamFile, readFastqFile, referenceFastaFile, chainFn)
    
    #If we do expectation maximisation we split here:
    if options.em:
        target.addChildTargetFn(learnModelFromSamFileTargetFn, 
                                args=(tempSamFile, readFastqFile, 
                                      referenceFastaFile, options))

    target.setFollowOnTargetFn(realignSamFile2TargetFn, args=(tempSamFile, 
              outputSamFile, readFastqFile, referenceFastaFile, options))


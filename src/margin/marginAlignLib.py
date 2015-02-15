import pysam, sys, os
from jobTree.src.bioio import reverseComplement, fastaRead, system, fastaWrite, \
cigarRead, logger, nameValue
from margin.utils import *
from cPecan import cPecanEm
from cPecan.cPecanEm import Hmm, SYMBOL_NUMBER
import numpy as np

def mergeChainedAlignedReads(chainedAlignedReads, refSequence, readSequence):
    """Makes a global alignment for the given chained reads.
    
    From doc on building pysam line
    a = pysam.AlignedRead()
    a.qname = "read_28833_29006_6945"
    a.seq="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.flag = 99
    a.rname = 0
    a.pos = 32
    a.mapq = 20
    a.cigar = ( (0,10), (2,1), (0,25) )
    a.mrnm = 0
    a.mpos=199
    a.isize=167
    a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
    a.tags = ( ("NM", 1),
               ("RG", "L1") )
    """
    cAR = pysam.AlignedRead()
    aR = chainedAlignedReads[0]
    cAR.qname = aR.qname
    
    #Parameters we don't and therefore set properly
    #cAR.flag = aR.flag
    #cAR.mapq = aR.mapq
    #cAR.mrnm = 0
    #cAR.mpos=0
    #cAR.isize=0
    #cAR.qual = "<" * len(readSequence)
    #cAR.tags = aR.tags 
    cAR.rnext = -1
    cAR.pos = 0
    cAR.is_reverse = aR.is_reverse
    if cAR.is_reverse:
        cAR.seq = reverseComplement(readSequence)
    else:
        cAR.seq = readSequence
    cAR.rname = aR.rname
    cigarList = []
    pPos = 0
    if cAR.is_reverse: #Iterate from the other end of the sequence
        pQPos = -(len(readSequence)-1)
    else:
        pQPos = 0
        
    for aR in chainedAlignedReads:
        assert cAR.is_reverse == aR.is_reverse
        #Add a deletion representing the preceding unaligned reference positions
        assert aR.pos >= pPos
        if aR.pos > pPos:
            cigarList.append((2, aR.pos - pPos))
            pPos = aR.pos 
    
        #Add an insertion representing the preceding unaligned read positions
        qPos = getAbsoluteReadOffset(aR, refSequence, readSequence)
        assert qPos >= pQPos
        if qPos > pQPos:
            cigarList.append((1, qPos - pQPos)) 
            pQPos = qPos
        
        #Add the operations of the cigar, filtering hard and soft clipping
        for op, length in aR.cigar:
            assert op in (0, 1, 2, 4, 5)
            if op in (0, 1, 2):
                cigarList.append((op, length))
            if op in (0, 2): #Is match or deletion
                pPos += length
            if op in (0, 1): #Is match or insertion
                pQPos += length
        
    #Now add any trailing deletions/insertions
    assert pPos <= len(refSequence)
    if pPos < len(refSequence):
        cigarList.append((2, len(refSequence) - pPos))
    
    if cAR.is_reverse:
        assert pQPos <= 1
        if pQPos < 1:
            cigarList.append((1, -pQPos + 1))
    else:
        assert pQPos <= len(readSequence)
        if pQPos < len(readSequence):
            cigarList.append((1, len(readSequence) - pQPos))
    
    #Check coordinates
    #print cAR.is_reverse, sum([ length for op, length in cigarList if op in (0, 2)]),  
    #len(refSequence), sum([ length for op, length in cigarList if op in (0, 1)]), 
    #len(readSequence), cAR.qname
    assert sum([ length for op, length in cigarList if op in (0, 2)]) == len(refSequence)
    assert sum([ length for op, length in cigarList if op in (0, 1)]) == len(readSequence)
    
    cAR.cigar = tuple(cigarList)
    
    return cAR

def chainFn(alignedReads, refSeq, readSeq, scoreFn=\
            lambda alignedRead, refSeq, readSeq : \
            len(list(AlignedPair.iterator(alignedRead, refSeq, readSeq))), maxGap=200):
    """Gets the highest scoring chain of alignments on either the forward or reverse 
    strand. Score is (by default) number of aligned positions.
    """
    def getStartAndEndCoordinates(alignedRead):
        """Gets the start and end coordinates in both the reference and query
        """
        alignedPairs = list(AlignedPair.iterator(alignedRead, refSeq, readSeq))
        return alignedPairs[0].refPos, alignedPairs[0].getSignedReadPos(), \
    alignedPairs[-1].refPos, alignedPairs[-1].getSignedReadPos()
    alignedReadToScores = dict([ (aR, scoreFn(aR, refSeq, readSeq)) for aR in alignedReads])
    alignedReadToCoordinates = dict([ (aR, getStartAndEndCoordinates(aR)) for \
                                     aR in alignedReads])
    alignedReadPointers = {}
    
    #Currently uses sloppy quadratic algorithm to find highest chain
    alignedReads = sorted(alignedReads, key=lambda aR : alignedReadToCoordinates[aR][0]) 
    #Sort by reference coordinate
    for i in xrange(len(alignedReads)):
        aR = alignedReads[i]
        rStart, qStart, rEnd, qEnd = alignedReadToCoordinates[aR]
        score = alignedReadToScores[aR]
        for j in xrange(i): #Look at earlier alignments in list
            aR2 = alignedReads[j]
            rStart2, qStart2, rEnd2, qEnd2 = alignedReadToCoordinates[aR2]
            assert rStart2 <= rStart
            if rStart > rEnd2 and qStart > qEnd2 and aR.is_reverse == aR2.is_reverse and \
            rStart - rEnd2 + qStart - qEnd2 <= maxGap and \
            score + alignedReadToScores[aR2] > alignedReadToScores[aR]: 
            #Conditions for a chain
                alignedReadToScores[aR] = score + alignedReadToScores[aR2]
                alignedReadPointers[aR] = aR2
    
    #Now find highest scoring alignment
    aR = sorted(alignedReads, key=lambda aR : alignedReadToScores[aR])[-1]
    
    #Construct chain of alignedReads
    chain = [ aR ]
    while aR in alignedReadPointers:
        aR = alignedReadPointers[aR]
        chain.append(aR)
    chain.reverse()
    
    return chain

def chainSamFile(samFile, outputSamFile, readFastqFile, referenceFastaFile, 
                 chainFn=chainFn):
    """Chains together the reads in the SAM file so that each read is covered by a 
    single maximal alignment
    """
    sam = pysam.Samfile(samFile, "r" )
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    readSequences = getFastqDictionary(readFastqFile) #Hash of names to sequences
    readsToAlignedReads = {}
    for aR in samIterator(sam): #Iterate on the sam lines and put into buckets by read
        if aR.qname not in readSequences:
            raise RuntimeError("Aligned read name: %s not in read sequences \
            names: %s" % (aR.qname, readSequences.keys()))
        key = (aR.qname,aR.rname)
        if key not in readsToAlignedReads:
            readsToAlignedReads[key] = []
        readsToAlignedReads[key].append(aR)
    #Now write out the sam file
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    
    #Chain together the reads
    chainedAlignedReads = []
    for readName, refID in readsToAlignedReads.keys():
        alignedReads = readsToAlignedReads[(readName, refID)]
        refSeq = refSequences[sam.getrname(refID)]
        readSeq = readSequences[readName]
        chainedAlignedReads.append(mergeChainedAlignedReads(chainFn(alignedReads, 
                                                refSeq, readSeq), refSeq, readSeq))
    chainedAlignedReads.sort() #Sort by reference coordinate
    for cAR in chainedAlignedReads:
        outputSam.write(cAR)
    sam.close()
    outputSam.close()
    
def learnModelFromSamFileTargetFn(target, samFile, readFastqFile, 
                                  referenceFastaFile, options):
    """Does expectation maximisation on sam file to learn the hmm for the sam file.
    """
    #Convert the read file to fasta
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    readSequences = getFastqDictionary(readFastqFile) #Hash of names to sequences
    
    reads = os.path.join(target.getGlobalTempDir(), "temp.fa")
    fH = open(reads, 'w')
    for name in readSequences.keys():
        seq = readSequences[name]
        fastaWrite(fH, name, seq)
        fastaWrite(fH, name + "_reverse", reverseComplement(seq))
    fH.close()
    
    #Get cigars file
    cigars = os.path.join(target.getGlobalTempDir(), "temp.cigar")
    fH = open(cigars, 'w')
    sam = pysam.Samfile(samFile, "r" )
    for aR in sam: #Iterate on the sam lines realigning them in parallel            
        #Because these are global alignments with reverse complement coordinates 
        #reversed the following should all be true
        assert aR.pos == 0
        assert aR.qstart == 0
        assert aR.qend == len(readSequences[aR.qname]) #aR.query)
        assert aR.aend == len(refSequences[sam.getrname(aR.rname)])
        assert len(aR.query) == len(readSequences[aR.qname])
        if aR.is_reverse: #Deal with reverse complements
            assert aR.query.upper() == reverseComplement(readSequences[aR.qname]).upper()
            aR.qname += "_reverse"
        else:
            assert aR.query.upper() == readSequences[aR.qname].upper()
            
        fH.write(getExonerateCigarFormatString(aR, sam) + "\n")
        #Exonerate format Cigar string, using global coordinates
        #fH.write(getGlobalAlignmentExonerateCigarFormatString(aR, sam, 
        #refSequences[sam.getrname(aR.rname)], readSequences[aR.qname]) + "\n")
    fH.close()
    
    unnormalisedOutputModel = os.path.join(target.getGlobalTempDir(), 
                                           "unnormalisedOutputModel.hmm")
    target.addChildTargetFn(cPecanEm.expectationMaximisationTrials, 
                            args=(" ".join([reads, referenceFastaFile ]), cigars, 
                                  unnormalisedOutputModel, options))
    
    #Now set up normalisation
    target.setFollowOnTargetFn(learnModelFromSamFileTargetFn2, 
                               args=(unnormalisedOutputModel, options))

def learnModelFromSamFileTargetFn2(target, unnormalisedOutputModel, options):
    hmm = Hmm.loadHmm(unnormalisedOutputModel)
    setHmmIndelEmissionsToBeFlat(hmm)
    #Normalise background emission frequencies, if requested to GC% given
    normaliseHmmByReferenceGCContent(hmm, 0.5)
    hmm.write(options.outputModel)
    
toMatrix = lambda e : map(lambda i : e[SYMBOL_NUMBER*i:SYMBOL_NUMBER*(i+1)], 
                          xrange(SYMBOL_NUMBER))
fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)
    
def normaliseHmmByReferenceGCContent(hmm, gcContent):
    """Normalise background emission frequencies to GC% given
    """
    for state in range(hmm.stateNumber):
        if state not in (2, 4): #Don't normalise GC content of insert states 
            #(as they don't have any ref bases!)
            n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * 
                                       state:(SYMBOL_NUMBER**2) * (state+1)])
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
            fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * 
            (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), 
                           range(SYMBOL_NUMBER))) #Normalise

def setHmmIndelEmissionsToBeFlat(hmm):
    """Set indel emissions to all be flat
    """
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
        [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2  

def modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate):
    #Normalise background emission frequencies, if requested to GC% given
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == \
                     i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), 
                     xrange(SYMBOL_NUMBER**2)))
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

def realignSamFileTargetFn(target, samFile, outputSamFile, readFastqFile, 
                           referenceFastaFile, options, chainFn=chainFn):
    """Chains and then realigns the resulting global alignments, using jobTree to 
    do it in parallel on a cluster.
    Optionally runs expectation maximisation.
    """
    #Chain the sam file
    tempSamFile = os.path.join(target.getGlobalTempDir(), "temp.sam")
    chainSamFile(samFile, tempSamFile, readFastqFile, referenceFastaFile, chainFn)
    
    #If we do expectation maximisation we split here:
    if options.em:
        target.addChildTargetFn(learnModelFromSamFileTargetFn, args=(tempSamFile, 
                                    readFastqFile, referenceFastaFile, options))

    #target.setFollowOnTargetFn(realignSamFile2TargetFn, args=(tempSamFile, outputSamFile, 
    #                                         readFastqFile, referenceFastaFile, options))
    options.hmmFile = options.outputModel if options.em else options.inputModel #This
    #setups the hmm to be used the realignment function
    
    target.setFollowOnTargetFn(paralleliseSamProcessingTargetFn, 
                               args=(tempSamFile, 
                                     referenceFastaFile, outputSamFile, 
                                     realignCigarTargetFn, realignSamFile3TargetFn,
                                     options))
    
def realignCigarTargetFn(target, exonerateCigarStringFile, referenceSequenceName, 
                         referenceSequence, querySequenceFile, 
                         outputCigarFile, options):
    #Temporary files
    tempRefFile = os.path.join(target.getLocalTempDir(), "ref.fa")
    tempReadFile = os.path.join(target.getLocalTempDir(), "read.fa")
    
    #Write the temporary reference file.
    fastaWrite(tempRefFile, referenceSequenceName, referenceSequence) 
    
    #For each cigar string
    for exonerateCigarString, (querySequenceName, querySequence) in \
    zip(open(exonerateCigarStringFile, "r"), fastaRead(querySequenceFile)):
        fastaWrite(tempReadFile, querySequenceName, querySequence)
        #Call to cPecanRealign
        loadHmm = nameValue("loadHmm", options.hmmFile)
        system("echo %s | cPecanRealign %s %s --diagonalExpansion=10 \
        --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s --matchGamma=%s >> %s" % \
               (exonerateCigarString[:-1], tempRefFile, tempReadFile, loadHmm, 
                options.gapGamma, options.matchGamma, outputCigarFile))

def realignSamFile3TargetFn(target, samFile, referenceFastaFile, 
                            outputSamFile, tempCigarFiles, options):
    #Setup input and output sam files
    sam = pysam.Samfile(samFile, "r" )
    
    #Replace the cigar lines with the realigned cigar lines
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    def cigarIterator():
        #Iterates over all the cigars in the temp files.
        for tempCigarFile in tempCigarFiles:
            for pA in cigarRead(open(tempCigarFile)):
                yield pA 
        yield None #This is put in to cause an error if there is fewer 
        #cigars than pairwise alignments
    for aR, pA in zip(samIterator(sam), cigarIterator()): #Iterate on the sam lines 
        #realigning them in parallel
        
        #Convert to sam line
        aR.cigar = tuple([ (op.type, op.length) for op in pA.operationList ])
        
        #Write out
        outputSam.write(aR)
    
    #Finish up
    sam.close()
    outputSam.close()
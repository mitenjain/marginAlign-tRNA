import pysam, sys, os
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, \
cigarReadFromString,PairwiseAlignment, fastaWrite, fastqWrite, logger, absSymPath

def pathToBaseNanoporeDir():
    """Returns path to base directory "marginAlign"
    """
    import marginAlign
    i = absSymPath(__file__)
    return os.path.split(os.path.split(i)[0])[0]

class AlignedPair:
    """Represents an aligned pair of positions.
    """
    def __init__(self, refPos, refSeq, readPos, isReversed, readSeq, pPair):
        assert refPos >= 0 and refPos < len(refSeq)
        self.refPos = refPos
        self.refSeq = refSeq
        assert readPos >= 0 and readPos < len(readSeq)
        self.readPos = readPos
        self.isReversed = isReversed
        self.readSeq = readSeq
        self.pPair = pPair #Pointer to the previous aligned pair
        
    def isMatch(self):
        return self.getRefBase().upper() == self.getReadBase().upper() and \
            self.getRefBase().upper() in "ACTG"
    
    def isMismatch(self):
        return self.getRefBase().upper() != self.getReadBase().upper() and \
            self.getRefBase().upper() in "ACTG" and self.getReadBase().upper() in "ACTG"
    
    def getRefBase(self):
        return self.refSeq[self.refPos]
    
    def getReadBase(self):
        if self.isReversed:
            return reverseComplement(self.readSeq[self.readPos]) 
        return self.readSeq[self.readPos]
    
    def getSignedReadPos(self):
        if self.isReversed:
            return -self.readPos
        return self.readPos
    
    def getPrecedingReadInsertionLength(self, globalAlignment=False):
        if self.pPair == None:
            if globalAlignment:
                if self.isReversed:
                    assert len(self.readSeq) - self.readPos - 1 >= 0
                    return len(self.readSeq) - self.readPos - 1
                return self.readPos
            return 0
        return self._indelLength(self.readPos, self.pPair.readPos)
    
    def getPrecedingReadDeletionLength(self, globalAlignment=False):
        if self.pPair == None:
            if globalAlignment:
                return self.refPos
            return 0
        return self._indelLength(self.refPos, self.pPair.refPos)
    
    @staticmethod
    def _indelLength(pos, pPos):
        length = abs(pPos - pos) - 1
        assert length >= 0
        return length

    @staticmethod
    def iterator(alignedRead, refSeq, readSeq):
        """Generates aligned pairs from a pysam.AlignedRead object.
        """
        readOffset = getAbsoluteReadOffset(alignedRead, refSeq, readSeq)
        pPair = None
        assert len(alignedRead.seq) <= len(readSeq)
        for readPos, refPos in alignedRead.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedRead.pos and refPos < alignedRead.aend
                if refPos >= len(refSeq): #This is masking an (apparently minor?) one 
                    #off error in the BWA sam files?
                    logger.critical("Detected an aligned reference position out of \
                    bounds! Reference length: %s, reference coordinate: %s" % \
                    (len(refSeq), refPos))
                    continue
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), 
                                 alignedRead.is_reverse, readSeq, pPair)
                if aP.getReadBase().upper() != alignedRead.query[readPos].upper():
                    logger.critical("Detected a discrepancy between the absolute read \
                    sequence and the aligned read sequence. Bases: %s %s, \
                    read-position: %s, is reversed: %s, absolute read offset: %s, \
                    length absolute read sequence %s, length aligned read sequence %s, \
                    length aligned read sequence plus soft clipping %s, read name: %s, \
                    cigar string %s" % (aP.getReadBase().upper(), 
                                        alignedRead.query[readPos].upper(), readPos, 
                                        alignedRead.is_reverse, readOffset, len(readSeq), 
                                        len(alignedRead.query), len(alignedRead.seq), 
                                        alignedRead.qname, alignedRead.cigarstring))
                assert aP.getReadBase().upper() == alignedRead.query[readPos].upper()
                pPair = aP
                yield aP

def getAbsoluteReadOffset(alignedRead, refSeq, readSeq):
    """Gets the absolute starting coordinate of the first non-clipped position in the read.
    """
    if alignedRead.cigar[0][0] == 5: #Translate the read position to the original 
        #coordinates by removing hard clipping
        readOffset = alignedRead.cigar[0][1]
    else:
        readOffset = 0
    if alignedRead.is_reverse: #SEQ is reverse complemented
        readOffset = -(len(readSeq) - 1 - readOffset)
    readOffset += alignedRead.qstart #This removes any soft-clipping
    return readOffset

def getExonerateCigarFormatString(alignedRead, sam):
    """Gets a complete exonerate like cigar-string describing the sam line
    """
    for op, length in alignedRead.cigar:
        assert op in (0, 1, 2, 4, 5)
    translation = { 0:"M", 1:"I", 2:"D" }
    cigarString = " ".join([ "%s %i" % (translation[op], length) for op, length in 
                            alignedRead.cigar if op in translation ]) 
    completeCigarString = "cigar: %s %i %i + %s %i %i + 1 %s" % (
    alignedRead.qname, 0, alignedRead.qend - alignedRead.qstart, 
    sam.getrname(alignedRead.rname), alignedRead.pos, alignedRead.aend, cigarString)
    pA = cigarReadFromString(completeCigarString) #This checks it's an okay cigar
    assert sum([ op.length for op in pA.operationList if op.type == \
                PairwiseAlignment.PAIRWISE_MATCH ]) == \
                len([ readPos for readPos, refPos in alignedRead.aligned_pairs if \
                     readPos != None and refPos != None ])
    return completeCigarString

def samToBamFile(samInputFile, bamOutputFile):
    """Converts a sam file to a bam file (sorted)
    """
    samfile = pysam.Samfile(samInputFile, "r" )
    bamfile = pysam.Samfile(bamOutputFile, "wb", template=samfile)
    for line in samfile:
        bamfile.write(line)

    samfile.close()
    bamfile.close()
    
def getFastaDictionary(fastaFile):
    """Returns a dictionary of the first words of fasta headers to their corresponding 
    fasta sequence
    """
    names = map(lambda x : x[0].split()[0], fastaRead(open(fastaFile, 'r')))
    assert len(names) == len(set(names)) #Check all the names are unique
    #Hash of names to sequences
    return dict(map(lambda x : (x[0].split()[0], x[1]), fastaRead(open(fastaFile, 'r')))) 

def getFastqDictionary(fastqFile):
    """Returns a dictionary of the first words of fastq headers to their corresponding 
    fastq sequence
    """
    names = map(lambda x : x[0].split()[0], fastqRead(open(fastqFile, 'r')))
    assert len(names) == len(set(names)) #Check all the names are unique
    #Hash of names to sequences
    return dict(map(lambda x : (x[0].split()[0], x[1]), fastqRead(open(fastqFile, 'r')))) 

def makeFastaSequenceNamesUnique(inputFastaFile, outputFastaFile):
    """Makes a fasta file with unique names
    """
    names = set()
    fileHandle = open(outputFastaFile, 'w')
    for name, seq in fastaRead(open(inputFastaFile, 'r')):
        while name in names:
            logger.critical("Got a duplicate fasta sequence name: %s" % name)
            name += "i"
        names.add(name)
        fastaWrite(fileHandle, name, seq)
    fileHandle.close()
    return outputFastaFile

def makeFastqSequenceNamesUnique(inputFastqFile, outputFastqFile):
    """Makes a fastq file with unique names
    """
    names = set()
    fileHandle = open(outputFastqFile, 'w')
    for name, seq, quals in fastqRead(open(inputFastqFile, 'r')):
        name = name.split()[0] #Get rid of any white space
        while name in names:
            logger.critical("Got a duplicate fastq sequence name: %s" % name)
            name += "i"
        names.add(name)
        fastqWrite(fileHandle, name, seq, quals)
    fileHandle.close()
    return outputFastqFile

def normaliseQualValues(inputFastqFile, outputFastqFile):
    """Makes a fastq with valid qual values
    """
    fileHandle = open(outputFastqFile, 'w')
    for name, seq, quals in fastqRead(open(inputFastqFile, 'r')):
        if quals == None:
            quals = [33] * len(seq)
        fastqWrite(fileHandle, name, seq, quals)
    fileHandle.close()
    return outputFastqFile

def samIterator(sam):
    """Creates an iterator over the aligned reads in a sam file, filtering out
    any reads that have no reference alignment.
    """
    for aR in sam:
        if aR.rname != -1:
            yield aR

def combineSamFiles(baseSamFile, extraSamFiles, outputSamFile):
    """Combines the lines from multiple sam files into one sam file
    """
    sam = pysam.Samfile(baseSamFile, "r" )
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    sam.close()
    for samFile in [ baseSamFile  ] + extraSamFiles:
        sam = pysam.Samfile(samFile, "r" )
        for line in sam:
            outputSam.write(line)
        sam.close()
    outputSam.close()
    
def paralleliseSamProcessingTargetFn(target, samFile, 
                            referenceFastaFile, outputFile, 
                            childTargetFn, followOnTargetFn, options):
    """Parallelise a computation over the alignments in a SAM file.
    """
    #Load reference sequences
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    
    tempOutputFiles = []
    childCount, totalSeqLength = 0, sys.maxint
    tempExonerateFile, tempQueryFile = None, None 
    tempExonerateFileHandle, tempQueryFileHandle = None, None
    refName = None
    
    #Read through the SAM file
    sam = pysam.Samfile(samFile, "r" )  
    
    def makeChild():
        #Add a child target to do the processing of a subset of the lines.
        if tempExonerateFile != None:
            tempExonerateFileHandle.close()
            tempQueryFileHandle.close()
            #Temporary cigar file to store the realignment
            tempOutputFiles.append(os.path.join(target.getGlobalTempDir(), 
                                               "tempOutput_%i.txt" % childCount))
            target.addChildTargetFn(childTargetFn,
                                    args=(tempExonerateFile, sam.getrname(aR.rname), 
                                          refSequences[sam.getrname(aR.rname)], 
                                          tempQueryFile, tempOutputFiles[-1], options))
    
    for aR, index in zip(samIterator(sam), xrange(sys.maxint)): 
        #Iterate on the sam lines realigning them in parallel
        if totalSeqLength > options.maxAlignmentLengthPerJob or \
        refName != sam.getrname(aR.rname):
            makeChild()
            tempExonerateFile = os.path.join(target.getGlobalTempDir(), 
                                             "tempExonerateCigar_%s.cig" % childCount)
            tempExonerateFileHandle = open(tempExonerateFile, 'w')
            tempQueryFile = os.path.join(target.getGlobalTempDir(), 
                                         "tempQueryCigar_%s.cig" % childCount)
            tempQueryFileHandle = open(tempQueryFile, 'w')
            childCount += 1
            totalSeqLength = 0
        
        tempExonerateFileHandle.write(getExonerateCigarFormatString(aR, sam) + "\n")
        fastaWrite(tempQueryFileHandle, aR.qname, aR.query)
        totalSeqLength += len(aR.query)
        refName = sam.getrname(aR.rname)
            
    makeChild()
    
    target.setFollowOnTargetFn(followOnTargetFn, args=(samFile, referenceFastaFile, \
                                                       outputFile, tempOutputFiles, options))
    #Finish up
    sam.close()
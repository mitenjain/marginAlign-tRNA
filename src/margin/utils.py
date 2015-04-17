import pysam, sys, os, random
from jobTree.src.bioio import fastaRead, fastqRead, \
cigarReadFromString,PairwiseAlignment, fastaWrite, fastqWrite, logger, absSymPath, reverseComplementChar

def pathToBaseNanoporeDir():
    """Returns path to base directory "marginAlign"
    """
    import marginAlign
    i = absSymPath(__file__)
    return os.path.split(os.path.split(os.path.split(i)[0])[0])[0]

def getFirstNonClippedPositionInRead(alignedSegment, readSeq):
    """Gets the coordinate of the first non-clipped position in the read relative to the 
    complete read sequence (including any hard clipped bases).
    If the alignment is on the reverse strand the coordinate is negative, e.g. the reverse strand coordinate of
    the 2nd position of the read sequence is -1 (0 based).
    """
    if alignedSegment.cigar[0][0] == 5: #Translate the read position to the original 
        #coordinates by removing hard clipping
        readOffset = alignedSegment.cigar[0][1]
    else:
        readOffset = 0
    if alignedSegment.is_reverse: #SEQ is reverse complemented
        readOffset = -(len(readSeq) - 1 - readOffset)
    readOffset += alignedSegment.query_alignment_start #This removes any soft-clipping
    return readOffset

def getLastNonClippedPositionInRead(alignedSegment, readSeq):
    """As getFirstNonClippedPositionInRead, but returns the last
    non-clipped position in the read, relative to the complete read sequence.
    """
    return getFirstNonClippedPositionInRead(alignedSegment, readSeq) + \
        alignedSegment.query_alignment_end - alignedSegment.query_alignment_start -1

def getExonerateCigarFormatString(alignedSegment, sam):
    """Gets a complete exonerate like cigar-string describing the sam line,
    with the cigar described with respect to the alignedSegment.query_sequence string,
    which includes softclip bases, but not hard-clipped bases.
    """
    for op, length in alignedSegment.cigar:
        assert op in (0, 1, 2, 4, 5)
    translation = { 0:"M", 1:"I", 2:"D" }
    cigarString = " ".join([ "%s %i" % (translation[op], length) for op, length in 
                            alignedSegment.cigar if op in translation ]) #Ignore soft clips
    completeCigarString = "cigar: %s %i %i + %s %i %i + 1 %s" % (
    alignedSegment.query_name, alignedSegment.qstart, alignedSegment.query_alignment_end, 
    sam.getrname(alignedSegment.reference_id), alignedSegment.reference_start, alignedSegment.reference_end, cigarString)
    ##Assertions
    pA = cigarReadFromString(completeCigarString) #This checks it's an okay cigar
    assert sum([ op.length for op in pA.operationList if op.type == \
                PairwiseAlignment.PAIRWISE_MATCH ]) == \
                len([ readPos for readPos, refPos in alignedSegment.aligned_pairs if \
                     readPos != None and refPos != None ])
    ##End assertions
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
    namesAndSequences = map(lambda x : (x[0].split()[0], x[1]), fastaRead(open(fastaFile, 'r')))
    names = map(lambda x : x[0], namesAndSequences)
    assert len(names) == len(set(names)) #Check all the names are unique
    return dict(namesAndSequences) #Hash of names to sequences

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

def samIterator(sam):
    """Creates an iterator over the aligned reads in a sam file, filtering out
    any reads that have no reference alignment.
    """
    for aR in sam:
        if aR.reference_id != -1:
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
                                    args=(tempExonerateFile, refName, 
                                          refSequences[refName], 
                                          tempQueryFile, tempOutputFiles[-1], options))
    
    for aR, index in zip(samIterator(sam), xrange(sys.maxint)): 
        #Iterate on the sam lines realigning them in parallel
        if totalSeqLength > options.maxAlignmentLengthPerJob or \
        refName != sam.getrname(aR.reference_id):
            makeChild()
            tempExonerateFile = os.path.join(target.getGlobalTempDir(), 
                                             "tempExonerateCigar_%s.cig" % childCount)
            tempExonerateFileHandle = open(tempExonerateFile, 'w')
            tempQueryFile = os.path.join(target.getGlobalTempDir(), 
                                         "tempQueryCigar_%s.fa" % childCount)
            tempQueryFileHandle = open(tempQueryFile, 'w')
            childCount += 1
            totalSeqLength = 0
        
        tempExonerateFileHandle.write(getExonerateCigarFormatString(aR, sam) + "\n")
        fastaWrite(tempQueryFileHandle, aR.query_name, aR.query_sequence) #This is the query sequence, including soft clipped bases, but excluding hard clip bases
        totalSeqLength += len(aR.query_sequence)
        refName = sam.getrname(aR.reference_id)
            
    makeChild()
    
    target.setFollowOnTargetFn(followOnTargetFn, args=(samFile, referenceFastaFile, \
                                                       outputFile, tempOutputFiles, options))
    #Finish up
    sam.close()
    
###The following code is used by the tests/plots

def getFastqDictionary(fastqFile):
    """Returns a dictionary of the first words of fastq headers to their corresponding 
    fastq sequence
    """
    namesAndSequences = map(lambda x : (x[0].split()[0], x[1]), fastqRead(open(fastqFile, 'r')))
    names = map(lambda x : x[0], namesAndSequences)
    assert len(names) == len(set(names)) #Check all the names are unique
    return dict(namesAndSequences) #Hash of names to sequences

class AlignedPair:
    """Represents an aligned pair of positions using absolute reference/read coordinates.
    
    Originally coded when I was figuring out pySam, hence is full of assertions and uses
    global coordinates.
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
        self.bases = set([ 'A', 'C', 'G', 'T' ])
        
    def isMatch(self):
        return self.getRefBase().upper() == self.getReadBase().upper() and \
            self.getRefBase().upper() in self.bases
    
    def isMismatch(self):
        return self.getRefBase().upper() != self.getReadBase().upper() and \
            self.getRefBase().upper() in self.bases and self.getReadBase().upper() in self.bases
    
    def getRefBase(self):
        return self.refSeq[self.refPos]
    
    def getReadBase(self):
        if self.isReversed:
            return reverseComplementChar(self.readSeq[self.readPos]) 
        return self.readSeq[self.readPos]
    
    def getSignedReadPos(self):
        if self.isReversed:
            return -self.readPos
        return self.readPos
    
    def getPrecedingReadInsertionLength(self, globalAlignment=False):
        #If global alignment flag is true any unaligned prefix/suffix insertion at the beginning
        #and end of the read sequence is interpreted as an insertion, rather than being ignored.
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
    def iterator(alignedSegment, refSeq, readSeq):
        """Generates aligned pairs from a pysam.AlignedSegment object.
        """
        readOffset = getFirstNonClippedPositionInRead(alignedSegment, readSeq)
        pPair = None
        assert len(alignedSegment.query_sequence) <= len(readSeq)
        for readPos, refPos in alignedSegment.aligned_pairs: #Iterate over the block
            if readPos != None and refPos != None:
                assert refPos >= alignedSegment.reference_start and refPos < alignedSegment.reference_end
                if refPos >= len(refSeq): #This is masking an (apparently minor?) one 
                    #off error in the BWA sam files?
                    logger.critical("Detected an aligned reference position out of \
                    bounds! Reference length: %s, reference coordinate: %s" % \
                    (len(refSeq), refPos))
                    continue
                aP = AlignedPair(refPos, refSeq, abs(readOffset + readPos), 
                                 alignedSegment.is_reverse, readSeq, pPair)
                if aP.getReadBase().upper() != alignedSegment.query_alignment_sequence[readPos].upper():
                    logger.critical("Detected a discrepancy between the absolute read \
                    sequence and the aligned read sequence. Bases: %s %s, \
                    read-position: %s, is reversed: %s, absolute read offset: %s, \
                    length absolute read sequence %s, length aligned read sequence %s, \
                    length aligned read sequence plus soft clipping %s, read name: %s, \
                    cigar string %s" % (aP.getReadBase().upper(), 
                                        alignedSegment.query_alignment_sequence[readPos].upper(), readPos, 
                                        alignedSegment.is_reverse, readOffset, len(readSeq), 
                                        len(alignedSegment.query_alignment_sequence), len(alignedSegment.query_sequence), 
                                        alignedSegment.query_name, alignedSegment.cigarstring))

                pPair = aP
                yield aP 
                
class ReadAlignmentStats:
    """Calculates stats of a given read alignment.
    Global alignment means the entire reference and read sequences (trailing indels).
    """
    def __init__(self, readSeq, refSeq, alignedRead, globalAlignment=False):
        self.matches = 0
        self.mismatches = 0
        self.ns = 0
        self.totalReadInsertionLength = 0
        self.totalReadInsertions = 0
        self.totalReadDeletionLength = 0
        self.totalReadDeletions = 0
        self.readSeq = readSeq
        self.refSeq = refSeq
        self.globalAlignment = globalAlignment
        
        #Now process the read alignment
        totalReadInsertionLength, totalReadDeletionLength = 0, 0
        aP = None
        for aP in AlignedPair.iterator(alignedRead, self.refSeq, self.readSeq): 
            if aP.isMatch():
                self.matches += 1
            elif aP.isMismatch():
                self.mismatches += 1
            else:
                self.ns += 1
            if aP.getPrecedingReadInsertionLength(self.globalAlignment) > 0:
                self.totalReadInsertions += 1
                totalReadInsertionLength += aP.getPrecedingReadInsertionLength(self.globalAlignment)
            if aP.getPrecedingReadDeletionLength(self.globalAlignment) > 0:
                self.totalReadDeletions += 1
                totalReadDeletionLength += aP.getPrecedingReadDeletionLength(self.globalAlignment)
        if self.globalAlignment and aP != None: #If global alignment account for any trailing indels
            assert len(self.refSeq) - aP.refPos - 1 >= 0
            if len(self.refSeq) - aP.refPos - 1 > 0:
                self.totalReadDeletions += 1
                self.totalReadDeletionLength += len(self.refSeq) - aP.refPos - 1

            if alignedRead.is_reverse:
                aP.readPos >= 0
                if aP.readPos > 0:
                    self.totalReadInsertions += 1
                    totalReadInsertionLength += aP.readPos
            else:
                assert len(self.readSeq) - aP.readPos - 1 >= 0
                if len(self.readSeq) - aP.readPos - 1 > 0:
                    self.totalReadInsertions += 1
                    totalReadInsertionLength += len(self.readSeq) - aP.readPos - 1

        assert totalReadInsertionLength <= len(self.readSeq)
        assert totalReadDeletionLength <= len(self.refSeq)
        self.totalReadInsertionLength += totalReadInsertionLength
        self.totalReadDeletionLength += totalReadDeletionLength
    
    @staticmethod
    def formatRatio(numerator, denominator):
        if denominator == 0:
            return float("nan")
        return float(numerator)/denominator
            
    def readCoverage(self):
        return self.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def referenceCoverage(self):
        return self.formatRatio(self.matches + self.mismatches, self.matches + self.mismatches + self.totalReadDeletionLength)
    
    def identity(self):
        return self.formatRatio(self.matches, self.matches + self.mismatches + self.totalReadInsertionLength)
    
    def mismatchesPerAlignedBase(self):
        return self.formatRatio(self.mismatches, self.matches + self.mismatches)
    
    def deletionsPerReadBase(self):
        return self.formatRatio(self.totalReadDeletions, self.matches + self.mismatches)
    
    def insertionsPerReadBase(self):
        return self.formatRatio(self.totalReadInsertions, self.matches + self.mismatches)
    
    def readLength(self):
        return len(self.readSeq)
    
    @staticmethod
    def getReadAlignmentStats(samFile, readFastqFile, referenceFastaFile, globalAlignment=True):
        """Gets a list of ReadAlignmentStats objects, one for each alignment in the same file
        """
        refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
        readSequences = getFastqDictionary(readFastqFile) #Hash of names to sequences
        sam = pysam.Samfile(samFile, "r")
        readsToReadCoverages = {}
        readAlignmentStats = map(lambda aR : ReadAlignmentStats(readSequences[aR.qname], \
            refSequences[sam.getrname(aR.rname)], aR, globalAlignment), samIterator(sam))
        sam.close()
        return readAlignmentStats
    
##Functions used in the creation of mutations - used in the testing of margin-caller

def mutateSequence(sequence, snpRate): #Does not preserve softmasking
    """Returns sequence with snpRate proportion of sites mutated and a list of those mutations
    """
    mutations = []
    mutatedSequence = list(sequence)
    for i in xrange(len(sequence)):
        if random.random() < snpRate:
            base = sequence[i]
            altBase = random.choice(list(set(("A", 'C', 'G', 'T')) - set(base.upper())))
            altBase = altBase if base.upper() == base else altBase.lower()
            mutations.append((i, altBase))
            mutatedSequence[i] = altBase
    return "".join(mutatedSequence), mutations

def mutateSequences(sequences, snpRate):
    """As mutateSequence, but for collection of sequences. Sequences is a dictionary 
    of sequences of names to sequences. Return value is a dictionary of names to mutated
    sequences and a list of those mutations, represented as triples of (sequenceName, position, alt).
    """
    mutatedSequences = {}; allMutations = [] #List of refSequenceName, position, altBase
    for name in sequences.keys():
        mutatedSequence, mutations = mutateSequence(sequences[name], snpRate)
        mutatedSequences[name] = mutatedSequence
        allMutations += map(lambda x : (name, x[0], x[1]), mutations) 
    return mutatedSequences, allMutations

import pysam, sys, os
from jobTree.src.bioio import reverseComplement, fastaRead, fastqRead, \
cigarReadFromString,PairwiseAlignment, fastaWrite, fastqWrite, logger, absSymPath

def pathToBaseNanoporeDir():
    """Returns path to base directory "marginAlign"
    """
    import marginAlign
    i = absSymPath(__file__)
    return os.path.split(os.path.split(i)[0])[0]

def getFirstNonClippedPositionInRead(alignedRead, readSeq):
    """Gets the coordinate of the first non-clipped position in the read relative to the 
    complete read sequence (including any hard clipped bases).
    If the alignment is on the reverse strand the coordinate is negative, e.g. the reverse strand coordinate of
    the 2nd position of the read sequence is -1 (0 based).
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

def getLastNonClippedPositionInRead(alignedRead, readSeq):
    """As getFirstNonClippedPositionInRead, but returns the last
    non-clipped position in the read, relative to the complete read sequence.
    """
    return getFirstNonClippedPositionInRead(alignedRead, readSeq) + alignedRead.qend- alignedRead.qstart -1

def getExonerateCigarFormatString(alignedRead, sam):
    """Gets a complete exonerate like cigar-string describing the sam line,
    with the cigar described with respect to the alignedRead.query_sequence string,
    which includes softclip bases, but not hard-clipped bases.
    """
    for op, length in alignedRead.cigar:
        assert op in (0, 1, 2, 4, 5)
    translation = { 0:"M", 1:"I", 2:"D" }
    cigarString = " ".join([ "%s %i" % (translation[op], length) for op, length in 
                            alignedRead.cigar if op in translation ]) #Ignore soft clips
    completeCigarString = "cigar: %s %i %i + %s %i %i + 1 %s" % (
    alignedRead.qname, alignedRead.qstart, alignedRead.qend, 
    sam.getrname(alignedRead.rname), alignedRead.pos, alignedRead.aend, cigarString)
    ##Assertions
    pA = cigarReadFromString(completeCigarString) #This checks it's an okay cigar
    assert sum([ op.length for op in pA.operationList if op.type == \
                PairwiseAlignment.PAIRWISE_MATCH ]) == \
                len([ readPos for readPos, refPos in alignedRead.aligned_pairs if \
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
    names = map(lambda x : x[0].split()[0], fastaRead(open(fastaFile, 'r')))
    assert len(names) == len(set(names)) #Check all the names in the sequence file are unique
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
                                    args=(tempExonerateFile, refName, 
                                          refSequences[refName], 
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
                                         "tempQueryCigar_%s.fa" % childCount)
            tempQueryFileHandle = open(tempQueryFile, 'w')
            childCount += 1
            totalSeqLength = 0
        
        tempExonerateFileHandle.write(getExonerateCigarFormatString(aR, sam) + "\n")
        fastaWrite(tempQueryFileHandle, aR.qname, aR.seq) #This is the query sequence, including soft clipped bases, but excluding hard clip bases
        totalSeqLength += len(aR.seq)
        refName = sam.getrname(aR.rname)
            
    makeChild()
    
    target.setFollowOnTargetFn(followOnTargetFn, args=(samFile, referenceFastaFile, \
                                                       outputFile, tempOutputFiles, options))
    #Finish up
    sam.close()
    
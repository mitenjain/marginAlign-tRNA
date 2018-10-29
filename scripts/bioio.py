#!/usr/bin/env python

#Copyright (C) 2006-2012 by Benedict Paten (benedictpaten@gmail.com)
#
#Released under the MIT license, see LICENSE.txt

import sys
import os
import re
import random
import math
import shutil
from argparse import ArgumentParser
from optparse import OptionParser, OptionContainer, OptionGroup
import subprocess
import array

DEFAULT_DISTANCE = 0.001


#########################################################
#########################################################
#########################################################
#system wrapper command
#########################################################
#########################################################
#########################################################

def system(command):
    logger.debug("Running the command: %s" % command)
    sts = subprocess.call(command, shell=True, bufsize=-1, stdout=sys.stdout, stderr=sys.stderr)
    if sts != 0:
        raise RuntimeError("Command: %s exited with non-zero status %i" % (command, sts))
    return sts


#########################################################
#########################################################
#########################################################
#fasta/fastq functions
#########################################################
#########################################################
#########################################################

def fastaNormaliseHeader(fastaHeader):
    """Removes white space which is treated weirdly by many programs.
    """
    i = fastaHeader.split()
    if len(i) > 0:
        return i[0]
    return ""

def fastaDecodeHeader(fastaHeader):
    """Decodes the fasta header
    """
    return fastaHeader.split("|")

def fastaEncodeHeader(attributes):
    """Decodes the fasta header
    """
    for i in attributes:
        assert len(str(i).split()) == 1
    return "|".join([ str(i) for i in attributes ])

def _getFileHandle(fileHandleOrFile, mode="r"):
    if isinstance(fileHandleOrFile, "".__class__):
        return open(fileHandleOrFile, mode)
    else:
        return fileHandleOrFile

def fastaRead(fileHandleOrFile):
    """iteratively a sequence for each '>' it encounters, ignores '#' lines
    """
    fileHandle = _getFileHandle(fileHandleOrFile)
    line = fileHandle.readline()
    while line != '':
        if line[0] == '>':
            name = line[1:-1]
            line = fileHandle.readline()
            seq = array.array('c')
            while line != '' and line[0] != '>':
                if line[0] != '#':
                    seq.extend([ i for i in line[:-1] if not i.isspace() ]) #The white-space check is to remove any annoying trailing characters.
                line = fileHandle.readline()
            for i in seq:
                #For safety and sanity I only allows roman alphabet characters in fasta sequences.
                if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'):
                    raise RuntimeError("Invalid FASTA character, ASCII code = \'%d\', found in input sequence %s" % (ord(i), name))
            yield name, seq.tostring()
        else:
            line = fileHandle.readline()
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()

def fastaWrite(fileHandleOrFile, name, seq, mode="w"):
    """Writes out fasta file
    """
    fileHandle = _getFileHandle(fileHandleOrFile, mode)
    assert seq.__class__ == "".__class__
    for i in seq:
        assert (i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-' #For safety and sanity I only allows roman alphabet characters in fasta sequences.
    fileHandle.write(">%s\n" % name)
    chunkSize = 100
    for i in xrange(0, len(seq), chunkSize):
        fileHandle.write("%s\n" % seq[i:i+chunkSize])
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()
        
def fastqRead(fileHandleOrFile):
    """Reads a fastq file iteratively
    """
    fileHandle = _getFileHandle(fileHandleOrFile)
    line = fileHandle.readline()
    while line != '':
        if line[0] == '@':
            name = line[1:-1]
            seq = fileHandle.readline()[:-1]
            plus = fileHandle.readline()
            if plus[0] != '+':
                raise RuntimeError("Got unexpected line: %s" % plus)
            qualValues = [ ord(i) for i in fileHandle.readline()[:-1] ]
            if len(seq) != len(qualValues):
                logger.critical("Got a mismatch between the number of sequence characters (%s) and number of qual values (%s) for sequence: %s, ignoring returning None" % (len(seq), len(qualValues), name))
                qualValues = None
            else:
                for i in qualValues:
                    if i < 33 or i > 126:
                        raise RuntimeError("Got a qual value out of range %s (range is 33 to 126)" % i)
            for i in seq:
                #For safety and sanity I only allows roman alphabet characters in fasta sequences.
                if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'):
                    raise RuntimeError("Invalid FASTQ character, ASCII code = \'%d\', found in input sequence %s" % (ord(i), name))
            yield name, seq, qualValues
        line = fileHandle.readline()
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()

def fastqWrite(fileHandleOrFile, name, seq, qualValues, mode="w"):
    """Writes out fastq file. If qualValues is None or '*' then prints a '*' instead.
    """
    fileHandle = _getFileHandle(fileHandleOrFile, mode)
    assert seq.__class__ == "".__class__
    for i in seq:
        if not ((i >= 'A' and i <= 'Z') or (i >= 'a' and i <= 'z') or i == '-'): #For safety and sanity I only allows roman alphabet characters in fasta sequences.
            raise RuntimeError("Invalid FASTQ character, ASCII code = \'%d\', char = '%s' found in input sequence %s" % (ord(i), i, name))
    if qualValues != None and qualValues != '*':
        if len(seq) != len(qualValues):
            raise RuntimeError("Got a mismatch between the number of sequence characters (%s) and number of qual values (%s) for sequence: %s " % (len(seq), len(qualValues), name))
        for i in qualValues:
            if i < 33 or i > 126:
                raise RuntimeError("Got a qual value out of range %s (range is 33 to 126)" % i)
        fileHandle.write("@%s\n%s\n+\n%s\n" % (name, seq, "".join([ chr(i) for i in qualValues ])))
    else:
        fileHandle.write("@%s\n%s\n+\n*\n" % (name, seq))
    if isinstance(fileHandleOrFile, "".__class__):
        fileHandle.close()

def fastaReadHeaders(fasta):
    """Returns a list of fasta header lines, excluding
    """
    headers = []
    fileHandle = open(fasta, 'r')
    line = fileHandle.readline()
    while line != '':
        assert line[-1] == '\n'
        if line[0] == '>':
            headers.append(line[1:-1])
        line = fileHandle.readline()
    fileHandle.close()
    return headers

def getRandomSequence(length=500):
    """Generates a random name and sequence.
    """
    fastaHeader = ""
    for i in xrange(int(random.random()*100)):
        fastaHeader = fastaHeader + random.choice([ 'A', 'C', '0', '9', ' ', '\t' ])
    return (fastaHeader, \
            "".join([ random.choice([ 'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'A', 'C', 'T', 'G', 'N' ]) for i in xrange((int)(random.random() * length))]))

def _expLength(i=0, prob=0.95):
    if random.random() >= prob:
        return _expLength(i+1)
    return i

def mutateSequence(seq, distance):
    """Mutates the DNA sequence for use in testing.
    """
    subProb=distance
    inProb=0.05*distance
    deProb=0.05*distance
    contProb=0.9
    l = []
    bases = [ 'A', 'C', 'T', 'G' ]
    i=0
    while i < len(seq):
        if random.random() < subProb:
            l.append(random.choice(bases))
        else:
            l.append(seq[i])
        if random.random() < inProb:
            l += getRandomSequence(_expLength(0, contProb))[1]
        if random.random() < deProb:
            i += int(_expLength(0, contProb))
        i += 1
    return "".join(l)

def reverseComplement(seq):
    seq = list(seq)
    seq.reverse()
    dNA = { 'A':'T', 'T':'A', 'C':'G', 'G':'C', 'a':'t', 't':'a', 'c':'g', 'g':'c' }
    def fn(i):
        if i in dNA:
            return dNA[i]
        return i
    return "".join([ fn(i) for i in seq ])

def main():
    pass

def _test():
    import doctest
    return doctest.testmod()

if __name__ == '__main__':
    _test()
    main()

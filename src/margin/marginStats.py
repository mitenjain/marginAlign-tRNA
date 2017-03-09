import os
import sys
from optparse import OptionParser
from jobTree.src.bioio import logger, setLoggingFromOptions, addLoggingOptions
from margin.utils import ReadAlignmentStats
import numpy

"""Simple script to calculate alignment statistics from a SAM file.
"""
    
def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: samFile, readFastqFile, referenceFastaFile [options]", 
                          version="%prog 0.1")
    
    #Options
    parser.add_option("--readIdentity", dest="readIdentity", 
                      help="Print readIdentity of alignments", 
                      default=False, action="store_true")
    
    parser.add_option("--alignmentIdentity", dest="alignmentIdentity", 
                      help="Print alignmentIdentity", 
                      default=False, action="store_true")
    
    parser.add_option("--readCoverage", dest="readCoverage", 
                      help="Print read coverage of alignments", 
                      default=False, action="store_true")
    
    parser.add_option("--mismatchesPerAlignedBase", dest="mismatchesPerAlignedBase", 
                      help="Print mismatches per aligned base", 
                      default=False, action="store_true")
    
    parser.add_option("--deletionsPerReadBase", dest="deletionsPerReadBase", 
                      help="Print deletions per base of alignments", 
                      default=False, action="store_true")
    
    parser.add_option("--insertionsPerReadBase", dest="insertionsPerReadBase", 
                      help="Print insertions per base of alignments", 
                      default=False, action="store_true")
    
    parser.add_option("--readLength", dest="readLength", 
                      help="Print read lengths of aligned reads", 
                      default=False, action="store_true")

    parser.add_option("--localAlignment", dest="localAlignment", 
                      help="Ignore unaligned prefix and suffix of each read in making calculation", 
                      default=False, action="store_true")
    
    parser.add_option("--printValuePerReadAlignment", dest="printValuePerReadAlignment", 
                      help="Prints the value of statistics for each read alignment", 
                      default=False, action="store_true")
    
    parser.add_option("--noStats", dest="noStats", 
                      help="Do not print stats (avg, median, min, max, mode) of desired statistic", 
                      default=False, action="store_true")
    
    addLoggingOptions(parser)
    
    #Parse the options/arguments
    options, args = parser.parse_args()
    
    #Setup logging
    setLoggingFromOptions(options)
    
    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    #Exit if the arguments are not what we expect
    if len(args) != 3:
        raise RuntimeError("Expected three arguments, got: %s" % " ".join(args))
    
    #Now do the stats calculation
    samFile, readFastqFile, referenceFastaFile = args
    
    readAlignmentStats = ReadAlignmentStats.getReadAlignmentStats(samFile, readFastqFile, 
                                             referenceFastaFile, globalAlignment=not options.localAlignment)
    
    def report(values, statisticName):
        if not options.noStats:
            print "Average" + statisticName, numpy.average(values)
            print "Median" + statisticName, numpy.median(values)
            print "Min" + statisticName, min(values)
            print "Max" + statisticName, max(values)
        if options.printValuePerReadAlignment:
            print "Values" + statisticName, "\t".join(map(str, values))
    
    if options.readIdentity:
        report(map(lambda rAS : rAS.readIdentity(), readAlignmentStats), "ReadIdentity")
    
    if options.alignmentIdentity:
        report(map(lambda rAS : rAS.alignmentIdentity(), readAlignmentStats), "AlignmentIdentity")
    
    if options.readCoverage:
        report(map(lambda rAS : rAS.readCoverage(), readAlignmentStats), "ReadCoverage")
    
    if options.mismatchesPerAlignedBase:
        report(map(lambda rAS : rAS.mismatchesPerAlignedBase(), readAlignmentStats), "MismatchesPerAlignedBase")
    
    if options.deletionsPerReadBase:
        report(map(lambda rAS : rAS.deletionsPerReadBase(), readAlignmentStats), "DeletionsPerReadBase")
    
    if options.insertionsPerReadBase:
        report(map(lambda rAS : rAS.insertionsPerReadBase(), readAlignmentStats), "InsertionsPerReadBase")

    if options.readLength:
        report(map(lambda rAS : rAS.readLength(), readAlignmentStats), "ReadLength")

if __name__ == '__main__':
    main()

import os
import sys
from optparse import OptionParser
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from jobTree.scriptTree.stack import Stack
from margin.mappers.last import LastChain, LastRealign, LastRealignEm
from margin.mappers.bwa import BwaChain, BwaRealign, BwaRealignEm
    
def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastqFile referenceFastaFile outputSamFile [options]", version="%prog 0.1")
    
    parser.add_option("--bamOutput", default=False, help="Output BAM instead of SAM", action="store_true")
    parser.add_option("--em", dest="em", help="Run expectation maximisation (EM) and output the trained model file in the given file")
    parser.add_option("--useModel", dest="model", help="Use hmm model instead of stock model")
    parser.add_option("--realign", dest="realign", help="Instead of taking the files as fastq files, assumes input files are existing SAM/BAM files")
    
    Stack.addJobTreeOptions(parser)
    options, args = parser.parse_args()
    setLoggingFromOptions(options)
    
    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    if len(args) != 3:
        raise RuntimeError("Expected three arguments, got: %s" % " ".join(args))
    
    #Pick the mapper
    mapper = LastRealign
    
    #This line invokes jobTree  
    i = Stack(mapper(readFastqFile=args[0], referenceFastaFile=args[1], 
                     outputSamFile=args[2], emptyHmmFile=None))
              
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from margin.marginAlign import *
    main()

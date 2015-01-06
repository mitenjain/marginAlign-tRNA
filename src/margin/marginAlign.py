import os
import sys
from optparse import OptionParser
from jobTree.src.bioio import getLogLevelString, isNewer, logger, setLoggingFromOptions
from jobTree.scriptTree.stack import Stack
from margin.mappers.last import LastChain, LastRealign, LastRealignEm
from margin.mappers.bwa import BwaChain, BwaRealign, BwaRealignEm
from margin.utils import pathToBaseNanoporeDir
    
def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastqFile referenceFastaFile outputSamFile [options]", 
                          version="%prog 0.1")
    
    #Options
    parser.add_option("--useModel", dest="inputModel", help="Use hmm model instead of stock model", 
                      default=os.path.join(pathToBaseNanoporeDir(), "margin", "mappers", "last_hmm_20.txt"))
    parser.add_option("--em", dest="outputModel", 
                      help="Run expectation maximisation (EM) and output the trained model file in the given file",
                      default=None)
    ##Most people would not want to use the following, but I put them here for debug purposes
    parser.add_option("--bwa", dest="bwa", help="Use BWA instead of LAST", 
                      default=False, action="store_true")
    parser.add_option("--noRealign", dest="noRealign", help="Don't run any realignment step", 
                      default=False, action="store_true")
    
    #Add the jobTree options
    Stack.addJobTreeOptions(parser)
    
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
    
    #Set the mapper
    if options.noRealign:
        mapper = BwaChain if options.bwa else LastChain
    elif options.outputModel != None:
        mapper = BwaRealignEm if options.bwa else LastRealignEm
    else:
        mapper = BwaRealign if options.bwa else LastRealign
    
    #This line invokes jobTree  
    i = Stack(mapper(readFastqFile=args[0], referenceFastaFile=args[1], outputSamFile=args[2], 
                     inputHmmFile=options.inputModel,
                     outputHmmFile=options.outputModel)).startJobTree(options) 
        
    #The return value of the jobtree script is the number of failed jobs. If we have any then
    #report this.       
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from margin.marginAlign import *
    main()

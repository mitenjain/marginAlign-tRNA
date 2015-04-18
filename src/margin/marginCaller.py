import os
import sys
from optparse import OptionParser
from jobTree.src.bioio import logger, setLoggingFromOptions
from jobTree.scriptTree.stack import Stack
from jobTree.scriptTree.target import Target
from margin.utils import pathToBaseNanoporeDir
from margin.marginCallerLib import marginCallerTargetFn
    
def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputSamFile referenceFastaFile outputVcfFile [options]", 
                          version="%prog 0.1")
    
    #Options
    parser.add_option("--noMargin", dest="noMargin", help="Do not marginalise over the read \
    alignments, rather use the input alignment to call the variants (this will be faster)", 
                      default=False, action="store_true")
    parser.add_option("--alignmentModel", default=os.path.join(pathToBaseNanoporeDir(), 
                                                          "src", "margin", "mappers", "last_hmm_20.txt"), 
                     help="The model to use in realigning the reads to the reference.")
    parser.add_option("--errorModel", default=os.path.join(pathToBaseNanoporeDir(), 
                                                          "src", "margin", "mappers", "last_hmm_20.txt"), 
                     help="The model to use in calculating the difference between the predicted true reference and the reads.")
    parser.add_option("--maxAlignmentLengthPerJob", default=7000000, 
                     help="Maximum total alignment length of alignments to include in one posterior prob calculation job.", 
                     type=int)
    parser.add_option("--threshold", default=0.3, 
                     help="The posterior probability threshold for a non-reference base above which to report a variant.", 
                     type=float)
    
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
    
    print options.errorModel
    print options.threshold
    
    #This line invokes jobTree  
    i = Stack(Target.makeTargetFn(fn=marginCallerTargetFn, args=(args[0], args[1], args[2], options))).startJobTree(options) 
        
    #The return value of the jobtree script is the number of failed jobs. If we have any then
    #report this.       
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from margin.marginCaller import *
    main()

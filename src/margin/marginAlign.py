#!/usr/bin/env python
import os
import sys
from optparse import OptionParser
from jobTree.src.bioio import logger, setLoggingFromOptions
from jobTree.scriptTree.stack import Stack
from margin.mappers.last import Last, LastChain, LastRealign
from margin.mappers.bwa import Bwa, BwaChain, BwaRealign
from margin.mappers.minimap2 import Minimap2, Minimap2Chain, Minimap2Realign
from margin.utils import pathToBaseNanoporeDir
import cPecan.cPecanEm
from cPecan.cPecanEm import addExpectationMaximisationOptions
    
def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastqFile referenceFastaFile outputSamFile [options]", 
                          version="%prog 0.1")
    
    #Options
    parser.add_option("--em", dest="em", 
                      help="Run expectation maximisation (EM)",
                      default=False, action="store_true")
    ##Most people would not want to use the following, but I put them here for debug purposes
    parser.add_option("--bwa", dest="bwa", help="Use BWA instead of LAST", 
                      default=False, action="store_true")
    parser.add_option("--noRealign", dest="noRealign", help="Don't run any realignment step", 
                      default=False, action="store_true")
    parser.add_option("--noChain", dest="noChain", help="Don't run any chaining step", 
                      default=False, action="store_true")
    parser.add_option("--gapGamma", dest="gapGamma", help="Set the gap gamma for the AMAP function", 
                      default=0.5, type=float)
    parser.add_option("--matchGamma", dest="matchGamma", help="Set the match gamma for the AMAP function", 
                      default=0.0, type=float)
    parser.add_option("--minimap2", dest="minimap2", help="Use minimap2 (map-ont) instead of LAST (primary only)", 
                      default=False, action="store_true")
    
    #Add the cPecan expectation maximisation options
    options = cPecan.cPecanEm.Options()
    options.inputModel = os.path.join(pathToBaseNanoporeDir(), "src", "margin", "mappers", "last_hmm_20.txt")
    options.modelType="fiveStateAsymmetric" #"threeStateAsymmetric"
    options.optionsToRealign="--diagonalExpansion=10 --splitMatrixBiggerThanThis=300" 
    options.randomStart = True
    options.trials = 3
    options.outputTrialHmms = True
    options.iterations = 100
    options.maxAlignmentLengthPerJob=700000
    options.maxAlignmentLengthToSample = 50000000
    #options.outputXMLModelFile = outputModel + ".xml"
    #options.updateTheBand = True
    #options.useDefaultModelAsStart = True
    #options.setJukesCantorStartingEmissions=0.3
    options.trainEmissions=True
    #options.tieEmissions = True
    addExpectationMaximisationOptions(parser, options)
    
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
        if options.noChain: # i.e. --noChain --noRealign
            if options.bwa:
                mapper = Bwa
            elif options.minimap2:
                mapper = Minimap2
            else:
                mapper = Last
        else: # i.e. --noRealign
            if options.bwa:
                mapper = BwaChain
            elif options.minimap2:
                mapper = Minimap2Chain
            else:
                mapper = LastChain
    else:
        if options.bwa:
	    mapper = BwaRealign
        elif options.minimap2:
	    mapper = Minimap2Realign
        else:
	    mapper = LastRealign
    
    #This line invokes jobTree  
    i = Stack(mapper(readFastqFile=args[0], referenceFastaFile=args[1], outputSamFile=args[2], 
                     options=options)).startJobTree(options) 
        
    #The return value of the jobtree script is the number of failed jobs. If we have any then
    #report this.       
    if i != 0:
        raise RuntimeError("Got failed jobs")

if __name__ == '__main__':
    from margin.marginAlign import *
    main()

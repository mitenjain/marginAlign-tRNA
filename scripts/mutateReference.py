import os, sys, argparse, time
from optparse import OptionParser
from margin.utils import mutateSequences, getFastaDictionary
from sonLib.bioio import fastaWrite

"""Generates a mutated reference genome and a list of the introduced mutations.
"""

def main():
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastaFile outputFastaFile outputMutationsFile [options]", 
                          version="%prog 0.1")

    parser.add_option("--snpRate", dest="snpRate", 
                      help="The probability of introducing a random different base at each position", 
                      default=0.2, type=float)

    #Parse the options/arguments
    options, args = parser.parse_args()
    
    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    
    #Exit if the arguments are not what we expect
    if len(args) != 3:
        raise RuntimeError("Expected three arguments, got: %s" % " ".join(args))
 
    #This call gets the mutated sequences and a list of mutations
    mutatedSequences, allMutations = mutateSequences(getFastaDictionary(args[0]), options.snpRate)
 
    #Write out the mutated sequences into the given file
    fH = open(args[1], 'w')
    for name in mutatedSequences:
        fastaWrite(fH, name, mutatedSequences[name])
    fH.close()
    
    #Write out mutations
    fH = open(args[2], 'w')
    for mutation in allMutations:
        fH.write("\t".join(map(str, mutation)) + "\n")
    fH.close()
    
if __name__ == '__main__':
    main()

import os, sys, argparse, time
from optparse import OptionParser
from margin.utils import makeFastaSequenceNamesUnique, makeFastqSequenceNamesUnique

"""Ensures that the first word of each fasta header is unique within the input file, 
outputting the result to a given output file.
"""

def main(myCommandLine=None):
    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastqFile outputFastqFile", 
                          version="%prog 0.1")

    #Parse the options/arguments
    options, args = parser.parse_args()

    #Print help message if no input
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    #Exit if the arguments are not what we expect
    if len(args) != 2:
        raise RuntimeError("Expected two arguments, got: %s" % " ".join(args))
 
    makeFastqSequenceNamesUnique(args[0], args[1])

if __name__ == '__main__':
    main()
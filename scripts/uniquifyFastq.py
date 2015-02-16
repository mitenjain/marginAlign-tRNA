import os, sys, argparse, time
from optparse import OptionParser
from margin.utils import makeFastaSequenceNamesUnique, makeFastqSequenceNamesUnique

##################################################################################
# Main
# Here is the main program
##################################################################################

def main(myCommandLine=None):
    # starting time
    t0 = time.time()

    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputFastqFile outputFastqFile", 
                          version="%prog 0.1")

    #Parse the options/arguments
    options, args = parser.parse_args()

    #Exit if the arguments are not what we expect
    if len(args) != 2:
        raise RuntimeError("Expected two arguments, got: %s" % " ".join(args))
 
    makeFastqSequenceNamesUnique(args[0], args[1])

    print >> sys.stderr, "\n", "Total time for the program %.3f" % (time.time()-t0), "s"

if (__name__ == "__main__"):
    main()
    raise SystemExit

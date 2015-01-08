import unittest
import os, sys, pysam

"""Basic system level tests that the marginAlign and marginVar scripts
work.
"""

def validateSam(samFile):
    """Checks if a sam file is valid
    """
    # Specify a flag to check if whether samfile had a header
    fileFlag = "None"
    
    # Check is samfile exists
    if not os.path.isfile(samFile) is True:
        print "No samfile found"
        sys.exit()

    # Try opening samfile
    samfile = pysam.Samfile(samFile, "r" )
    # Check if samfile contains a proper header, and set fileFlag accordingly
    if samfile.header:
        if "HD" and "SQ" in samfile.header:
            fileFlag = "success"

    samfile.close()

    if fileFlag != "None":
        print "Samfile validation successful"
    else:
        print "Samfile validation failed"


class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testMarginAlignNoEm(self):

        print "Running marginAlign"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --jobTree testJobTree")

        # Validate samfile
        validateSam("./tests/test.sam")

        # Clean up
        os.system("rm -rf tests/test.sam testJobTree")


    def testMarginAlignEm(self):

        print "Running marginAlign and writing output.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --em tests/output.hmm --jobTree testJobTree")

        # Validate samfile
        validateSam("./tests/test.sam")

        # Check if output.hmm files were written
        if os.path.isfile("./tests/output.hmm"):
            print "Output hmm file successful"
        else:
            print "No output hmm file"

        # Clean up
        os.system("rm -rf tests/test.sam testJobTree tests/output.hmm*")


    def testMarginAlignReadModel(self):

        # Check if input.hmm exists
        if os.path.isfile("./tests/input.hmm"):
            print "Input hmm file found"
        else:
            print "No input hmm file"
            sys.exit()

        print "Running marginAlign and reading model from input.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --useModel tests/input.hmm --jobTree testJobTree")

        # Validate samfile
        validateSam("./tests/test.sam")

        # Clean up
        os.system("rm -rf tests/test.sam testJobTree")


if __name__ == '__main__':
    unittest.main()

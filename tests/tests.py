import unittest
import os, sys, numpy, pysam

"""Basic system level tests that the marginAlign and marginVar scripts
work.
"""

class TestCase(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def validateSam(self, samFile):
        """Checks if a sam file is valid
        """
        # Specify a flag to check if whether samfile had a header
        fileFlag = "None"
        # Check if samfile exists
        self.assertTrue(os.path.isfile(samFile))
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

    def testMarginAlignNoEm(self):
        print "Running marginAlign with no EM"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")
        # Clean up
        print "Cleaning files"
        os.system("rm -rf tests/test.sam testJobTree")

    def testMarginAlignEm(self):
        print "Running marginAlign and writing output.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --em --outputModel tests/output.hmm --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")
        # Check if output.hmm files were written
        self.assertTrue(os.path.isfile("./tests/output.hmm"))
        # Clean up
        print "Cleaning files"
        os.system("rm -rf tests/test.sam testJobTree tests/output.hmm*")

    def testMarginAlignReadModel(self):
        # Check if input.hmm exists
        self.assertTrue(os.path.isfile("./tests/input.hmm"))
        print "Running marginAlign and reading model from input.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --inputModel tests/input.hmm --jobTree testJobTree")
       # Validate samfile
        self.validateSam("./tests/test.sam")
        # Clean up
        print "Cleaning files"
        os.system("rm -rf tests/test.sam testJobTree")

    def testMarginCaller(self):
        # Validate input.sam
        self.validateSam("./tests/input.sam")
        print "Running marginCaller"
        os.system("./marginCaller tests/input.sam tests/reference.fa tests/output.vcf --jobTree testJobTree")
        self.assertTrue(os.path.isfile("./tests/output.vcf"))
        # Clean up
        print "Cleaning files"
      # os.system("rm -rf tests/output.vcf testJobTree")

    def testMarginCallerNoMarginalize(self):
        # Validate input.sam
        print "Running marginCaller without marginalize option"
        self.validateSam("./tests/input.sam")
        os.system("./marginCaller tests/input.sam tests/reference.fa tests/output.vcf --noMargin --jobTree testJobTree")
        # Clean up
        print "Cleaning files"
        os.system("rm -rf tests/output.vcf testJobTree")


if __name__ == '__main__':
    unittest.main()

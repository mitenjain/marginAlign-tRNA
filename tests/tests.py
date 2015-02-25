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
        # Clean up
        print "Cleaning files"
        os.system("rm -rf tests/test.sam testJobTree tests/output.hmm* tests/output.vcf")

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

    def testMarginAlignNoEm1Ref(self):
        print "\nRunning marginAlign for 1 reference, and a pre-trained model (no EM)"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")

    def testMarginAlignEm1Ref(self):
        print "\nRunning marginAlign for 1 reference, with EM and writing output.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --em --outputModel tests/output.hmm --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")
        # Check if output.hmm files were written
        self.assertTrue(os.path.isfile("./tests/output.hmm"))

    def testMarginAlignReadModel1Ref(self):
        # Check if input.hmm exists
        self.assertTrue(os.path.isfile("./tests/input.hmm"))
        print "\nRunning marginAlign using for 1 reference, and reading model from user-provided model input.hmm"
        os.system("./marginAlign tests/reads.fq tests/reference.fa tests/test.sam --inputModel tests/input.hmm --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")

    def testMarginCaller1Ref(self):
        print "\nRunning marginCaller using for 1 reference, and marginalization"
        # Validate input.sam
        self.validateSam("./tests/input1Ref.sam")
        os.system("./marginCaller tests/input1Ref.sam tests/reference.fa tests/output.vcf --jobTree testJobTree")
        self.assertTrue(os.path.isfile("./tests/output.vcf"))

    def testMarginCallerNoMarginalize1Ref(self):
        print "\nRunning marginCaller for 1 reference, with no marginalize option"
        # Validate input.sam
        self.validateSam("./tests/input1Ref.sam")
        os.system("./marginCaller tests/input1Ref.sam tests/reference.fa tests/output.vcf --noMargin --jobTree testJobTree")

    def testMarginAlignNoEm2Refs(self):
        print "\nRunning marginAlign for 2 references, and a pre-trained model (no EM)"
        os.system("./marginAlign tests/reads.fq tests/references.fa tests/test.sam --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")

    def testMarginAlignEm2Refs(self):
        print "\nRunning marginAlign for 2 references, with EM and writing output.hmm"
        os.system("./marginAlign tests/reads.fq tests/references.fa tests/test.sam --em --outputModel tests/output.hmm --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")
        # Check if output.hmm files were written
        self.assertTrue(os.path.isfile("./tests/output.hmm"))

    def testMarginAlignReadModel2Refs(self):
        # Check if input.hmm exists
        self.assertTrue(os.path.isfile("./tests/input.hmm"))
        print "\nRunning marginAlign using for 2 references, and reading model from user-provided model input.hmm"
        os.system("./marginAlign tests/reads.fq tests/references.fa tests/test.sam --inputModel tests/input.hmm --jobTree testJobTree")
        # Validate samfile
        self.validateSam("./tests/test.sam")

    def testMarginCaller2Refs(self):
        print "\nRunning marginCaller using for 2 references, and marginalization"
        # Validate input.sam
        self.validateSam("./tests/input2Refs.sam")
        os.system("./marginCaller tests/input2Refs.sam tests/references.fa tests/output.vcf --jobTree testJobTree")
        self.assertTrue(os.path.isfile("./tests/output.vcf"))

    def testMarginCallerNoMarginalize2Refs(self):
        print "\nRunning marginCaller for 2 references, with no marginalize option"
        # Validate input.sam
        self.validateSam("./tests/input2Refs.sam")
        os.system("./marginCaller tests/input2Refs.sam tests/references.fa tests/output.vcf --noMargin --jobTree testJobTree")



if __name__ == '__main__':
    unittest.main()

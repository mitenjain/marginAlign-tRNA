import unittest
import time
import os, sys, numpy, pysam
from margin.utils import calculateIdentity, pathToBaseNanoporeDir
from cPecan.cPecanEm import Hmm
from sonLib.bioio import system, parseSuiteTestOptions, logger


"""Basic system level tests that the marginAlign and marginCaller scripts
work.
"""

class TestCase(unittest.TestCase):
    def setUp(self):
        self.marginAlign = self.getFile("marginAlign") #Path to margin align binary
        self.readFastqFile1 = self.getFile("tests/reads.fq")
        self.referenceFastaFile1 = self.getFile("tests/references.fa")
        self.outputSamFile = self.getFile("tests/test.sam")
        self.inputHmmFile = self.getFile("tests/input.hmm")
        self.outputHmmFile = self.getFile("tests/output.hmm")
        self.jobTree = self.getFile("tests/testJobTree")
        unittest.TestCase.setUp(self)

    def getFile(self, file):
        return os.path.join(pathToBaseNanoporeDir(), file)

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        # Clean up
        system("rm -rf %s %s %s" % (self.outputSamFile, self.outputHmmFile, self.jobTree))

    def validateSam(self, samFile, readFastqFile, referenceFastaFile):
        """Checks if a sam file is valid.
        """
        # Check if samfile exists
        self.assertTrue(os.path.isfile(samFile))
        #The call calculate identity will run a lot of internal consistency checks
        #as it calculates the alignment identity.
        return calculateIdentity(samFile, readFastqFile, referenceFastaFile, globalAlignment=True)
    
    def checkHmm(self, hmmFile):
        Hmm.loadHmm(hmmFile) #This performs a bunch of internal consistency checks 
       
    def runMarginAlign(self, readFastqFile, referenceFastaFile, args=""):
        startTime = time.time()
        system("\t".join([ self.marginAlign, readFastqFile,
                         referenceFastaFile, self.outputSamFile, "--jobTree=%s" % self.jobTree, args ]))
        runTime = time.time() - startTime
        identity = self.validateSam(self.outputSamFile, readFastqFile, referenceFastaFile)
        logger.info("Ran marginAlign with args: %s for with reference: %s and reads: %s. Got identity: %s Took: %s seconds" % \
                    (args, readFastqFile, referenceFastaFile, identity, runTime))

    def testMarginAlignDefaults(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1)
    
    def testMarginAlignNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noChain")
    
    def testMarginAlignNoRealign(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noRealign")
    
    def testMarginAlignNoRealignNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noRealign --noChain")

    def testMarginAlignEm(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--em --outputModel %s" % self.outputHmmFile)
        self.checkHmm(self.outputHmmFile)
    
    def testMarginAlignEmNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noChain --em --outputModel %s" % self.outputHmmFile)
        self.checkHmm(self.outputHmmFile)

    def testMarginAlignLoadCustomInputModel(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--inputModel %s" % self.inputHmmFile)

    def testMarginAlignBwa(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--bwa")
        
    def testMarginAlignBwaNoRealign(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--bwa --noRealign")

def main():
    parseSuiteTestOptions()
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()
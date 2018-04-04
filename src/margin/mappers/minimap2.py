from margin.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class Minimap2(AbstractMapper):
    def run(self, args="-ax map-ont"):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because BWA builds these crufty index files, copy to a temporary directory
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile))
        system("minimap2 %s %s %s > %s" % (args, localReferenceFastaFile, self.readFastqFile, self.outputSamFile))

class Minimap2Chain(Minimap2):
    def run(self):
        Minimap2.run(self)
        self.chainSamFile()

class Minimap2Realign(Minimap2):
    def run(self):
        Minimap2.run(self)
        self.realignSamFile()


from sonLib.bioio import system
import os, sys, glob, time
from optparse import OptionParser
from margin.utils import samToBamFile
import pysam

class Fastaseq():
    """
    fasta reader
    """
    def __init__(self):
        self.id = None
        self.seq = ''
        self.length = ''
      
    @staticmethod 
    def readline(linein):
        seqobj = Fastaseq()
        for line in linein:
            if len(line)==0: 
                print >> sys.stderr, 'empty line'
                continue
            if line.startswith('>'):
                if seqobj.id is None:
                    seqobj.id = line.rstrip()
                    continue
                else:
                    yield seqobj
                    seqobj = Fastaseq()
                    seqobj.id = line.rstrip()
            else:
                seqobj.seq += line.rstrip('\n\r').upper()
        yield seqobj

def CustomTrackAssemblyHub(samFile, outputDir, hubTag, referenceFastaFile):
    """
    creates mapping.sorted.bam, mapping.sorted.bam.bai, moves them files to a folder;
    and creates trackDb.txt
    """
    parentFolder = outputDir + hubTag + "/"
    # Create output folder
    if not os.path.exists(parentFolder):
        os.mkdir(parentFolder)

    hubFastaDir = referenceFastaFile.split("/")[-1].split(".fa")[0]
    outFolderReferenceFiles = parentFolder + hubFastaDir + "/"
    outFolderBamFiles = outFolderReferenceFiles + "bamFiles/"

    # Create hierarchical reference and bamfile folders
    if not os.path.exists(outFolderReferenceFiles):
        os.mkdir(outFolderReferenceFiles)
    if not os.path.exists(outFolderBamFiles):
        os.mkdir(outFolderBamFiles)
    
    # Check and create bam, sorted bam, and indexed bam files
    bamFile = samFile.split(".sam")[0] + ".bam"
    sortedbamFile = outFolderBamFiles + samFile.split("/")[-1].split(".sam")[0] + ".sorted"
    if not os.path.isfile(bamFile):
        try:
            samToBamFile(samFile, bamFile)
        except:
            print >> sys.stderr, "Invalid Sam file"
            sys.exit()

    pysam.sort(bamFile, sortedbamFile)
    pysam.index(sortedbamFile + ".bam")


    genomes = open(os.path.join(parentFolder, "genomes.txt"), "w")
    if referenceFastaFile.endswith(".fa") or referenceFastaFile.endswith(".fasta"):
        header = os.path.split(referenceFastaFile)[-1].split(".fa")[0]
        system("cp %s %s" % (referenceFastaFile, parentFolder + header + "/"))

        # Create 2bit referenceFastaFile
        oldreferenceFastaFile = os.path.split(referenceFastaFile)[-1]
        newreferenceFastaFile = os.path.join(parentFolder, header, oldreferenceFastaFile)
        ref2bit = os.path.join(parentFolder, header, os.path.split(referenceFastaFile)[-1].split(".fa")[0] + ".2bit")
        system("./scripts/faToTwoBit %s %s" % (newreferenceFastaFile, ref2bit))
        
        # Get reference length for coordinates
        fastaFile = open(newreferenceFastaFile, "r")
        for seq in Fastaseq.readline(fastaFile):
            id = seq.id.split(" ")[0].replace(">", "")
            coord = len(seq.seq)
        fastaFile.close()
        
        # Fasta referenceFastaFile name without .fasta
        genomes.write("genome " + header + "\n")
        genomes.write("trackDb " + header + "/trackDb.txt\n")
        genomes.write("groups groups.txt\n")
        genomes.write("description " + header + "\n")
        genomes.write("twoBitPath " + header + "/" + header + ".2bit\n")
        genomes.write("organism " + header + "\n")
        genomes.write("defaultPos " + id + ":1-" + str(coord) + "\n")
        genomes.write("\n")
    
    track_label = 1
    tracks = open(os.path.join(outFolderReferenceFiles, "trackDb.txt"), "w")
    label = os.path.split(samFile)[-1]
    readType = os.path.split(samFile)[-1]
    tracks.write("track " + str(track_label) + "_\n")
    tracks.write("longLabel " + readType + "\n")
    shortLabel = readType
    tracks.write("shortLabel " + readType + "_" + shortLabel + "\n")
    tracks.write("priority 10\n")
    tracks.write("visibility pack\n")
    tracks.write("colorByStrand 150,100,30 230,170,40\n")
    tracks.write("color 150,100,30\n")
    tracks.write("altColor 230,170,40\n")
    tracks.write("bigDataUrl bamFiles/" + os.path.split(samFile)[-1].split(".sam")[0] + ".sorted.bam\n")
    tracks.write("type bam\n")
    tracks.write("group " + readType + "\n")
    tracks.write("html assembly\n\n")
    track_label += 1
    tracks.close()

    genomes.close()

    groups = open(os.path.join(parentFolder, "groups.txt"), "w")
    groups.write("name samFile\n")
    groups.write("label readType\n")
    groups.write("priority 1\n")
    groups.write("defaultIsClosed 0\n\n")
    groups.close()
    
    
    hubFile = open(os.path.join(parentFolder, "hub.txt"), "w")
    hubFile.write("hub marginAlignHub\n")
    hubFile.write("shortLabel " + samFile + "\n")
    hubFile.write("longLabel marginAlign Assembly Hub from UCSC Nanopore Group\n")
    hubFile.write("genomesFile genomes.txt\n")
    hubFile.write("email miten@soe.ucsc.edu")
    hubFile.close()
    
    
    htmlFile = open(parentFolder + "/index.html", "w")
    htmlFile.write("<!DOCTYPE HTML PUBLIC\"-//W3C//DTD HTML 3.2 Final//EN\">\n")
    htmlFile.write("<html>")
    htmlFile.write("<head>")
    htmlFile.write("<title>Index of Assembly Hub</title>")
    htmlFile.write("</head>")
    htmlFile.write("<body>")
    htmlFile.write("<h1>Index of Assembly Hub</h1>")

    htmlFile.write("<table><tr><th><img src=\"/icons/blank.gif\" alt=\"[ICO]\"></th><th><a href=\"?C=N;O=D\">Name</a></th><th><a href=\"?C=M;O=A\">Last modified</a></th><th><a href=\"?C=S;O=A\">Size</a></th><th><a href=\"?C=D;O=A\">Description</a></th></tr><tr><th colspan=\"5\"><hr></th></tr>")
    htmlFile.write("<tr><td valign=\"top\"><img src=\"/icons/back.gif\" alt=\"[DIR]\"></td><td><a href=\"../\">Parent Directory</a></td><td>&nbsp;</td><td align=\"right\">  - </td><td>&nbsp;</td></tr>")
    htmlFile.write("<tr><td valign=\"top\"><img src=\"/icons/text.gif\" alt=\"[TXT]\"></td><td><a href=\"hub.txt\">hub.txt</a></td><td align=\"right\">13-April-2015 12:00  </td><td align=\"right\">186 </td><td>&nbsp;</td></tr>")
    htmlFile.write("<tr><td valign=\"top\"><img src=\"/icons/text.gif\" alt=\"[TXT]\"></td><td><a href=\"genomes.txt\">genomes.txt</a></td><td align=\"right\">13-April-2015 12:00  </td><td align=\"right\">806 </td><td>&nbsp;</td></tr>")
    htmlFile.write("<tr><td valign=\"top\"><img src=\"/icons/text.gif\" alt=\"[TXT]\"></td><td><a href=\"groups.txt\">groups.txt</a></td><td align=\"right\">13-April-2015 12:00  </td><td align=\"right\">448 </td><td>&nbsp;</td></tr>")

    htmlFile.write("<tr><td valign=\"top\"><img src=\"/icons/folder.gif\" alt=\"[DIR]\"></td><td><a href=" + str(hubFastaDir) + ">" + hubFastaDir + "/</a></td><td align=\"right\">13-April-2015 12:00  </td><td align=\"right\">  - </td><td>&nbsp;</td></tr>")
    htmlFile.write("<tr><th colspan=\"5\"><hr></th></tr>")
    htmlFile.write("</table>")
    htmlFile.write("</body></html>")
    
    htmlFile.close()

##################################################################################
# Main
# Here is the main program
##################################################################################

def main(myCommandLine=None):
    # starting time
    t0 = time.time()

    #Parse the inputs args/options
    parser = OptionParser(usage="usage: inputSam/BamFile outputDir hubTag \
                                 referenceFastaFile", version="%prog 0.1")

    #Parse the options/arguments
    options, args = parser.parse_args()

    #Exit if the arguments are not what we expect
    if len(args) != 4:
        raise RuntimeError("Expected four arguments, got: %s" % " ".join(args))

    samFile = args[0]
    outputDir = args[1]
    hubTag = args[2]
    referenceFastaFile = args[3]
    
    print >> sys.stderr, samFile, outputDir, hubTag, referenceFastaFile

    CustomTrackAssemblyHub(samFile, outputDir, hubTag, referenceFastaFile)

    print >> sys.stderr, "\n", "Total time for the program %.3f" % (time.time()-t0), "s"

if (__name__ == "__main__"):
    main()
    raise SystemExit

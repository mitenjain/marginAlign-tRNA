The marginAlign package can be used to align reads to a reference genome and call single nucleotide variations (SNVs). It is specifically tailored for Oxford Nanopore Reads.

The package comes with two programs: marginAlign, a short read aligner, and marginVar, a program to call single nucleotide variations.

### Requirements
* git
* python 2.7


### Installation
To install the code run:

    git clone git://github.com/mitenjain/nanopore.git
    cd nanopore
    git pull
    git submodule update --init
    make

This will build the code. Two executables: "marginAlign" and "marginVar" are in the base directory
of the package. Place these binaries on your path if you wish to use them outside of the base directory.

### Testing
To test the installation run:

    make test
    
This will run demo sequences through marginAlign and marginVar.
    
### Updating the installation
To update a marginAlign installation, from the base directory type:

    git pull
    git submodule update --init
    make clean
    make

### Running marginAlign

To access the help message type:

    marginAlign --help

This will give you a full list of options. Most are related to the (jobTree)[https://github.com/benedictpaten/jobTree] options, which control how the script is executed. 

To align a FASTQ file ("input.fastq") to a reference fasta file ("reference.fasta") and output the alignment in SAM format with marginAlign ("output.sam") do:

    marginAlign input.fastq reference.fasta output.sam

To enable EM training do, putting the trained model file in "output.hmm" do:

    marginAlign input.fastq reference.fasta output.sam --em output.hmm

To use a pretrained model "input.hmm" do:

    marginAlign input.fastq reference.fasta output.sam --useModel input.hmm

### Running marginVar

To call single nucleotide variations from an existing alignment ("input.sam") sam file and the output ("output.vcf"):

    marginVar input.sam output.vcf 

To NOT marginalise over the read alignments do:

    marginVar input.sam output.vcf --noMargin

### Citing marginAlign/marginVar

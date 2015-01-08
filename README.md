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
of the package. Place these binaries on your path if you wish to use them without referring to their absolute
path on the filesystem.

### Creating a virtual environment to handle dependencies
marginAlign uses numpy and pysam. The system python coule be used if have these dependencies are present. Otherwise, a virtual environment can be created by running:

    virtualenv --no-site-packages --distribute env && source env/bin/activate && pip install -r requirements.txt

Creation of a virtual environment requires that machine has an existing installation of pip [https://pip.pypa.io/en/latest/index.html] and virtualenv [https://virtualenv.pypa.io/en/latest/index.html].

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

### jobTree

Both marginAlign and marginVar are (jobTree)[https://github.com/benedictpaten/jobTree] scripts. This allows you to parallelise the execution of the scripts either using a single multi-processor machine or a cluster. The most important thing to note is that a jobTree script creates a "jobTree", a directory of files that manages the running of the job. As a result both scripts take a --jobTree DIRECTORY argument (by default "jobTree" in the current working directory if not specified). This directory should not exist prior to executing the script. It will be not be deleted after the script has finished (as it can be used to interogate the details of a run), and should be manually deleted.

### Running marginAlign

To access the help message type:

    marginAlign --help

or just:

    marginAlign

This will give you a full list of options. Most are related to the (jobTree)[https://github.com/benedictpaten/jobTree] options, which control how the script is executed. 

To align a FASTQ file ("input.fastq") to a reference fasta file ("reference.fasta") and output the alignment in SAM format with marginAlign ("output.sam") using the "./jobTree" jobTree to manage the batch do:

    marginAlign input.fastq reference.fasta output.sam --jobTree ./jobTree

After executing the "./jobTree" directory should be deleted.

To enable EM training do, putting the trained model file in "output.hmm" do:

    marginAlign input.fastq reference.fasta output.sam --em --outputModel output.hmm --jobTree ./jobTree

To use a different model to the default one do:

    marginAlign input.fastq reference.fasta output.sam --inputModel input.hmm --jobTree ./jobTree

### Running marginVar

To call single nucleotide variations from an existing alignment ("input.sam") sam file and the output ("output.vcf"):

    marginVar input.sam output.vcf --jobTree ./jobTree

To NOT marginalise over the read alignments do:

    marginVar input.sam output.vcf --noMargin --jobTree ./jobTree

### Citing marginAlign/marginVar

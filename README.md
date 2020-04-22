# NGSBuilder
Pipeline code and associated python scripts for generation of HIV-consensus sequences from Illumina MiSeq

N.B. there are multiple lines in the build_contig_from_nextera.sh file that are system-specific and should be changed to suit the file-hierachy of the machine. These include the path to both auxiliary sequence files and path to dependencies.

# Dependencies
Pipeline requires the following standard dependencies to work
 - bwa
 - samtools
 - trimmomatic
 - picard
 - python
 
The following dependencies are required in the building step:
  - SPAdes (Bankevich et al. 2012)
  - SHIVER (Wymant et al. 2018)
  
The following pre-steps are required:
  - an indexed hg19/hg38 file for BWA
    -- should also contain 3 HIV strains for reference as well
  - a set of HIV reference sequences for SHIVER
  - a fasta file containing PCR primer sequences
  - a fasta file containing MiSeq adapter sequences
  - a fasta file combining the two above
  
# Usage:
The pipeline usage is as follows:
$build_contig_from_nextera.sh $WORKING_DIR $FASTQ_FOLDER $PATH_TO_HIV_REFERENCES

# Auxiliary Scripts
The pipeline contains 3 auxiliary python scripts that are optional for handling output sequence data
 - cullConsensi.py: grab all consensus sequences from folders containing fasta files with consensus sequences from builder
 - rename_seqs.py: input yaml file of python dictionary to rename sequences
 - bff_parse.py: re-generate sequence from SHIVER base frequencies file using ambiguous bases and custom read-depth parameters
 
 
 
 All bugs/issues/questions may be directed at michaeljbale93@gmail.com

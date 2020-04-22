#!/usr/bin/env bash

###CONSTANTS###

SHIVER_CONFIG_FILE="/home/balemj/NGSTools/full_length_seq_pipeline/shiver-master/config_bwa.sh"
AUX_SEQUENCE_FOLDER="/home/balemj/SGSScript/dir_assem/dir_NFL_aux_files"
AUX_NEXTERA_PRIMERS="nextera.fa"	
AUX_PCR_PRIMERS="pcr_primers.fas"	
#AUX_HIV_REFERENCES="AlignSeq.nt.fasta"	
AUX_ALL_ADAPTERS="all_adapters.fas"	
AUX_HUMAN_GENOME="hg_ref/hg19.fa"	
PATH_TO_SPADES="/home/balemj/NGSTools/full_length_seq_pipeline/SPAdes-3.13.1-Linux/bin"
PATH_TO_SHIVER="/home/balemj/NGSTools/full_length_seq_pipeline/shiver-master"
PATH_TO_TRIMMOMATIC="/mnt/nasapps/production/trimmomatic/0.38/bin"
PATH_TO_BWA="/mnt/nasapps/development/bwa/0.7.17/bin"
PATH_TO_SAMTOOLS="/mnt/nasapps/production/samtools/1.8/bin"
PATH_TO_PICARD="/mnt/nasapps/production/picard/2.18.26/bin"

WORKING_DIR=$2
AUX_HIV_REFERENCES=$3
############################


#SETUP#

cd ${1}
cd $WORKING_DIR

echo "CREATING TEMPORARY WORKSPACE IN FOLDER " $WORKING_DIR
mkdir raw_reads trimmed_reads new_reads temp_aux_files
cp ${AUX_SEQUENCE_FOLDER}/${AUX_ALL_ADAPTERS} temp_aux_files
r1="read1_index_${WORKING_DIR}.fq.gz"
r2="read2_index_${WORKING_DIR}.fq.gz"
mv $r1 raw_reads
mv $r2 raw_reads


######QC BLOCK#######################
echo "TRIMMING ADAPTER SEQUENCES WITH TRIMMOMATIC"

${PATH_TO_TRIMMOMATIC}/trimmomatic PE \
    raw_reads/${r1} \
    raw_reads/${r2} \
    trimmed_reads/${r1/.fq.gz}.paired.trimmed.fq.gz \
    trimmed_reads/${r1/.fq.gz}.unpaired.trimmed.fq.gz \
    trimmed_reads/${r2/.fq.gz}.paired.trimmed.fq.gz \
    trimmed_reads/${r2/.fq.gz}.unpaired.trimmed.fq.gz \
    ILLUMINACLIP:temp_aux_files/${AUX_ALL_ADAPTERS}:2:10:7:1:true \
    MINLEN:50 \
    SLIDINGWINDOW:4:20

r1="${r1/.fq.gz}.paired.trimmed.fq.gz"
r2="${r2/.fq.gz}.paired.trimmed.fq.gz"


echo "REMOVING CONTAMINANT READS WITH BWA"

${PATH_TO_BWA}/bwa mem \
    ${AUX_SEQUENCE_FOLDER}/${AUX_HUMAN_GENOME} \
    trimmed_reads/${r1} \
    trimmed_reads/${r2} \
    | \
    ${PATH_TO_SAMTOOLS}/samtools view \
    -Sb \
    -f 4 \
    -o new_reads/reads12_index_${WORKING_DIR}.trimmed.unmapped.bam \
    -U new_reads/reads12_index_${WORKING_DIR}.trimmed.mapped.bam

${PATH_TO_SAMTOOLS}/samtools sort \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.mapped.bam \
    -o new_reads/reads12_index_${WORKING_DIR}.trimmed.sort_mapped.bam

${PATH_TO_SAMTOOLS}/samtools index \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.sort_mapped.bam

${PATH_TO_SAMTOOLS}/samtools view \
    -Sb \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.sort_mapped.bam \
    chrB \
    chrC \
    chrH \
    -o new_reads/reads12_index_${WORKING_DIR}.trimmed.HIVMapped.bam \
    -U new_reads/reads12_index_${WORKING_DIR}.trimmed.HGMapped.bam 

${PATH_TO_SAMTOOLS}/samtools view \
    -Sb \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.HGMapped.bam \
    -q 21 \
    -o new_reads/reads12_index_${WORKING_DIR}.trimmed.Contam.bam \
    -U new_reads/reads12_index_${WORKING_DIR}.trimmed.lowQ.bam

${PATH_TO_SAMTOOLS}/samtools merge \
    new_reads/out.bam \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.unmapped.bam \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.HIVMapped.bam \
    new_reads/reads12_index_${WORKING_DIR}.trimmed.lowQ.bam

java -jar ${PATH_TO_PICARD}/picard.jar \
    SamToFastq \
    I=new_reads/out.bam \
    F=new_reads/${r1/.gz} \
    F2=new_reads/${r2/.gz} \
    VALIDATION_STRINGENCY=SILENT


for i in new_reads/*fq
do
    /home/balemj/SGSScript/dir_script/make_unique_fastq.py $i ${i/.fq}.forDN.unique.fq
    gzip ${i/.fq}.forDN.unique.fq
done

r1="${r1/.fq.gz}.forDN.unique.fq"
r2="${r2/.fq.gz}.forDN.unique.fq"

rm -rf trimmed_reads
rm new_reads/*.bam

#####################



########DE NOVO BUILD##############

echo "BUILD CONTIGS WITH SPAdes"

${PATH_TO_SPADES}/spades.py -k 21,33,55,77,99,127 -1 new_reads/${r1}.gz -2 new_reads/${r2}.gz -o temp_${WORKING_DIR}_SPAdes_OUTPUT_DIR

mv temp_${WORKING_DIR}_SPAdes_OUTPUT_DIR/contigs.fasta .
rm -rf temp_${WORKING_DIR}_SPAdes_OUTPUT_DIR

###################################

#############SHIVER BUILD###############

echo "REFINE BUILD WITH SHIVER"

echo "INITIALIZING SHIVER"

${PATH_TO_SHIVER}/shiver_init.sh shiver_init_dir ${SHIVER_CONFIG_FILE} ${AUX_HIV_REFERENCES} ${AUX_SEQUENCE_FOLDER}/${AUX_NEXTERA_PRIMERS} ${AUX_SEQUENCE_FOLDER}/${AUX_PCR_PRIMERS}

echo "BEGIN SHIVER STEP 1: SHIVER_ALIGN_CONTIGS.SH"

${PATH_TO_SHIVER}/shiver_align_contigs.sh shiver_init_dir ${SHIVER_CONFIG_FILE} contigs.fasta index_${WORKING_DIR}_shiverBuild

echo "BEGIN SHIVER STEP 2: SHIVER_MAP_READS.SH"
${PATH_TO_SHIVER}/shiver_map_reads.sh shiver_init_dir ${SHIVER_CONFIG_FILE} contigs.fasta index_${WORKING_DIR}_shiverBuild \
    index_${WORKING_DIR}_shiverBuild.blast index_${WORKING_DIR}_shiverBuild_cut_wRefs.fasta new_reads/${r1}.gz new_reads/${r2}.gz

#############################

############CLEAN WORKSHPACE################

echo 'BUILD FINISHED; CLEANING WORKSPACE'
rm -rf temp*
rm contigs.fasta
mkdir shiver_files
mv index_${WORKING_DIR}_shiverBuild* shiver_files
mv shiver_files/index_${WORKING_DIR}_shiverBuild_consensus_MinCov_15_30.fasta .
mv shiver_files/index_${WORKING_DIR}_shiverBuild_BaseFreqs.csv .
mv raw_reads/* .
rm -rf */
rm *fq
cd ..


#TODO: Make log code:
#read loss
#conensus length with and wo N's
#    TODO: make gen_consensus.py
#ratio of lower:upper
#number of mixture bases
#Avg read depth

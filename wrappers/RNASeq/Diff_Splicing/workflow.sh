#!/bin/sh
#SBATCH -N 1 #Number of nodes
#SBATCH -n 12 # number of cores
#SBATCH --mem 48G # memory pool for all cores
#SBATCH -t 03-12:00 # time (D-HH:MM)
#SBATCH -o mRNA-IKAROS.%N.%j.out
#SBATCH -e mRNA-IKAROS.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --account=rrg-majewski-ab
#SBATCH --mail-user=erm.ba@yandex.com
#
### Global variables
#
RNA_HOME="/project/6007495/barekeer/projects/bareke/milot_ikaros"
RNA_REF_FASTA="/project/6007495/singularity/1.1.0/builtin/genomes/Mmusculus/mm10/seq/mm10.fa"
RNA_REF_INDEX="/project/6007495/singularity/1.1.0/builtin/genomes/Mmusculus/mm10/hisat2/mm10"
RNA_REF_GTF="/project/6007495/singularity/1.1.0/builtin/genomes/Mmusculus/mm10/rnaseq/ref-transcripts.gtf"
RNA_DATA_DIR="${RNA_HOME}/reads"
RNA_DATA_TRIM_DIR="${RNA_HOME}/reanalysis/reads/trimmed"
RNA_REFS_DIR="${RNA_HOME}/reanalysis/refs/adaptors"
RNA_RUN_DIR="${RNA_HOME}/reanalysis/runs"
RNA_ALIGN_DIR="${RNA_HOME}/reanalysis/reads/bams"
RNA_EXPRESSION_DIR="${RNA_HOME}/reanalysis/output/expression"
#
### Create trimming folder
mkdir -p $RNA_DATA_TRIM_DIR
### Create reference folder
mkdir -p $RNA_REFS_DIR
#cd $RNA_REFS_DIR
#wget http://genomedata.org/rnaseq-tutorial/illumina_multiplex.fa
### Trim raw reads
mkdir -p $RNA_RUN_DIR
cd $RNA_RUN_DIR
##
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/WT1-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/WT1-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/WT1-OP9
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/WT3-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/WT3-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/WT3-OP9
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/WT4-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/WT4-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/WT4-OP9
#
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/IK1-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/IK1-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/IK1-OP9
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/IK2-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/IK2-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/IK2-OP9
flexbar --adapter-min-overlap 7 --adapter-trim-end RIGHT --adapters $RNA_REFS_DIR/illumina_multiplex.fa --pre-trim-left 13 --max-uncalled 300 --min-read-length 25 --threads 12 --zip-output GZ --reads $RNA_DATA_DIR/IK3-OP9_R1.fastq.gz --reads2 $RNA_DATA_DIR/IK3-OP9_R2.fastq.gz --target $RNA_DATA_TRIM_DIR/IK3-OP9
### Post QC
#cd $RNA_DATA_TRIM_DIR
#fastqc *.fastq.gz
##
### Align trimmed data
mkdir -p $RNA_ALIGN_DIR
cd $RNA_ALIGN_DIR
#
hisat2 -p 12 --rg-id=WT1-OP9 --rg SM:WT --rg LB:WT1-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-ACTGAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/WT1-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/WT1-OP9_2.fastq.gz -S ./WT1-OP9.sam
hisat2 -p 12 --rg-id=WT3-OP9 --rg SM:WT --rg LB:WT3-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/WT3-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/WT3-OP9_2.fastq.gz -S ./WT3-OP9.sam
hisat2 -p 12 --rg-id=WT4-OP9 --rg SM:WT --rg LB:WT4-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-CTGACA.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/WT4-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/WT4-OP9_2.fastq.gz -S ./WT4-OP9.sam
#
hisat2 -p 12 --rg-id=IK1-OP9 --rg SM:IK --rg LB:IK1-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-TGACAC.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/IK1-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/IK1-OP9_2.fastq.gz -S ./IK1-OP9.sam
hisat2 -p 12 --rg-id=IK2-OP9 --rg SM:IK --rg LB:IK2-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-GACACT.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/IK2-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/IK2-OP9_2.fastq.gz -S ./IK2-OP9.sam
hisat2 -p 12 --rg-id=IK3-OP9 --rg SM:IK --rg LB:IK3-OP9_IKAROS --rg PL:ILLUMINA --rg PU:CXX1234-ACACTG.1 -x $RNA_REF_INDEX --dta --rna-strandness RF -1 $RNA_DATA_TRIM_DIR/IK3-OP9_1.fastq.gz -2 $RNA_DATA_TRIM_DIR/IK3-OP9_2.fastq.gz -S ./IK3-OP9.sam
##
### SAM to BAM conversion
samtools sort -@ 12 -o WT1-OP9.bam WT1-OP9.sam
samtools sort -@ 12 -o WT3-OP9.bam WT3-OP9.sam 
samtools sort -@ 12 -o WT4-OP9.bam WT4-OP9.sam 
samtools sort -@ 12 -o IK1-OP9.bam IK1-OP9.sam
samtools sort -@ 12 -o IK2-OP9.bam IK2-OP9.sam
samtools sort -@ 12 -o IK3-OP9.bam IK3-OP9.sam
### Index BAM files
find *.bam -exec echo samtools index {} \; | sh
### BAM read counting
# mkdir bam_readcount
# cd bam_readcount
### Create faidx indexed reference sequence file for use with mpileup 
# samtools faidx $RNA_REF_FASTA
### Create mpileup files
# samtools mpileup -f $RNA_REF_FASTA $RNA_ALIGN_DIR/WT.bam $RNA_ALIGN_DIR/IK.bam
### BAM counts
# bam-readcount -f $RNA_REF_FASTA $RNA_ALIGN_DIR/WT.bam 2>/dev/null 1>WT_bam-readcounts.txt
# bam-readcount -f $RNA_REF_FASTA $RNA_ALIGN_DIR/IK.bam 2>/dev/null 1>IK_bam-readcounts.txt
### Parse read-counts per base
# cat WT_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "WT Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";'
# cat IK_bam-readcounts.txt | perl -ne '@data=split("\t", $_); @Adata=split(":", $data[5]); @Cdata=split(":", $data[6]); @Gdata=split(":", $data[7]); @Tdata=split(":", $data[8]); print "IK Counts\t$data[0]\t$data[1]\tA: $Adata[1]\tC: $Cdata[1]\tT: $Tdata[1]\tG: $Gdata[1]\n";'
#
###
mkdir -p ${RNA_EXPRESSION_DIR}/stringtie/ref_only/
cd ${RNA_EXPRESSION_DIR}/stringtie/ref_only/
#
stringtie -p 12 -G $RNA_REF_GTF -e -B -o IK1-OP9/transcripts.gtf -A IK1-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/IK1-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -e -B -o IK2-OP9/transcripts.gtf -A IK2-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/IK2-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -e -B -o IK3-OP9/transcripts.gtf -A IK3-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/IK3-OP9.bam
#
stringtie -p 12 -G $RNA_REF_GTF -e -B -o WT1-OP9/transcripts.gtf -A WT1-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/WT1-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -e -B -o WT3-OP9/transcripts.gtf -A WT3-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/WT3-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -e -B -o WT4-OP9/transcripts.gtf -A WT4-OP9/gene_abundances.tsv $RNA_ALIGN_DIR/WT4-OP9.bam
#
stringtie_expression_matrix.pl --expression_metric=TPM --result_dirs='IK1-OP9,IK2-OP9,IK3-OP9,WT1-OP9,WT3-OP9,WT4-OP9' --transcript_matrix_file=transcript_tpm_all_samples.tsv --gene_matrix_file=gene_tpm_all_samples.tsv
stringtie_expression_matrix.pl --expression_metric=FPKM --result_dirs='IK1-OP9,IK2-OP9,IK3-OP9,WT1-OP9,WT3-OP9,WT4-OP9' --transcript_matrix_file=transcript_fpkm_all_samples.tsv --gene_matrix_file=gene_fpkm_all_samples.tsv
stringtie_expression_matrix.pl --expression_metric=Coverage --result_dirs='IK1-OP9,IK2-OP9,IK3-OP9,WT1-OP9,WT3-OP9,WT4-OP9' --transcript_matrix_file=transcript_coverage_all_samples.tsv --gene_matrix_file=gene_coverage_all_samples.tsv
##
### Ref-only
mkdir -p ${RNA_EXPRESSION_DIR}/htseq_counts
cd ${RNA_EXPRESSION_DIR}/htseq_counts
#
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/WT1-OP9.bam $RNA_REF_GTF > WT1-OP9_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/WT3-OP9.bam $RNA_REF_GTF > WT3-OP9_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/WT4-OP9.bam $RNA_REF_GTF > WT4-OP9_gene.tsv
#
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/IK1-OP9.bam $RNA_REF_GTF > IK1-OP9_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/IK2-OP9.bam $RNA_REF_GTF > IK2-OP9_gene.tsv
htseq-count --format bam --order pos --mode intersection-strict --stranded reverse --minaqual 1 --type exon --idattr gene_id $RNA_ALIGN_DIR/IK3-OP9.bam $RNA_REF_GTF > IK3-OP9_gene.tsv
##
join WT1-OP9_gene.tsv WT3-OP9_gene.tsv | join - WT4-OP9_gene.tsv | join - IK1-OP9_gene.tsv | join - IK2-OP9_gene.tsv | join - IK3-OP9_gene.tsv > gene_read_counts_table_all.tsv
echo "GeneID WT1-OP9 WT3-OP9 WT4-OP9 IK1-OP9 IK2-OP9 IK3-OP9" > header.txt
cat header.txt gene_read_counts_table_all.tsv | grep -v "__" | perl -ne 'chomp $_; $_ =~ s/\s+/\t/g; print "$_\n"' > gene_read_counts_table_all_final.tsv
rm -f gene_read_counts_table_all.tsv header.txt
head gene_read_counts_table_all_final.tsv
##
#### Ref-guided
#
mkdir -p ${RNA_EXPRESSION_DIR}/stringtie/ref_guided/
cd ${RNA_EXPRESSION_DIR}/stringtie/ref_guided/
#
stringtie -p 12 -G $RNA_REF_GTF -l IK1-OP9 -o IK1-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK1-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -l IK2-OP9 -o IK2-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK2-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -l IK3-OP9 -o IK3-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK3-OP9.bam
#
stringtie -p 12 -G $RNA_REF_GTF -l WT1-OP9 -o WT1-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT1-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -l WT3-OP9 -o WT3-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT3-OP9.bam
stringtie -p 12 -G $RNA_REF_GTF -l WT4-OP9 -o WT4-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT4-OP9.bam
#
#### Do novo
#
mkdir -p ${RNA_EXPRESSION_DIR}/stringtie/de_novo/
cd ${RNA_EXPRESSION_DIR}/stringtie/de_novo/
#
stringtie -p 12 -l IK1-OP9 -o IK1-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK1-OP9.bam
stringtie -p 12 -l IK2-OP9 -o IK2-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK2-OP9.bam
stringtie -p 12 -l IK3-OP9 -o IK3-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK3-OP9.bam
#
stringtie -p 12 -l WT1-OP9 -o WT1-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT1-OP9.bam
stringtie -p 12 -l WT3-OP9 -o WT3-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT3-OP9.bam
stringtie -p 12 -l WT4-OP9 -o WT4-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT4-OP9.bam
################## Merging GTF files ##################
#
### Merge Ref-guided ###
#
cd ${RNA_EXPRESSION_DIR}/stringtie/ref_guided/
ls -1 *OP9*/transcripts.gtf > assembly_GTF_list.txt
# cat assembly_GTF_list.txt
stringtie --merge -p 12 -o stringtie_merged.gtf -G $RNA_REF_GTF assembly_GTF_list.txt
gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
awk '{if($3=="transcript") print}' stringtie_merged.gtf | cut -f 1,4,9 > stringtie_merged.txt
awk '{if($3=="transcript") print}' gffcompare.annotated.gtf | cut -f 1,4,9 > gffcompare.annotated.txt
##
### Merge De novo
cd ${RNA_EXPRESSION_DIR}/stringtie/de_novo/
ls -1 *OP9*/transcripts.gtf > assembly_GTF_list.txt
#cat assembly_GTF_list.txt
stringtie --merge -p 12 -o stringtie_merged.gtf assembly_GTF_list.txt
gffcompare -r $RNA_REF_GTF -o gffcompare stringtie_merged.gtf
# cat gffcompare.stats
##
##### Differential (Expression) Splicing ####
##
# Ballgown tables : guided merged GTF
#
cd ${RNA_EXPRESSION_DIR}/stringtie/
mkdir -p ref_guided_merged
cd ref_guided_merged
#
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o IK1-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK1-OP9.bam
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o IK2-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK2-OP9.bam
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o IK3-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK3-OP9.bam
#
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o WT1-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT1-OP9.bam
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o WT3-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT3-OP9.bam
stringtie -p 12 -G ../ref_guided/stringtie_merged.gtf -e -B -o WT4-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT4-OP9.bam
##
mkdir -p ${RNA_EXPRESSION_DIR}/de/ballgown/ref_guided_merged/
cd ${RNA_EXPRESSION_DIR}/de/ballgown/ref_guided_merged/
#
printf "\"ids\",\"type\",\"path\"\n\"WT1-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/WT1-OP9\"\n\"WT3-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/WT3-OP9\"\n\"WT4-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/WT4-OP9\"\n\"IK1-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/IK1-OP9\"\n\"IK2-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/IK2-OP9\"\n\"IK3-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/IK3-OP9\"\n" > WT_vs_IK.csv
##
#
# Ballgown table: De novo merged
#
cd ${RNA_EXPRESSION_DIR}/stringtie/
mkdir -p de_novo_merged
cd de_novo_merged
#
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o IK1-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK1-OP9.bam
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o IK2-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK2-OP9.bam
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o IK3-OP9/transcripts.gtf $RNA_ALIGN_DIR/IK3-OP9.bam
#
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o WT1-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT1-OP9.bam
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o WT3-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT3-OP9.bam
stringtie -p 12 -G ../de_novo/stringtie_merged.gtf -e -B -o WT4-OP9/transcripts.gtf $RNA_ALIGN_DIR/WT4-OP9.bam
##
mkdir -p ${RNA_EXPRESSION_DIR}/de/ballgown/de_novo_merged/
cd ${RNA_EXPRESSION_DIR}/de/ballgown/de_novo_merged/
#
printf "\"ids\",\"type\",\"path\"\n\"WT1-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/WT1-OP9\"\n\"WT3-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/WT3-OP9\"\n\"WT4-OP9\",\"WT\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/WT4-OP9\"\n\"IK1-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/IK1-OP9\"\n\"IK2-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/IK2-OP9\"\n\"IK3-OP9\",\"IK\",\"${RNA_EXPRESSION_DIR}/stringtie/de_novo_merged/IK3-OP9\"\n" > WT_vs_IK.csv
##
## How many genes have at least one potentially novel transcript assembled? ##
cd ${RNA_EXPRESSION_DIR}/stringtie/de_novo/
# 
grep "j" gffcompare.stringtie_merged.gtf.tmap | cut -f 1 | sort | uniq > ${RNA_EXPRESSION_DIR}/de/ballgown/de_novo_merged/genes_with_novel_transcript_assembled.txt
# How many ?
# wc -l genes_with_novel_transcript_assembled.txt
grep -w "u" gffcompare.stringtie_merged.gtf.tmap | sort -n -k 10 | column -t > ${RNA_EXPRESSION_DIR}/de/ballgown/de_novo_merged/transcripts_in_intergenic_regions__candidate_novel_regions_of_transcription.txt
#
############## Annotate all individual splice junctions #######################
cd $RNA_ALIGN_DIR
#
## Merge IK and WT replicates into merged BAM files
picard -Xmx16g MergeSamFiles OUTPUT=WT.bam INPUT=WT1-OP9.bam INPUT=WT3-OP9.bam INPUT=WT4-OP9.bam
picard -Xmx16g MergeSamFiles OUTPUT=IK.bam INPUT=IK1-OP9.bam INPUT=IK2-OP9.bam INPUT=IK3-OP9.bam
samtools index WT.bam
samtools index IK.bam
##
regtools junctions extract WT.bam > WT.junctions.bed
regtools junctions annotate WT.junctions.bed $RNA_REF_FASTA $RNA_REF_GTF > WT.junctions.anno.bed
#
regtools junctions extract IK.bam > IK.junctions.bed
regtools junctions annotate IK.junctions.bed $RNA_REF_FASTA $RNA_REF_GTF > IK.junctions.anno.bed
#
## Pull out any junctions from either sample that appear to involve novel exon skipping, acceptor site usage, or donor site usage
#
mkdir -p ${RNA_EXPRESSION_DIR}/splice_events
cd ${RNA_EXPRESSION_DIR}/splice_events/
#
grep -P -w "NDA|A|D" ${RNA_ALIGN_DIR}/WT.junctions.anno.bed | perl -ne 'chomp; @l=split("\t",$_); if ($l[4] > 5){print "$_\n"}' > WT.splice.events.txt
grep -P -w "NDA|A|D" ${RNA_ALIGN_DIR}/IK.junctions.anno.bed | perl -ne 'chomp; @l=split("\t",$_); if ($l[4] > 5){print "$_\n"}' > IK.splice.events.txt
#
#### Organize illustrative GTF files ####
cd ${RNA_EXPRESSION_DIR}/stringtie/ref_only/
stringtie_filter_gtf.pl --expression_metric=FPKM --result_dirs='WT1-OP9,WT3-OP9,WT3-OP9,IK1-OP9,IK2-OP9,IK3-OP9' --input_gtf_file='${RNA_REF_GTF}' --filtered_gtf_file='${RNA_EXPRESSION_DIR}/stringtie/ref_only/stringtie_merged.filtered.gtf' --exp_cutoff=0 --min_sample_count=2
#
cd ${RNA_EXPRESSION_DIR}/stringtie/ref_guided_merged/
stringtie_filter_gtf.pl --expression_metric=FPKM --result_dirs='WT1-OP9,WT3-OP9,WT3-OP9,IK1-OP9,IK2-OP9,IK3-OP9' --input_gtf_file='${RNA_EXPRESSION_DIR}/stringtie/ref_guided/stringtie_merged.gtf' --filtered_gtf_file='${RNA_EXPRESSION_DIR}/stringtie/ref_guided/stringtie_merged.filtered.gtf' --exp_cutoff=0 --min_sample_count=2
#
## Rename the GTF files for visualization ###
#
mkdir -p ${RNA_EXPRESSION_DIR}/stringtie/visualization/
cd ${RNA_EXPRESSION_DIR}/stringtie/visualization/
#
cat ${RNA_REF_GTF} | perl -ne 'chomp; @l=split("\t", $_); print "$_\n" unless ($l[2] eq "gene");' > reference.gtf
cp ${RNA_EXPRESSION_DIR}/stringtie/ref_only/stringtie_merged.filtered.gtf ref_only.gtf
cp ${RNA_EXPRESSION_DIR}/stringtie/ref_guided/stringtie_merged.filtered.gtf ref_guided.gtf
cp ${RNA_EXPRESSION_DIR}/stringtie/de_novo/stringtie_merged.gtf de_novo.gtf
## Identify some candidate novel transcripts to visualize -- Sashimi plot ?? ###
#

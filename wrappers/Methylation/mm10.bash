#!/bin/bash
set -e
set -o pipefail
#
while [ "$#" -gt 0 ]; do
case "$1" in
-s) sample="$2"; shift 2;;
-1) read1="$2"; shift 2;;
-2) read2="$2"; shift 2;;
-o) outfolder="$2"; shift 2;;
-b) build="$2"; shift 2;;
-m) mincov="$2"; shift 2;;
-t) threads="$2"; shift 2;;

--sample=*) sample="${1#*=}"; shift 1;;
--read1=*) read1="${1#*=}"; shift 1;;
--read2=*) read2="${1#*=}"; shift 1;;
--outfolder=*) outfolder="${1#*=}"; shift 1;;
--build=*) build="${1#*=}"; shift 1;;
--mincov=*) mincov="${1#*=}"; shift 1;;
--threads=*) threads="${1#*=}"; shift 1;;

--sample|--read1|--read2|--outfolder|--build|--mincov|--threads) echo "$1 requires an argument" >&2; exit 1;;

-*) echo "unknown option: $1" >&2; exit 1;;
*) die "unrecognized argument: $1"; shift 1;;
esac
done
#
BASEDIR=`readlink -f "${0%/*}"`
#
refmain="/project/6007495/singularity/1.1.0/builtin/genomes/Mmusculus/${build}"
reference="${refmain}/bismark/${build}.fa"
refbismark="${refmain}/bismark"
refbiscuit="${refmain}/biscuit/${build}.fa"
outdir="${outfolder}/${sample}"
# prepare temporary folder
mkdir -p ${outdir}
#mkdir -p $outdir/QC_WGBS
#
cd $outdir
rm -rf *
############################################################## create symbolic links to read files in local outdir #####################################
ln -s $read1 ${sample}_R1.fastq.gz
ln -s $read2 ${sample}_R2.fastq.gz
### trim reads
trim_galore --paired --retain_unpaired --clip_R1 6 --clip_R2 6 --fastqc ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
############################################################## map reads to the reference genome #######################################################
rm -rf ${sample}_R1.fastq.gz
rm -rf ${sample}_R2.fastq.gz
mv ${sample}_R1_val_1.fq.gz ${sample}_R1.fastq.gz
mv ${sample}_R2_val_2.fq.gz ${sample}_R2.fastq.gz
gunzip ${sample}_R1.fastq.gz
gunzip ${sample}_R2.fastq.gz
#########################################################################################################################################################
bismark --bam --nucleotide_coverage --multicore 8 ${refbismark} -1 ${sample}_R1.fastq -2 ${sample}_R2.fastq
mv ${sample}_R1_bismark_bt2_pe.bam ${sample}_bismark_bt2_pe.bam
mv ${sample}_R1_bismark_bt2_pe.nucleotide_stats.txt ${sample}_bismark_bt2_pe.nucleotide_stats.txt
mv ${sample}_R1_bismark_bt2_PE_report.txt ${sample}_bismark_bt2_PE_report.txt 
bismark_methylation_extractor --no_overlap --gzip --bedGraph --no_header --comprehensive --multicore 5 --cytosine_report --ignore 3 --ignore_r2 3 --ignore_3prime 3 --ignore_3prime_r2 3 --buffer_size 40G --genome_folder ${refbismark} ${sample}_bismark_bt2_pe.bam
rm -rf ${sample}_R1.fastq
rm -rf ${sample}_R2.fastq
############################################################# QC component #############################################################################
samtools  sort -@ ${threads} -T ${outdir}/${sample} -o ${sample}_bismark_bt2_pe.coordsrt.bam ${sample}_bismark_bt2_pe.bam
samtools index -@ ${threads} ${sample}_bismark_bt2_pe.coordsrt.bam
rm -rf ${sample}_bismark_bt2_pe.bam
############################################################# operate methyldackel
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --CHG --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --CHH --minDepth ${mincov} -o ${sample}-bismethyl-${build}
#########################################
# generate TDF for visualization
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' ${sample}-bismethyl-${build}_CpG.bedGraph > ${sample}-${build}_viz.bedGraph
python ${BASEDIR}/stretch_bed.py ${sample}-${build}_viz.bedGraph ${sample}-${build}_CpG.igv $sample
igvtools toTDF ${sample}-${build}_CpG.igv ${sample}-${build}_CpG.tdf ${refbiscuit/.fa/.chrom.sizes}
rm -rf igv.log
# generate multiqc report
multiqc .

#!/bin/bash
#
# equivalent to bswalkermm -- current mouse release
#
set -e
set -x
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
refbiscuit="${refmain}/biscuit/${build}.fa"
refbismark="${refmain}/bismark"
outdir="${outfolder}/${sample}"
# prepare temporary folder
mkdir -p ${outdir}
#mkdir -p $outdir/QC_WGBS
#
cd $outdir
### create symbolic links to read files in local outdir
ln -s $read1 ${sample}_R1.fastq.gz
ln -s $read2 ${sample}_R2.fastq.gz
### trim reads
trim_galore --paired --gzip --retain_unpaired --fastqc ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
### map reads to the reference genome
rm -rf ${sample}_R1.fastq.gz
rm -rf ${sample}_R2.fastq.gz
mv ${sample}_R1_val_1.fq.gz ${sample}_R1.fastq.gz
mv ${sample}_R2_val_2.fq.gz ${sample}_R2.fastq.gz
biscuit align -t ${threads} ${refbiscuit} ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz | samtools sort -@ ${threads} -T ${outdir}/ -O bam -o ${sample}_biscuit.bam
samtools index ${sample}_biscuit.bam
samtools flagstat ${sample}_biscuit.bam > ${sample}_biscuit.flagstat
bismark_methylation_extractor --no_overlap --gzip --bedGraph --no_header --comprehensive --multicore 4 --cytosine_report --buffer_size 48G --genome_folder ${refbismark} ${sample}_biscuit.bam
cp ${sample}_biscuit.bismark.cov.gz ${sample}-${build}.cov.gz
gunzip ${sample}-${build}.cov.gz
############## QC component ####################
samtools  sort -@ ${threads} -T ${outdir}/${sample} -o ${sample}_bismark_bt2_pe.coordsrt.bam ${sample}_biscuit.bam
samtools index -@ ${threads} ${sample}_bismark_bt2_pe.coordsrt.bam
samtools flagstat ${sample}_bismark_bt2_pe.coordsrt.bam > ${sample}_biscuit.flagstat
samtools  sort -n -@ ${threads} -T ${outdir}/${sample} -o ${sample}_bismark_bt2_pe.querynam.bam ${sample}_bismark_bt2_pe.bam
samtools index -@ ${threads} ${sample}_bismark_bt2_pe.querynam.bam
####
### operate methyldackel
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --CHG --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam -@ ${threads} --CHH --minDepth ${mincov} -o ${sample}-bismethyl-${build}
#########################################
# generate TDF for visualization
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' ${sample}-bismethyl-${build}_CpG.bedGraph > ${sample}-${build}_viz.bedGraph
python ${BASEDIR}/stretch_bed.py ${sample}-${build}_viz.bedGraph ${sample}-${build}_CpG.igv $sample
igvtools toTDF ${sample}-${build}_CpG.igv ${sample}-${build}_CpG.tdf ${refbiscuit/.fa/.chrom.sizes}
rm -rf igv.log
####################################################################################################
unset DISPLAY
qualimap  bamqc -nt 8 -bam ${sample}_bismark_bt2_pe.coordsrt.bam -gd MOUSE -outdir ${outdir} -outformat HTML --java-mem-size=40G
############## extract DNA methylation as well as genetic information ##############################
biscuit pileup -o ${sample}_biscuit-${build}.vcf -q ${threads} ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam
bgzip ${sample}_biscuit-${build}.vcf
tabix -p vcf ${sample}_biscuit-${build}.vcf.gz
biscuit vcf2bed -e -k ${mincov} -t cg ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.cg.bed
biscuit vcf2bed -e -k ${mincov} -t snp ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.snp.bed
biscuit vcf2bed -e -k ${mincov} -t ch ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.ch.bed
tabix -h  ${sample}_biscuit-${build}.vcf.gz -R ${sample}-biscuit-${build}.snp.bed| grep -v lambda > ${sample}-biscuit-${build}.snp.vcf
bgzip -@ ${threads} ${sample}-biscuit-${build}.snp.vcf
tabix -f -p vcf ${sample}-biscuit-${build}.snp.vcf.gz
#### CGMap
cgmaptools convert bam2cgmap -b ${sample}_bismark_bt2_pe.coordsrt.bam -g ${reference} -o ${sample}
cgmaptools snv -i ${sample}.ATCGmap.gz -m bayes -v ${sample}.bayes.vcf -o ${sample}.bayes.snv --bayes-dynamicP
cgmaptools snv -i ${sample}.ATCGmap.gz -m binom -o ${sample}.binom.snv
gawk '{if(/^#/){print} else {print "chr"$0;}}' ${sample}.bayes.vcf > ${sample}.bayes2.vcf
cgmaptools asm -r ${refbiscuit} -b ${sample}_bismark_bt2_pe.coordsrt.bam -l ${sample}.bayes2.vcf > ${sample}.asm
cgmaptools mbin -i ${sample}.CGmap.gz -c 10 --CXY 5 -B 5000 -f png -t ${sample} -p ${sample} > ${sample}.5000bin.data
cgmaptools mstat -i ${sample}.CGmap.gz -c 10 -f png -t ${sample} -p ${sample} > ${sample}.stat.data
cgmaptools oac bin -i ${sample}.ATCGmap.gz -B 5000 -f png -p ${sample} -t ${sample} > ${sample}.oac_5000bin.data
cgmaptools mec bin -i ${sample}.CGmap.gz -B 5000 -f png -p ${sample} -t ${sample} > ${sample}.mec_5000bin.data
cgmaptools mec stat -i ${sample}.CGmap.gz -p ${sample} -f png > ${sample}.mec_stat.data
# viewbs

coverage2cytosine -CX -o ${sample}-${build}.tab --genome_folder ${refbismark}/ ${sample}-${build}.cov
mv ${sample}-${build}.tab.CX_report.txt ${sample}-${build}.tab
bgzip ${sample}-${build}.tab
tabix -p vcf ${sample}-${build}.tab.gz
ViewBS MethCoverage --reference ${refbiscuit} --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-MethCoverage
for ${chromosome} in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
do
ViewBS BisNonConvRate --minDepth 5 --height 20 --width 20 --chrom ${chromosome} --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir}/BisConvRate --prefix ${sample}-BisNonConvRate-${chromosome}
done
ViewBS GlobalMethLev --minDepth 5 --height 20 --width 20 --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-GlobalMethLev
ViewBS MethLevDist --minDepth 5 --height 20 --width 20 --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-MethLevDist --binMethLev 0.05
ViewBS MethGeno --height 20 --width 20 --genomeLength ${refbiscuit}.fai --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-MethGeno --context CG
ViewBS MethGeno --height 20 --width 20 --genomeLength ${refbiscuit}.fai --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-MethGeno --context CHG
ViewBS MethGeno --height 20 --width 20 --genomeLength ${refbiscuit}.fai --sample ${sample}-${build}.tab.gz,${sample} --outdir ${outdir} --prefix ${sample}-MethGeno --context CHH
# generate multiqc report
multiqc .

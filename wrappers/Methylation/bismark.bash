#!/bin/bash
#
### 
bismark --bam --nucleotide_coverage --multicore 8 ${refbismark} -1 ${sample}_R1.fastq -2 ${sample}_R2.fastq
mv ${sample}_R1_bismark_bt2_pe.bam ${sample}_bismark_bt2_pe.bam
mv ${sample}_R1_bismark_bt2_pe.nucleotide_stats.txt ${sample}_bismark_bt2_pe.nucleotide_stats.txt
mv ${sample}_R1_bismark_bt2_PE_report.txt ${sample}_bismark_bt2_PE_report.txt 
bismark_methylation_extractor --no_overlap --gzip --bedGraph --no_header --comprehensive --multicore 8 --cytosine_report --ignore 3 --ignore_r2 3 --ignore_3prime 3 --ignore_3prime_r2 3 --buffer_size 48G --genome_folder ${refbismark} ${sample}_bismark_bt2_pe.bam
rm -rf ${sample}_R1.fastq
rm -rf ${sample}_R2.fastq
############## QC component ####################
samtools  sort -@ ${threads} -T ${outdir}/${sample} -o ${sample}_bismark_bt2_pe.coordsrt.bam ${sample}_bismark_bt2_pe.bam
samtools index -@ ${threads} ${sample}_bismark_bt2_pe.coordsrt.bam
samtools flagstat ${sample}_bismark_bt2_pe.coordsrt.bam > ${sample}_biscuit.flagstat
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
qualimap  bamqc -nt 10 -bam ${sample}_bismark_bt2_pe.coordsrt.bam -gd HUMAN -outdir ${outdir} -outformat HTML --java-mem-size=48G
############## extract DNA methylation as well as genetic information ##############################
biscuit pileup -o ${sample}_biscuit-${build}.vcf -q ${threads} ${refbiscuit} ${sample}_bismark_bt2_pe.coordsrt.bam
bgzip ${sample}_biscuit-${build}.vcf
tabix -p vcf ${sample}_biscuit-${build}.vcf.gz
# extract CpG beta values from the VCF file
biscuit vcf2bed -e -k ${mincov} -t cg ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.cg.bed
biscuit vcf2bed -e -k ${mincov} -t snp ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.snp.bed
biscuit vcf2bed -e -k ${mincov} -t ch ${sample}_biscuit-${build}.vcf.gz | awk '{printf "%s\t%s\t%s\t%4.0f\t%4.0f\n", $1, $2, $2, $4*$5, $5}' > ${sample}-biscuit-${build}.ch.bed
#
tabix -h  ${sample}_biscuit-${build}.vcf.gz -R ${sample}-biscuit-${build}.snp.bed| grep -v lambda > ${sample}-biscuit-${build}.snp.vcf
bgzip -@ ${threads} ${sample}-biscuit-${build}.snp.vcf
tabix -f -p vcf ${sample}-biscuit-${build}.snp.vcf.gz
#
#### CGMap
cgmaptools convert bam2cgmap -b ${sample}_bismark_bt2_pe.coordsrt.bam -g ${reference} -o ${sample}
#cgmaptools snv -i ${sample}.ATCGmap.gz -m bayes -v ${sample}.bayes.vcf -o ${sample}.bayes.snv --bayes-dynamicP
#cgmaptools snv -i ${sample}.ATCGmap.gz -m binom -o ${sample}.binom.snv
#gawk '{if(/^#/){print}else{print "chr"$0;}}' ${sample}.bayes.vcf > ${sample}.bayes2.vcf
#cgmaptools asm -r ${reference} -b ${sample}_bismark_bt2_pe.coordsrt.bam -l ${sample}.bayes2.vcf > ${sample}.asm
cgmaptools mbin -i ${sample}.CGmap.gz -c 10 --CXY 5 -B 1000 -f png -t ${sample} -p ${sample} > ${sample}.1000bin.data
cgmaptools mstat -i ${sample}.CGmap.gz -c 5 -f png -t ${sample} -p ${sample} > ${sample}.stat.data
cgmaptools oac bin -i ${sample}.ATCGmap.gz -B 1000 -f png -p ${sample} -t ${sample} > ${sample}.oac_1000bin.data
cgmaptools mec bin -i ${sample}.CGmap.gz -B 1000 -f png -p ${sample} -t ${sample} > ${sample}.mec_1000bin.data
cgmaptools mec stat -i ${sample}.CGmap.gz -p ${sample} -f png > ${sample}.mec_stat.data

#!/bin/bash
#
biscuit align -t ${threads} ${refbiscuit} ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz | samtools sort -@ ${threads} -T ${outdir}/ -O bam -o ${sample}_biscuit.bam
samtools index ${sample}_biscuit.bam
samtools flagstat ${sample}_biscuit.bam > ${sample}_biscuit.flagstat
# extract DNA methylation as well as genetic information
biscuit pileup -o ${sample}_biscuit-${build}.vcf -q ${threads} ${refbiscuit} ${sample}_biscuit.bam
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
####
### operate methyldackel
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --CHG --minDepth ${mincov} -o ${sample}-bismethyl-${build}
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --CHH --minDepth ${mincov} -o ${sample}-bismethyl-${build}
#### merge contexts
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --mergeContext --minDepth ${mincov} -o ${sample}-bismethyl-merged-${build}
####
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --mergeContext --minDepth ${mincov} --fraction -o ${sample}-bismethyl-merged-${build}
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --mergeContext --minDepth ${mincov} --counts -o ${sample}-bismethyl-merged-${build}
MethylDackel extract ${refbiscuit} ${sample}_biscuit.bam -@ ${threads} --mergeContext --minDepth ${mincov} --logit -o ${sample}-bismethyl-merged-${build}
#####################################################################################################################################################
#########################################
# generate mbias plot
MethylDackel mbias -@ ${threads} ${refbiscuit} ${sample}_biscuit.bam ${sample}-biscuit-${build}
#########################################
# generate TDF for visualization
awk -F"\t" '{print $1"\t"$2"\t"$3"\t"$4}' ${sample}-bismethyl-${build}_CpG.bedGraph > ${sample}-${build}_viz.bedGraph
python ${BASEDIR}/stretch_bed.py ${sample}-${build}_viz.bedGraph ${sample}-${build}_CpG.igv $sample
igvtools toTDF ${sample}-${build}_CpG.igv ${sample}-${build}_CpG.tdf ${refbiscuit/.fa/.chrom.sizes}

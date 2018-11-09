#!/bin/sh
set -eu

## this is a wrapper for mapping single end and paired end data

STAR=`which STAR` # adapt to your needs
GENDIR=$codebase/$build/star

BASE=`basename $1`
DNAME=`dirname $1`
ZCAT="--readFilesCommand zcat"
THREADS=16
OPTS_P1="--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 100000000000 --genomeLoad LoadAndKeep --seedSearchStartLmax 8 --outFilterMultimapNmax 100 --outFilterMismatchNoverLmax 0.5 --outFileNamePrefix $DNAME/${BASE}"

if [ "$#" -ge 3  ]; then
   echo "illegal number of parameters"
   exit 1
fi

if [ "$#" -eq 2 ]; then
CMD="$STAR --runThreadN $THREADS --outBAMsortingThreadN 10 $OPTS_P1 $ZCAT  --genomeDir $GENDIR --readFilesIn $1 $2"
echo $CMD
$CMD
fi

if [ "$#" -eq 1 ]; then
CMD="$STAR --runThreadN $THREADS --outBAMsortingThreadN 10 $OPTS_P1 $ZCAT  --genomeDir $GENDIR --readFilesIn $1"
echo $CMD
$CMD
fi
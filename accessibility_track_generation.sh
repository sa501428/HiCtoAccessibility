#!/bin/bash

# dependencies
# samtools 1.13
# UCSC bedGraphToBigWig utility

#usage ./dnase_track_generation.sh [merged_dedup.bam] -c [chr] -t [threads] -g [chrom.sizes] -o [output file name]
# -c: list of chromsomes separated by "|"
# -t: number of threads to use
# -g: path to chrom.sizes file
# -o: output file name [bigWig format]
dedupBam=$1
bwPath="/gpfs0/work/suhas/scripts/p-e_paper_2021/accessibility"
threads=24
chr="chr1|chr10|chr11|chr12|chr13|chr14|chr15|chr16|chr17|chr18|chr19|chr2|chr20|chr21|chr22|chr3|chr4|chr5|chr6|chr7|chr8|chr9|chrX"
outFile="output.bw"

while getopts "g:t:c:ho:" opt; do
    case $opt in
        g) chromSizes=$OPTARG ;;
        h) printHelpAndExit 0;;
	c) chr=$OPTARG ;;
	t) threads=$OPTARG ;;
	o) outFile=$OPTARG ;;
        [?]) printHelpAndExit 1;;
    esac
done

## sort merged dedup bam by chr, retaining only true Hi-C contacts and excluding duplicates
## this is currently only compatible with merged_dedup.bam files from single end Juicer2 runs
samtools view -u -d "rt:0" -d "rt:1" -d "rt:2" -d "rt:3" -d "rt:4" -d "rt:5" -@ 24 -F 0x400 -q 1 $dedupBam |  samtools sort -@ 24 -m 6G -o reads.sorted.bam
samtools index -@ 24 reads.sorted.bam

## generate 1bp resolution pileup track in bedgraph format (temporary)
date
# if only a single chromosome is specified
if [[ $chr != *"|"* ]]; then
        samtools view -@ ${threads} -h reads.sorted.bam $chr | samtools sort -@ ${threads} -n -m 1G -O sam | awk -v chr=$chr '{for (i=12; i<=NF; i++) {if ($i ~ /^ip/) {split($i, ip, ":"); locus[ip[3]]++; break}}}END{for (i in locus) {              print chr "\t" i-1 "\t" i "\t" locus[i]}}' | sort -k1,1 -k2,2n --parallel=${threads} -S 6G > tmp.bedgraph

# if multiple chromosomes are specified, parallelize over chromosomes
else
        export SHELL=$(type -p bash)
        doit () {
                samtools view -@ 2 -h reads.sorted.bam $1 | samtools sort -n -m 1G -O sam | awk -v chr=$1 '{for (i=12; i<=NF; i++) {if ($i ~ /^ip/) {split($i, ip, ":"); locus[ip[3]]++; break}}}END{for (i in locus) {print chr "\t" i-                        1 "\t" i "\t" locus[i]}}' | sort -k1,1 -k2,2n -S 6G
        }

        export -f doit
        echo $chr | tr "|" "\n" | parallel -j $threads --will-cite --joblog temp.log -k doit | sort -k1,1 -k2,2n -S 6G > tmp.bedgraph
fi
date

## convert the temporary bedgraph output file to final bigWig format using the UCSC bedGraphToBigWig utility
$bwPath/bedGraphToBigWig tmp.bedgraph $chromSizes $outFile
rm tmp.bedgraph
rm reads.sorted.bam*
date

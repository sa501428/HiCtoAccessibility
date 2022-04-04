#!/bin/bash

mnd=$1
chrom_sizes=$2
threads=8

export SHELL=$(type -p bash)
export mnd=$mnd
doit () {
	awk -v chr=$1 'BEGIN{OFS="\t"}$2==chr{c[$3]++}$6==chr{c[$7]++}END{for (i in c) {print chr, i-1, i, c[i]}}' $mnd | sort -k2,2n
    }
export -f doit

awk '{print $1}' $chrom_sizes | parallel -j $threads --will-cite --joblog temp.log -k doit {} > tmp.bedgraph

sort -k1,1 -k2,2n -S6G tmp.bedgraph > merged30.bedgraph

bedGraphToBigWig merged30.bedgraph $chrom_sizes inter_30.bw
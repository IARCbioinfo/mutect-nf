#!/bin/bash
for v in `ls *.vcf | grep -v -E "$1|$2"`; do cat $v | mergeBed > "$v.bed"; done
for v in `ls *.vcf.gz | grep -v -E "$1|$2"`; do zcat $v | mergeBed > "$v.bed"; done
beds=''
for b in `ls *.bed`; do beds=$beds' '$b; done
bedops -m $beds > regionsI.bed
awk '{print $1"\t"($2)"\t"($3)}' regionsI.bed | sed -r 's/-[0-9]?[0-9]/1/g' | bedops -m /dev/stdin > regions.bed
#bgzip regions.bed
#tabix -p bed regions.bed.gz

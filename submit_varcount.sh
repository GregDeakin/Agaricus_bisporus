#!/bin/bash

#Assemble contigs using Bowtie
#$ -S /bin/bash
#$ -cwd
#$ -l virtual_free=8G

REF=$1
SAMPLE=$2

bowtie2 -p 8 --no-unal -x $REF -U $SAMPLE -S $REF.sam
samtools view -S -b $REF.sam >$REF.bam
samtools sort $REF.bam $REF.sorted
samtools mpileup -o $REF.pile.vcf -v -t DPR -u -f $REF.txt $REF.sorted.bam
bcftools call -Ov -v -m $REF.pile.vcf > $REF.vars.vcf

grep ^[^#] $REF.vars.vcf|awk -F"\t" '{print $1"\t"$10}'|awk -F":" '{print $1}'|awk -F"/" '{print $1"\t"$2}' > $REF.var.txt

viruses=$(awk -F"\t" '{print $1}' <$REF.var.txt|sort|uniq )
for v in $viruses
do
    x=$( grep -P "$v\t" $f.var.txt |awk -F"\t" '{if ($2!=$3) {print $1}}'|wc -l )
    echo $REF $v $x
done > $REF.var_count.txt

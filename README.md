# Agaricus_bisporus

##RNA-seq virus pipeline 
AB=~/projects/Ab_virome

###Trimming and joining with usearch (trimmomatic may be better)
```shell
counter=0
S=200
for f in $AB/MVX_Seq/Original_reads/*
do counter=$((counter+1))
    if (( $counter % 2 == 0 ))
    then
        R2=$f
        S=$((S+1))
        $METAGENOMICS/scripts/ujoin.sh $R1 $R2 ${S}.joined.fq $AB/MVX/joined 1 200
    fi
    R1=$f
done
```

###Align with bowtie to remove contaminants
```shell

bowtie2-build contaminants.fasta contaminants

S=200
for f in $AB/MVX/joined/*
  do S=$((S+1))
  $AB/scripts/bowtie.sh $f /home/deakig/projects/ab_virome/supp/contaminants  $AB/MVX/cleaned $S.cleaned
done
```

###Assemble with trinity
```shell
for f in $AB/MVX/cleaned/*
do
  S=$(echo $f|grep -ohP "[0-9]{2,}")
  $AB/scripts/assembly.sh $f trinity_$S
done
```

###Blast and filter for non-hits
Make blast db of known sequences which I'm not interested in (this is a second contaminants filter for things found within the RNA-seq data which I'm not interested in) 
```shell
makeblastdb -dbtype nucl -in combi -out combi
```

####Return no hit fasta file 
Finds all lines between the given fasta header and the next fasta header inclusive - it's a bit converluted as I haven't found a good way to remove the second header (and as slow as a slow thing).
```shell
blastn -db $AB/supp/combi -query $AB/MVX/assembled/201.fa -outfmt 7|grep -B2 "^# 0"|grep Query > X
while read s; do
  echo $s|awk -F" " '{print $3}'|xargs -I pattern sed -n "/pattern/,/>/p" $AB/MVX/assembled/201.fa|head -n -1 
done <X>nohits.fa
rm X
```

###Virus variant discovery
Finds virus snps

```
f=138
cl=201

bowtie2-build $AB/variants/$f/$f.txt $f
bowtie2 -p 8 --no-unal -x $f -U $AB/MVX2/cleaned/$cl.cleaned.fq -S $f.sam
samtools view -S -b $f.sam >$f.bam
samtools sort $f.bam $f.sorted
samtools mpileup -o $f.pile.vcf -v -t DPR -u -f $f.txt $f.sorted.bam
bcftools call -Ov -v -m $f.pile.vcf > $f.vars.vcf

grep ^[^#] $f.vars.vcf|awk -F"\t" '{print $1"\t"$10}'|awk -F":" '{print $1}'|awk -F"/" '{print $1"\t"$2}' > $f.var.txt

viruses=$(awk -F"\t" '{print $1}' <$f.var.txt|sort|uniq )
for v in $viruses
do
    x=$( grep -P "$v\t" $f.var.txt |awk -F"\t" '{if ($2!=$3) {print $1}}'|wc -l )
    echo $f $v $x
done > $f.var_count.txt
```

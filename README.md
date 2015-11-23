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
Finds all lines between the given fasta header and the next fasta header inclusive - it's a bit converluted as I haven't found a good way to remove the second header.
```shell
blastn -db $AB/supp/combi -query $AB/MVX/assembled/201.fa -outfmt 7|grep -B2 "^# 0"|grep Query > X
while read s; do
  echo $s|awk -F" " '{print $3}'|xargs -I pattern sed -n "/pattern/,/>/p" $AB/MVX/assembled/201.fa|head -n -1 
done <X>nohits.fa
rm X
```


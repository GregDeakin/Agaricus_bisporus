

# filter adapters, phix and rrna
for FR in $PROJECT_FOLDER/MVX3/fastq/*_R1.fq.gz; do
  RR=$(sed 's/_R1/_R2/' <<< $FR);   
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c MEGAFILT  \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/adapters/truseq.fa \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/phix_174.fa \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/ribokmers.fa.gz \
  $PROJECT_FOLDER/MVX3/filtered   $FR   $RR   true; 
done

# Agaricus bisporus gene removal (bbduk)
for FR in $PROJECT_FOLDER/MVX3/filtered/*_R1.fq.gz.filtered.fq.gz; do
  RR=$(sed 's/_R1/_R2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/contaminants/contaminants.fasta \
  $PROJECT_FOLDER/MVX3/cleaned \
  $FR \
  $RR \
  tossjunk \
  k=31 \
  t=2
done

# Known virus gene removal (bbduk)
for FR in $PROJECT_FOLDER/MVX3/cleaned/*_R1.fq.gz; do
  RR=$(sed 's/_R1/_R2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c filter -p bbduk \
  $PROJECT_FOLDER/MVX3/temp/all_vir.fa \
  $PROJECT_FOLDER/MVX3/cleaner \
  $FR \
  $RR \
  tossjunk \
  k=31 \
  t=2
done

# Human contaminant removal (BBMap)
for FR in $PROJECT_FOLDER/MVX3/cleaner/*_R1.fq.gz.filtered.fq.gz.filtered.fq.gz; do
  RR=$(sed 's/_R1/_R2/' <<< $FR)
  $PROJECT_FOLDER/metatranscriptomics_pipeline/scripts/PIPELINE.sh -c filter -p bbmap \
  $PROJECT_FOLDER/metatranscriptomics_pipeline/common/resources/contaminants/bbmap_human/hg19_main_mask_ribo_animal_allplant_allfungus.fa \
  $PROJECT_FOLDER/MVX3/cleanest \
  $FR \
  $RR \
  minid=0.95 \
  maxindel=3 \
  bwr=0.16 \
  bw=12 \
  quickmatch \
  fast \
  minhits=2 \
  t=2
done

for FR in $PROJECT_FOLDER/MVX3/cleaner/*_R1.fq.gz; do
  RR=$(sed 's/_R1/_R2/' <<< $FR)
  $PROJECT_FOLDER/metagenomics_pipeline/scripts/PIPELINE.sh -c normalise -p bbnorm \
  $PROJECT_FOLDER/MVX3/corrected \
  $FR \
  $RR  \
  t=2 \
  target=100 \
  min=2 \
  ecc=t \
  passes=1 \
  bits=16 prefilter
done

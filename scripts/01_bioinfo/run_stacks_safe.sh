#!/bin/bash

set -e

POPMAP="popmap-filogenia.txt"
IN_DIR="trimmed_data"
OUT_DIR="results_stacks"

mkdir -p $OUT_DIR

echo "--- Iniciando Pipeline: $(date) ---"

# Loop ustacks com trava individual
id=1
for sample in $(cut -f1 $POPMAP); do
    echo "Processando $id: $sample"
    ustacks -f $IN_DIR/${sample}.fastq.gz -o $OUT_DIR -i $id -p 6 -t gzfastq -d \
      -m 6 -M 4 -N 6 --name $sample --force-diff-len || { echo "Falha na amostra $sample"; exit 1; }
    let "id++"
done && \
cstacks -P $OUT_DIR -M $POPMAP -p 6 -n 3 && \
sstacks -P $OUT_DIR -M $POPMAP -p 6 && \
tsv2bam -P $OUT_DIR -M $POPMAP -t 6 && \
gstacks -P $OUT_DIR -M $POPMAP -t 6 && \

mkdir -p populations_R05 populations_R08

populations -P $OUT_DIR -M $POPMAP -O ./populations_R05 \
-t 32 -r 0.5 --min-mac 3 --max-obs-het 0.6 --structure --vcf

populations -P $OUT_DIR -M $POPMAP -O ./populations_R08 \
-t 32 -r 0.8 --min-mac 3 --max-obs-het 0.6 --structure --vcf

echo "--- Finalizado com sucesso: $(date) ---"

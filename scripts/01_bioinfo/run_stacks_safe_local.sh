#!/bin/bash

POPMAP="popmap-CB_CM_modVACB01_VACM33.txt"
IN_DIR="TEMP_CB_CM_CL"
OUT_DIR="results_stacks_mod_CB_CM"

echo "--- Iniciando Pipeline: $(date) ---"

# Loop ustacks com trava individual
id=1
for sample in $(cut -f1 $POPMAP); do
    echo "Processando $id: $sample"
    ustacks -f $IN_DIR/${sample}.fastq.gz -o $OUT_DIR -i $id -p 2 -t gzfastq -d -m 4 -M 2 -N 4 --name $sample --force-diff-len || { echo "Falha na amostra $sample"; exit 1; }
    let "id++"
done && \
cstacks -P $OUT_DIR -M $POPMAP -p 2 -n 2 && \
sstacks -P $OUT_DIR -M $POPMAP -p 2 && \
tsv2bam -P $OUT_DIR -M $POPMAP -t 2 && \
gstacks -P $OUT_DIR -M $POPMAP -t 2 && \
populations -P $OUT_DIR -M $POPMAP -p 4 -r 0.8 --vcf

echo "--- Finalizado com sucesso: $(date) ---"

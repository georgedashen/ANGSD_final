mkdir alignments  

for((i=35;i<=49;i++))  
do  
STAR --runMode alignReads \
     --runThreadN 8 \
     --genomeDir STARindex \
     --readFilesIn SRR87699${i}_1_val_1.fq SRR87699${i}_2_val_2.fq \
     --outFileNamePrefix alignments/SRR878699${i}. \
     --outSAMtype BAM SortedByCoordinate
done  
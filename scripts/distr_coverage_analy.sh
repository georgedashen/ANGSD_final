ln -s /athena/angsd/scratch/zhc4002/angsd/final  
spack load -r singularity@2.6.0  
RSEQC_IMAGE="/athena/angsd/scratch/simg/rseqc-3.0.1.simg"  
BED_FILE="./Homo_Sapiens.GRCh38.bed"  

for((i=35;i<=49;i++))  
do  
singularity exec $RSEQC_IMAGE read_distribution.py -r $BED_FILE -i ./alignments/SRR87699${i}.Aligned.sortedByCoord.out.bam > ./alignment_qc/SRR87699${i}.read_distribution.txt  

singularity exec $RSEQC_IMAGE geneBody_coverage.py -r $BED_FILE -i ./alignments/SRR87699${i}.Aligned.sortedByCoord.out.bam -o ./alignment_qc/SRR87699${i}
done
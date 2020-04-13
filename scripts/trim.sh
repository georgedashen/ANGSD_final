mkdir trim_fastq  

for((i=35;i<=49;i++))  
do  
trim_galore --paired --retained --stringency 10 -o /.trim_fastq ./raw_fastq/SRR87699${i}_1.fastq ./raw_fastq/87699${i}_2.fastq  
done  
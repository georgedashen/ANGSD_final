mkdir FastQC  
fastqc -o ./FastQC *fastq --extract  

cd FastQC  
mkdir multiqc  
multiqc ./ -o ./multiqc  
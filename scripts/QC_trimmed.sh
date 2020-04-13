mkdir FastQC_trimmed  
fastqc -o ../FastQC_trimmed $(ls | grep val) --extract  

cd ../FastQC_trimmed  
mkdir multiqc  
multiqc ./ -o ./multiqc  
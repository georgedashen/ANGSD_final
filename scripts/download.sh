prefetch $(cat SRR_Acc_List.txt)  

for(i=35;i<=49;i++)  
do  
fasterq-dump --split-3 SRR87699${i}  
rm SRR87699${i}.sra  
done  
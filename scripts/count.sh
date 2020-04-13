spack load subread

featureCounts -p -B -T 8 -t exon -g gene_id -a Homo_sapiens.GRCh38.gtf -o counts.txt ./alignments/*bam
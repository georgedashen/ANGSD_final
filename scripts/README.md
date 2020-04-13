Scripts used in this project for raw RNA-seq data processing.

QC_original.sh: FastQC and MultiQC for raw RNA-seq data

QC_trimmed.sh: FastQC and MultiQC for trimmed reads data

STAR_align.sh: implement STAR for alignment, index building not included

bedGenerate.sh: generate BED file from gtf annotation

count.sh: featureCounts

distr_coverage_analy.sh: RSeQC for alignment QC

download.sh: download all RNA-seq data from a list of accession number of all samples by prefetch, then fastq-dump to split them

trim.sh: use trim_galore to trim low quality bases and potential adapters

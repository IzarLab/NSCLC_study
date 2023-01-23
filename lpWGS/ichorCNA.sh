# Author: Somnath Tagore, Ph.D. Title: Low-pass WGS analysis
# Script Name: ichorCNA.sh
# Last Updated: 02/19/2022

# fastq to sam
/home/ubuntu/anaconda3/bin/bowtie2 -x /home/ubuntu/wgs/hg38index  -1 /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001_R2_001.fastq.gz -2 /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001_R1_001.fastq.gz -S /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sam

# sam to unsorted bam
/home/ubuntu/anaconda3/bin/samtools  view -bS /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sam -o /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.bam

# unsorted sam to sorted bam
/home/ubuntu/anaconda3/bin/samtools sort /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.bam -o /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sorted.bam
/home/ubuntu/anaconda3/bin/samtools index /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sorted.bam

# bam to wig
/home/ubuntu/anaconda3/bin/readCounter --window 1000000 --quality 20 \
--chromosome "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY" \
/home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sorted.bam > /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sorted.wig

# CNA analysis using ichorCNA
sudo Rscript /home/ubuntu/wgs/ichorCNA/scripts/runIchorCNA.R --id KRAS4_S113_L001_tumor \
  --WIG /home/ubuntu/wgs/KRAS4_S113_L001/KRAS4_S113_L001.sorted.wig --ploidy "c(2,3)" --normal "c(0.5,0.6,0.7,0.8,0.9)" --maxCN 5 \
  --gcWig /home/ubuntu/wgs/ichorCNA/inst/extdata/gc_hg38_1000kb.wig \
  --mapWig /home/ubuntu/wgs/ichorCNA/inst/extdata/map_hg38_1000kb.wig \
  --centromere /home/ubuntu/wgs/ichorCNA/inst/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt \
  --normalPanel /home/ubuntu/wgs/ichorCNA/inst/extdata/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds \
  --includeHOMD False --chrs "c(1:22, \"X\")" --chrTrain "c(1:22)" \
  --estimateNormal True --estimatePloidy True --estimateScPrevalence True \
  --scStates "c(1,3)" --txnE 0.9999 --txnStrength 10000 --outDir /home/ubuntu/wgs/KRAS4_S113_L001


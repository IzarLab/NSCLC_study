conda create -n luca-cnv
conda activate luca-cnv
mamba install -n luca-cnv -c conda-force jupyter numpy scipy pandas pip seaborn plotly tqdm
conda env export > ./luca-cnv.env.yml


# conda create -n luca-cnv-numbat
# mamba install -n luca-cnv-numbat -c bioconda -c conda-forge \
# -c bioconda -c conda-forge \
#     cellsnp-lite samtools  r-essentials r-base 'r-seurat<5.0.0' \
#     r-irkernel Jupyter r-devtools \
#     bioconductor-ggtree bioconductor-genomicranges  r-gifski r-ggpubr \
#     bioconductor-milor r-beeswarm bioconductor-destiny bioconductor-mast \
#     bioconductor-clusterprofiler bioconductor-org.hs.eg.db bioconductor-ucell \
#     bioconductor-deseq2 r-matrix.utils bioconductor-fgsea
# conda activate numbat && conda env export > luca-cnv-numbat.env.yml && conda deactivate
mamba env create -f numbat.env.yml
sudo apt update && sudo apt install unzip moreutils libgomp1 libcairo2-dev

# Fix errors in irlba
R -e 'install.packages("Matrix", type = "source")'
R -e 'install.packages("irlba", type = "source")'

NUMBAT="/home/ubuntu/InstallTemp/numbat"  # Local dir for numbat requirements
mkdir -p "${NUMBAT}"
wget -O "${NUMBAT}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz" \
    https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
wget -O "${NUMBAT}/1000G_hg38.zip" \
    http://pklab.med.harvard.edu/teng/data/1000G_hg38.zip
wget -O "${NUMBAT}/Eagle_v2.4.1.tar.gz" \
    https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/Eagle_v2.4.1.tar.gz
tar -xvzf "${NUMBAT}/Eagle_v2.4.1.tar.gz"
gunzip "${NUMBAT}/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
unzip "${NUMBAT}/1000G_hg38.zip"
rm "${NUMBAT}/Eagle_v2.4.1.tar.gz"; rm "${NUMBAT}/1000G_hg38.zip"
ln -s "${NUMBAT}/Eagle_v2.4.1/eagle" "${CONDA_PREFIX}/bin/eagle"  # IMPORTANT!

conda activate luca-cnv-numbat
Rscript -e 'install.packages("numbat", dependencies=TRUE, repos="http://cran.us.r-project.org")'
Rscript -e 'install.packages("devtools")'

# Rscript "${CONDA_PREFIX}/lib/R/library/numbat/bin/pileup_and_phase.R" to run preprocessing script

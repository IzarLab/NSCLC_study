#!/bin/sh
#SBATCH --job-name="process_dataset_brainmets"
#SBATCH --partition=single
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem=16gb

module load devel/java_jdk/1.18
module load system/singularity/3.11.3

nextflow -c /gpfs/bwfor/work/ws/hd_hl269-brainmets/nextflow_configs/config_exploratory.yml \
    run kbestak/nf_mcmicro -r unstitch_restitch \
    --input_cycle /gpfs/bwfor/work/ws/hd_hl269-brainmets/uint16_data/dataset/samplesheet_dataset.csv \
    -profile singularity \
    --marker_sheet /gpfs/bwfor/work/ws/hd_hl269-brainmets/uint16_data/dataset/marker_sheet.csv \
    --skip_registration \
    --extract_channel \
    --unstitch_restitch \
    --pixel_size 0.4977523 \
    --mesmer_model "whole-cell" \
    --nuclear_channel 0 \
    --outdir /gpfs/bwfor/work/ws/hd_hl269-brainmets/uint16_data/dataset \
    --spotiflow_prob_threshold 0.6 \
    --spot_intensity_threshold 150 \
    -resume

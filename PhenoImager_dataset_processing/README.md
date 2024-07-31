# NSCLC Brainmets
This repository contains code snippets and scripts used for processing the matched primary and metastasis dataset.

## Image processing

The whole dataset was processed with a branch of nf-core/mcmicro, an adaptation of the nf-core (Ewels et al. 2020) expansion of MCMICRO (Schapiro et al. 2022) that can be found here: https://github.com/kbestak/nf_mcmicro/tree/unstitch_restitch

Image processing was performed on the bwForCluster Helix with the ./image_processing/process_dataset.sh script, ./image_processing/config_exploratory.yml config, and ./image_processing/samplesheet_dataset.csv samplesheet file.

The pipeline consists of the following steps:

* Image unstitching into 4 "tiles" (as Mesmer was breaking on these images)
* Segmentation with Mesmer
* Segmentation mask restitching
    * label assignment based on centroids
    * filtering by area
* MCQUANT - creates cell by marker matrix (mean intensities) and region props from masks
* Channel extraction - the cGAS channel is extracted from the original multichannel image
* Spot detection with Spotiflow (probability threshold 0.5)
* Spot2cell - spot assignment to cells
* combine spot2cell and mcquant outputs

## Cell type calling

Cell type calling was performed on the bwForCluster Helix with the ./cell_type_calling/cell_type_calling.ipynb Jupyter notebook from output files of the above pipeline using Scimap.

Detailed explanation can be found in the Jupyter notebook.

## Plotting and statistical analysis

Statistical analysis was done in R and RStudio 


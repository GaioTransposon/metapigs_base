#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=120:00:00
#PBS -l mem=70g
#PBS -N run_guppy.sh
#PBS -M daniela.gaio@uts.edu.au

cd /shared/homes/s1/pig_microbiome/phy_10M_202109/PS_temp

# compute per-sample alpha diversity with various diversity metrics
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fpd plate_*.fastq.gz/* > all.alphadiv

# cluster the samples: performs squash clustering
#/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy squash plate_*.fastq.gz/*jplace

# edge PCA to explore variation in community composition among samples: performs edge principal components
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_ plate_*.fastq.gz/*jplace

# run guppy fat: makes trees with edges fattened in proportion to the number of reads
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_ plate_*.fastq.gz/*jplace
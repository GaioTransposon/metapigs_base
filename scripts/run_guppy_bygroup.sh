#!/bin/bash
#PBS -l ncpus=20
#PBS -l walltime=120:00:00
#PBS -l mem=100g
#PBS -N run_guppy_bygroup.sh
#PBS -M daniela.gaio@uts.edu.au


cd /shared/homes/s1/pig_microbiome/phy_10M_202109/PS_temp


# edge PCA to explore variation in community composition among samples - by group
for f in /shared/homes/s1/pig_microbiome/phy_10M_202109/guppy_groups/*.txt
do N=$(basename $f)
#cat ../guppy_groups/$N
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_$N `cat ../guppy_groups/$N`
done


# run guppy fat - by group
for f in /shared/homes/s1/pig_microbiome/phy_10M_202109/guppy_groups/*.txt
do N=$(basename $f)
#cat ../guppy_groups/$N
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_$N `cat ../guppy_groups/$N`
done




#!/bin/bash
#PBS -l ncpus=10
#PBS -l walltime=120:00:00
#PBS -l mem=50g
#PBS -N run_guppy.sh
#PBS -M daniela.gaio@student.uts.edu.au


cd /shared/homes/s1/pig_microbiome/phy_10M/PS_temp



# edge PCA to explore variation in community composition among samples - by group
for f in /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/*.txt # or: /shared/homes/12705859/phylosift_metapigs_20200225/*.txt
do N=$(basename $f)
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy epca --prefix pca_$N `cat $N`
done


# run guppy fat - by group
for f in /shared/homes/s1/pig_microbiome/phy_10M/guppy_groups/*.txt # or: /shared/homes/12705859/phylosift_metapigs_20200225/*.txt
do N=$(basename $f)
/shared/homes/s1/pig_microbiome/phylosift_v1.0.1/bin/guppy fat --prefix fat_$N `cat $N`
done




#!/bin/bash
#PBS -l ncpus=4
#PBS -l walltime=48:00:00
#PBS -l mem=24g

export OMP_NUM_THREADS=4
cd $PBS_O_WORKDIR


# run guppy epca by group
for f in /shared/homes/12705859/phylosift_metapigs/*.txt # or: /shared/homes/12705859/phylosift_metapigs_20200225/*.txt
do N=$(basename $f)
/shared/homes/12705859/phylosift_v1.0.1/bin/guppy epca --prefix $N `cat $N`
done


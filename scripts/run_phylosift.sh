
#!/bin/bash
#PBS -l ncpus=28
#PBS -l walltime=120:00:00
#PBS -l mem=150g
#PBS -N run_phylosift.sh
#PBS -M daniela.gaio@student.uts.edu.au


cd /shared/homes/12705859/phylosift_202105

# run phylosift
while read  first  second
do
/shared/homes/12705859/phylosift_v1.0.1/bin/phylosift all \
      --threads 24 --disable_updates \
      --chunks 1 --chunk_size 100000000 \
      --paired $first $second \
      --output out
done < "all_reads.txt"




# compute per-sample alpha diversity with various metrics
#~/phylosift_v1.0.1/bin/guppy fpd *.jplace > all.alphadiv

# cluster the samples 
#~/phylosift_v1.0.1/bin/guppy squash *.jplace > all.clust

# edge PCA to explore variation in community composition among samples
#~/phylosift_v1.0.1/bin/guppy epca --prefix pca *.jplace


# run guppy fat
#for f in /shared/homes/12705859/phylosift_metapigs_20200225/plate_*.jplace
#do
#~/phylosift_v1.0.1/bin/guppy fat --prefix fat *.jplace
#done
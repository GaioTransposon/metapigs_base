


# SortMeRNA - extract 16S using bact. database


export PATH=/shared/homes/12705859/sortmerna/sortmerna-4.0.0-Linux/bin/:$PATH

for i in `ls /shared/homes/s1/pig_microbiome/sortmerna_16S/plate_*.gz`; do N=$(basename $i); sortmerna --ref /shared/homes/12705859/sortmerna-2.1b/rRNA_databases/silva-bac-16s-id90.fasta --reads $i -fastx -blast 1 --num_alignments 1 --workdir /shared/homes/12705859/sortmerna_aligned/$N; done
done



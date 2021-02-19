#script name: sortmerna_filter.sh


#!/bin/bash
#PBS -q i3q
#PBS -l ncpus=56
#PBS -l walltime=120:00:00
#PBS -l mem=200g
#PBS -N sortmerna_filter.sh
#PBS -M daniela.gaio@student.uts.edu.au


# STEPS: 
#0. add sample name to fastq files before "@"
#1. extract from .fq files all headers and create a file "fast_header"
#2. join fastq_header to blast output as last column
#3. filter blast based on : e-value, %id, align length
#4. extract the last column to make a list of reads to extract
#5. use list to filter the fastq files based on the list

#blast files here : 
#/shared/homes/12705859/sortmerna_aligned/all_aligned_blast

#fq files here: 
#/shared/homes/12705859/sortmerna_aligned/all_aligned_fq


############

#Step 0: 
# copy all original fq files into new dir: 
mkdir ~/sortmerna_aligned/all_aligned_fq/original
cp plate*.fq original/.

# insert filename in header
awk -i inplace '/@/{sub("@","&"FILENAME"_");sub(/\.fq/,x)}1' plate*.fq
# headers of fastq file now look like this: 
#@plate_3_E11_S277.rrna.gz_aligned.fq_A00152:61:HF7TGDSXX:1:1101:17409:14700 1:N:0:TACTTACC+TTCGTGGC

############

#Step 1: 
cd /shared/homes/12705859/sortmerna_aligned/all_aligned_fq
for f in `ls /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/plate*.fq`
do
filename=$(basename $f)
N="${filename%.*}"
cat $f | grep "@" > header_$N
done

############

#Step 2: 
cd /shared/homes/12705859/sortmerna_aligned/all_aligned_fq
for f in `ls /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/plate*.fq`
do
filename=$(basename $f)
N="${filename%.*}"
echo ../all_aligned_blast/$N.blast header_$N
paste -d' ' ../all_aligned_blast/$N.blast header_$N > ../pasted/$N
done

############

#Step 3: 
cd /shared/homes/12705859/sortmerna_aligned/pasted
mkdir filtered
for p in `ls /shared/homes/12705859/sortmerna_aligned/pasted/*`
do
filename=$(basename $p)
P="${filename%.*}"
awk -v x=1e-30 '$11 <= x'  $p | awk -v x=80 '$3 >= x'  | awk -v x=100 '$4 >= x'  > filtered/$P
done

############

#Step 4: 
cd /shared/homes/12705859/sortmerna_aligned/pasted/filtered
cat plate* | cut -f 12 | awk -F@ '{ print $NF }' >> all_filtered_reads_headers.txt
#cat all_filtered_reads_headers.txt | wc -l
#32419310

############


#Step 5: 

# test: 
#conda activate py_3.5
#filterbyname.sh in=/shared/homes/12705859/sortmerna_aligned/all_aligned_fq/plate_3_G7_S247.rrna.fq.gz_aligned.fq names=/shared/homes/12705859/sortmerna_aligned/pasted/filtered/#all_filtered_reads_headers.txt out=test_filt.fq include=true 
#Input is being processed as unpaired
#Time:                         	111.738 seconds.
#Reads Processed:       91376 	0.82k reads/sec
#Bases Processed:      13797k 	0.12m bases/sec
#Reads Out:          47940
#Bases Out:          7238940
# it works!


# it was run as nano run_filterbyname.sh (129323.hpcnode0) (it does use all the cpus available)

source activate py_3.5

cd /shared/homes/12705859/sortmerna_aligned/all_aligned_fq
rm -r filtered_fq
mkdir filtered_fq

for fq in `ls /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/plate*.fq`
do
filename=$(basename $fq)
N="${filename%.*}"
filterbyname.sh in=$fq names=/shared/homes/12705859/sortmerna_aligned/pasted/filtered/all_filtered_reads_headers.txt out=filtered_fq/$N include=true 
done


# check if the number of reads extracted from the original match with the number of reads from the list: 
#cat /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/plate_3_G7_S247.rrna.fq.gz_aligned.fq | grep "@" | wc -l
#91376
#cat /shared/homes/12705859/sortmerna_aligned/pasted/filtered/all_filtered_reads_headers.txt | grep "plate_3_G7" | wc -l
#47940
#cat /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/filtered_fq/plate_3_G7_S247.rrna.fq.gz_aligned | grep "@" | wc -l
#47940



#Step 6: 
#cd /shared/homes/12705859/sortmerna_aligned/all_aligned_blast
#rm -r aligned_blast_with_header
#mkdir aligned_blast_with_header
#for blast in `ls plate*.blast`
#do
#filename=$(basename $blast)
#B="${filename%.*}"
#awk '{print FILENAME"\t"$0}' $blast > aligned_blast_with_header/$B
#done

#cd aligned_blast_with_header
#rm all_blast.tsv
#cat plate* >> all_blast.tsv
# keep only cols of interest (to reduce file size)
#cat all_blast.tsv | cut -f 1,3,4,5,12 > all_blast_essential.tsv
# add header
#echo -e "Sample\tTargetID\tSeqIdentity\tAlignLength\tEvalue" | cat - all_blast_essential.tsv > all_blast_essential_wHeader.tsv








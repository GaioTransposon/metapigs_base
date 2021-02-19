#script name: RDP_Krona.sh

#!/bin/bash
#PBS -q i3q
#PBS -l ncpus=56
#PBS -l walltime=120:00:00
#PBS -l mem=200g
#PBS -N RDP_Krona
#PBS -M daniela.gaio@student.uts.edu.au


# STEPS: 
#1. feed filtered fq files to RDP classifier
#2. split classification output into moms vs piglets using lists 
#3. run Krona

############

#Step 1:
cd /shared/homes/12705859/sortmerna_aligned/all_aligned_fq/filtered_fq
# run RDP classify
for s in plate*aligned; do java -Xmx1g -jar /shared/homes/12705859/RDP_classifier/rdp_classifier_2.13/dist/classifier.jar     classify     -g 16srrna      -b `basename ${s}`.bootstrap      -h `basename ${s}`.hier.tsv      -o `basename ${s}`.class.tsv      ${s}; done

########################################################################################################################for s in plate_7_G*aligned; do java -Xmx1g -jar /shared/homes/12705859/RDP_classifier/rdp_classifier_2.13/dist/classifier.jar     classify     -g 16srrna      -b `basename ${s}`.bootstrap      -h `basename ${s}`.hier.tsv      -o `basename ${s}`.class.tsv      ${s}; done

############

#Step 2: 
# AIM: merge files based on (mothers or piggies) list, so the merged RDP classification files can be fed to Krona: 
# expected in folder: moms_list.txt + piggies_list.txt from Github repo metapigs_base/middle_dir 

# move class files to own folder
mkdir classified_out
mv plate*.class.tsv classified_out/.
cd classified_out

# rename class files to include only sample name 
for f in plate_*.tsv; do
   mv -- "$f" "${f/_S*./.}"
done
# remove the .tsv extension
rename -- .tsv '' *.tsv

# concatenate class files based on lists:
{ xargs cat < moms_list.txt ; } > classifications_moms.tsv
{ xargs cat < piggies_list.txt ; } > classifications_piggies.tsv

# create a concatenated class file of all samples: 
cat plate* > classifications_all.tsv

##############

#Step 3: 
singularity exec docker://quay.io/biocontainers/krona:2.7.1--pl526_1 ktImportRDP -o moms_krona.html -c classifications_moms.tsv
singularity exec docker://quay.io/biocontainers/krona:2.7.1--pl526_1 ktImportRDP -o piggies_krona.html -c classifications_piggies.tsv
singularity exec docker://quay.io/biocontainers/krona:2.7.1--pl526_1 ktImportRDP -o all_krona.html -c classifications_all.tsv

# now move html files to local and view it! 

# use classifications_all.tsv as input for RDP_analyze.R


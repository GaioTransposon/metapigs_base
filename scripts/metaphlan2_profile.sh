
#!/bin/bash
#PBS -q workq
#PBS -l select=1:ncpus=8:mem=32gb

source .bashrc
source activate blah

for f in ~/mock_communities/*.fastq.gz
do
  metaphlan2.py $f --input_type fastq --nproc 4 > ${f%.fastq.gz}_profile.txt
done
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
# generate species only merged table: 
grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
# visualise:
~/miniconda3/envs/blah/bin/hclust2.py -i mock_communities_merged_abundance_table_species.txt -o abundance_heatmap_species.png --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300

for f in ~/protexin/*.fastq.gz
do
  metaphlan2.py $f --input_type fastq --nproc 4 > ${f%.fastq.gz}_profile.txt
done
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
# generate species only merged table: 
grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
# visualise:
~/miniconda3/envs/blah/bin/hclust2.py -i mock_communities_merged_abundance_table_species.txt -o abundance_heatmap_species.png --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300

for f in ~/coli_guard/*.fastq.gz
do
  metaphlan2.py $f --input_type fastq --nproc 4 > ${f%.fastq.gz}_profile.txt
done
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
# generate species only merged table: 
grep -E "(s__)|(^ID)" merged_abundance_table.txt | grep -v "t__" | sed 's/^.*s__//g' > merged_abundance_table_species.txt
# visualise:
~/miniconda3/envs/blah/bin/hclust2.py -i mock_communities_merged_abundance_table_species.txt -o abundance_heatmap_species.png --ftop 25 --f_dist_f braycurtis --s_dist_f braycurtis --cell_aspect_ratio 0.5 -l --flabel_size 6 --slabel_size 6 --max_flabel_len 100 --max_slabel_len 100 --minv 0.1 --dpi 300










last_dir=$(basename $(dirname $PWD))
current_dir=$(basename $PWD)
cp -a /data/project/${last_dir}/${current_dir}/otu_table.biom .

export PATH=/data/miniconda3/envs/qiime1/bin/:$PATH

filter_otus_from_otu_table.py --min_count_fraction 0.001 -i  otu_table.biom  -o otu_table.new.biom

biom convert -i otu_table.new.biom -o otu_table.txt --to-tsv --header-key taxonomy

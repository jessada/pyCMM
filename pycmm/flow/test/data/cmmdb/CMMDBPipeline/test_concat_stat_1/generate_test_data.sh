#!/bin/bash
source $PYCMM/bash/cmm_functions.sh


input_dir="/proj/b2011117/nobackup/private/projects_output/cmmdb/mutstat/SweGen/first_release/data_out"
input_prefix="$input_dir/SweGen"

output_prefix="test_concat_stat_1"

# generate data for chromosome 1
input_stat="$input_prefix"_1.stat
output_stat="$output_prefix"_1.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome 2
input_stat="$input_prefix"_2.stat
output_stat="$output_prefix"_2.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome 7
input_stat="$input_prefix"_7.stat
output_stat="$output_prefix"_7.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome 12
input_stat="$input_prefix"_12.stat
output_stat="$output_prefix"_12.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome 20
input_stat="$input_prefix"_20.stat
output_stat="$output_prefix"_20.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome X
input_stat="$input_prefix"_X.stat
output_stat="$output_prefix"_X.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome MT
input_stat="$input_prefix"_MT.stat
output_stat="$output_prefix"_MT.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"

# generate data for chromosome X
input_stat="$input_prefix"_Y.stat
output_stat="$output_prefix"_Y.stat
head -1 "$input_stat" > $output_stat
sed -n 14700,14705p  "$input_stat" >> "$output_stat"
sed -n 45,50p  "$input_stat" >> "$output_stat"
sed -n 167600,167650p  "$input_stat" >> "$output_stat"


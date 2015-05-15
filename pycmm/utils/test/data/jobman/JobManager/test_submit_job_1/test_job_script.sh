#!/bin/bash 

start_pos=$1
end_pos=$2
out_file=$3
for (( i=$start_pos; i<=$end_pos; i++ ))
do
   echo "Welcome $i times" >> $out_file
done


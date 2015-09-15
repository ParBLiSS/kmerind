#!/bin/bash

#file
filename=$1;
echo $filename;

## store everything as arrays/

#get the container type
containers=( `grep "RANK" $filename | cut -d , -f 2` );

#get the index type
indices=( `grep "RANK" $filename | cut -d , -f 3 | cut -d ' ' -f 2` );

#get the proc count
procs=(`grep "test header" $filename | cut -d / -f 2 | cut -d ' ' -f 1`)

#get min
build_dur_mins=(`grep -P "test\tdur_min\t" $filename | cut -d [ -f 2 | cut -d , -f 1`)
#get max
build_dur_maxs=(`grep -P "test\tdur_max\t" $filename | cut -d [ -f 2 | cut -d , -f 1`)
#get mean
build_dur_means=(`grep -P "test\tdur_mean\t" $filename | cut -d [ -f 2 | cut -d , -f 1`)
#get stdev
build_dur_stdevs=(`grep -P "test\tdur_stdev\t" $filename | cut -d [ -f 2 | cut -d , -f 1`)

#get min
count_dur_mins=(`grep -P "test\tdur_min\t" $filename | cut -d [ -f 2 | cut -d , -f 6`)
#get max
count_dur_maxs=(`grep -P "test\tdur_max\t" $filename | cut -d [ -f 2 | cut -d , -f 6`)
#get mean
count_dur_means=(`grep -P "test\tdur_mean\t" $filename | cut -d [ -f 2 | cut -d , -f 6`)
#get stdev
count_dur_stdevs=(`grep -P "test\tdur_stdev\t" $filename | cut -d [ -f 2 | cut -d , -f 6`)

#get min
query_dur_mins=(`grep -P "test\tdur_min\t" $filename | cut -d [ -f 2 | cut -d , -f 7`)
#get max
query_dur_maxs=(`grep -P "test\tdur_max\t" $filename | cut -d [ -f 2 | cut -d , -f 7`)
#get mean
query_dur_means=(`grep -P "test\tdur_mean\t" $filename | cut -d [ -f 2 | cut -d , -f 7`)
#get stdev
query_dur_stdevs=(`grep -P "test\tdur_stdev\t" $filename | cut -d [ -f 2 | cut -d , -f 7`)


for ((i=0;i<${#containers[@]};++i)); do
  echo "${containers[$i]}_${indices[$i]},${procs[$i]},${build_dur_mins[$i]},${build_dur_maxs[$i]},${build_dur_means[$i]},${build_dur_stdevs[$i]},${count_dur_mins[$i]},${count_dur_maxs[$i]},${count_dur_means[$i]},${count_dur_stdevs[$i]},${query_dur_mins[$i]},${query_dur_maxs[$i]},${query_dur_means[$i]},${query_dur_stdevs[$i]}"

done;

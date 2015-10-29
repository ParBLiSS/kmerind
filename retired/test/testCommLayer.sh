#!/bin/sh

set rootdir=/home/tpan/build/bliss

for iter in `seq 1 100`
do
	for nt in `seq 1 4`
	do

		for count in 1 10 100 1000 10000
		do
			echo "RUN: bin/testCommLayer $nt $count"	
			${rootdir}/bin/testCommLayer $nt $count 
		
			for np in `seq 2 4`
			do
				echo "RUN: mpirun -np $np bin/testCommLayer $nt $count"	
				mpirun -np $np ${rootdir}/bin/testCommLayer $nt $count
# 2>&1 comm-layer-run.${np}.${nt}.${count}.${iter}.log
		
			done
		done
	done
done

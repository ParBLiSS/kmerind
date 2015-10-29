#!/bin/sh

for iter in `seq 1 100`
do

		for count in 1 10 100 1000 10000
		do
			bin/testDistrMap $count &> distr-map-run.1.${count}.${iter}.log
		
			for np in `seq 1 4`
			do
	
				mpirun -np $np bin/testDistrMap $count &> distr-map-run.${np}.${count}.${iter}.log
		
			done
		done
done

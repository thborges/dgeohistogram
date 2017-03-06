#!/bin/bash

for i in `seq 10 $max`
do
	echo "Execution $i"
	./main areafs fix 100 100 2 datasets/hidrografia.shp 0.20
done


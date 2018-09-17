#!/bin/bash

#remove if exist result folder
#rm -rf "/home/isabelladfn/experiment"

#create main folder to write results
#mkdir "/home/isabelladfn/experiment"

#array with lines from input folder, the commands are here
mapfile -t myArray < input.txt

#read file and write in array each line
while IFS= read -r line; do echo "$line"; done < input.txt

#for (( i = 0 ; i < 32 ; i++))
#do
#  echo "Element [$i]: ${myArray[$i]}"
#done

#start experiment
#--------------------------
for (( i=0; i<24; i++ ))
do
	mkdir /home/isabelladfn/experiment/"$((i+40))"
	echo "Size,ARE,STD,SUM,Method,Split.Qnt,Name" >> /home/isabelladfn/dgeohistogram/data.csv
	for (( j=0; j<10; j++ ))
	do
		echo "------------------"
		echo "Execution $((i+40))-$j"
		echo "------------------"
		eval "${myArray[$i]}"
		mkdir "/home/isabelladfn/experiment/$((i+40))/$((i+40))-$j"
		cp -r /home/isabelladfn/dgeohistogram/histogram /home/isabelladfn/experiment/"$((i+40))"/"$((i+40))-$j"/histogram
		cp -r /home/isabelladfn/dgeohistogram/mbrs /home/isabelladfn/experiment/"$((i+40))"/"$((i+40))-$j"/mbrs
		make clean_files
		echo "------------------"
		echo
		if [ "$j" == 9 ]; then
			cp /home/isabelladfn/dgeohistogram/data.csv /home/isabelladfn/experiment/"$((i+40))"/data.csv
			rm /home/isabelladfn/dgeohistogram/data.csv
		fi
	done
done
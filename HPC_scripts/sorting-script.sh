#!/bin/bash

#Requires plate-date to be entered as first argument

for rowLetter in B C D E F G H I J K L M N O
do
	mkdir ./output-data/$1_${rowLetter}_row
	mv ./output-data/$1_${rowLetter}*.csv ./output-data/$1_${rowLetter}_row
done


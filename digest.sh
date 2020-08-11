#!/bin/bash

for prot in $@
do
	python3 main.py $prot
	echo "$prot done"
done

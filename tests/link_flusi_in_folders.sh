#!/bin/bash

# create a link to flusi in all folders
# useful if you intent to re-generate some tests

for dir in *
do
	if [ -d $dir ]; then
		echo "found directory" $dir
		cd $dir
		ln -s ../flusi
		cd ..
	fi
done

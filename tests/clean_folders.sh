#!/bin/bash

# clean the folders you find here (reduces the size of unit test tarball)

for dir in *
do
	if [ -d $dir ]; then
		echo "found directory" $dir
		cd $dir
		rm -f succ* *.rigidsolver decomp* divu.t mask_volume.t runtime_control.ini flusi *.h5 *.xmf return
		cd ..
	fi
done

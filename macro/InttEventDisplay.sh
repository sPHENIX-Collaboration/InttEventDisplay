#!/bin/bash
#echo "filepath " $1
#echo "ncluster " $2
#echo "savepicture " $3
SaveDirectory=/sphenix/u/mfujiwara/Workspace/tutorials/inttgitclone/work/pictures/
if [ -e "$1" ]; then
    root -l -n 'Loadfile.C("'$1'",'$2','$3','$4'"'$SaveDirectory'")'
else
    echo "$1 is not found"
fi

#!/bin/bash
#echo "filepath " $1
#echo "ncluster " $2
#echo "savepicture " $3
if [ -e "$1" ]; then
    root -l 'Loadfile.C("'$1'",'$2','$3')'
else
    echo "$1 is not found"
fi

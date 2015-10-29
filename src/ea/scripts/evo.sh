#!/bin/bash -f

#evoScript for running several independent runs in implicit
#encoding in python
#to launch from implicit_genome.py folder
usage="Usage:`basename $0` propertiesFile nbRuns"

if [ "$#" -lt 2 ] ; then
    echo $usage
    exit 1;
fi
mkdir logs
params="$1";

let "nbrun= $2";

out_global_dir=`echo $1 | awk 'BEGIN { FS= "[/.]"} ; { print $2}';`;
echo $out_global_dir

mkdir $out_global_dir;

i=1
    while [ $i -le ${nbrun} ]; do
	rm -rf logs/$out_global_dir
	mkdir logs/$out_global_dir
	python implicit_genome.py ${params}
	
	out_dir=`printf "out-%02d" $i`
	rm -rf $out_dir
        #gzip a faire plus tard dans le programme (mtn je peux pas -> experiences en route)
	#penser a supprimer l affichage texte du simulateur
	
	mv logs/$out_global_dir "$out_global_dir/$out_dir"
	
	i=$[i+1]

    done
#gzip -r "$out_global_dir/$out_dir"
mkdir logs

exit 0;


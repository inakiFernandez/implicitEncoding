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
#echo $out_global_dir
rm -rf logs/$out_global_dir

mkdir logs/$out_global_dir

out_parallel_file=`echo $out_global_dir.par`
touch $out_parallel_file
i=1
    while [ $i -le ${nbrun} ]; do
	out_dir=`printf "out-%02d" $i`
	rm -rf $out_dir
	mkdir logs/$out_global_dir/$out_dir
	#python implicit_genome.py ${params} $out_dir
	echo "python implicit_genome.py ${params} $out_dir" >> $out_parallel_file

        #gzip a faire plus tard dans le programme (mtn je peux pas -> experiences en route)
	#penser a supprimer l affichage texte du simulateur
	
	#mv logs/$out_global_dir/$out_dir "$out_global_dir/$out_dir"
	
	i=$[i+1]

    done
#gzip -r "$out_global_dir/$out_dir"

#mv $out_parallel_file $out_global_dir
parallel -j 4 < $out_parallel_file
mv logs/$out_global_dir $out_global_dir

mkdir logs




exit 0;


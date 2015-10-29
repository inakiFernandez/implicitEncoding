#!/bin/bash -f

#expScript for running statistical computations on the results of foraging RoboRobo MONEE 

usage="Usage:`basename $0` configFile nbRuns itemLogName"

if [ "$#" -lt 3 ] ; then
    echo $usage
    exit 1;
fi


params="-l $1";

logfile="$3";

let "nbrun= $2";

#i=1
#out_files="";

#while [ $i -le ${nbrun} ]; do

#	new_dir=`printf "out-%02d/" $i`
#	gunzip -r ${new_dir}; #or zcat later instead ?
#	out_files="${out_files}${new_dir}${logfile} ";
	
	#rm -rf $out_dir
	#mv logs $out_dir
	#mv $controller $out_dir
#	i=$[i+1]

#done
#echo ${out_files};


out_global_dir=`echo $1 | awk 'BEGIN { FS= "[/.]"} ; { print $2}';`;
i=1
new_dir=`printf "out-%02d/" $i`;
out_file="${out_global_dir}/${new_dir}${logfile}";


cat $out_file | awk '{ print $1}' >> global1.item.log;
i=2

while [ $i -le ${nbrun} ]; do

	new_dir=`printf "out-%02d/" $i`
	out_file="${out_global_dir}/${new_dir}${logfile}";
	echo $out_file;
	cat $out_file | awk '{ print $1}' >> tmp1.log;

	paste -d " " global1.item.log tmp1.log >> tmp.global1.log;	
	mv tmp.global1.log global1.item.log;
	rm tmp1.log;

	i=$[i+1]
done

mv global1.item.log $out_global_dir;

i=1
new_dir=`printf "out-%02d/" $i`;
out_file="${out_global_dir}/${new_dir}${logfile}";


cat $out_file | awk '{ print $2}' >> global2.item.log;
i=2

while [ $i -le ${nbrun} ]; do

	new_dir=`printf "out-%02d/" $i`
	out_file="${out_global_dir}/${new_dir}${logfile}";

	cat $out_file | awk '{ print $2}' >> tmp2.log;

	paste -d " " global2.item.log tmp2.log >> tmp.global2.log;	
	mv tmp.global2.log global2.item.log;
	rm tmp2.log;

	i=$[i+1]
done

mv global2.item.log $out_global_dir;

exit 0;


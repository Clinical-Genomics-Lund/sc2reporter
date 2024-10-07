#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "No run specified"
	exit
fi

run="${1%/}"
run="${run##*/}"

machine="${run:7:8}"

case $machine in
	NB501697)
		machinepath="/data/NextSeq1"
		;;
	NB501699)
		machinepath="/data/NextSeq2"
		;;
	A00681*)
		machinepath="/fs2/seqdata/NovaSeq"
		;;
	A01932*)
		machinepath="/fs2/seqdata/NovaSeq"
		;;
esac

echo $machinepath

export found=0
export not_found=0

while read i ; do
	if [ -f ${machinepath}/${run##*/}/Data/Intensities/BaseCalls/*${i}*R1*gz ] ; then
		let "found = found + 1"
		# echo $i exists
	else
		let "not_found = not_found + 1"
		echo $i has not been demultiplexed
	fi
done < <(grep sars ${machinepath}/${run}/SampleSheet.csv | cut -d',' -f1)

echo $found found
echo $not_found not demultiplexed


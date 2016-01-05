#!/bin/bash

set -e 


for thisDir in data/*;do
	echo $thisDir
	virusFile=$thisDir/*.fasta
	echo $virusFile
	for ii in $thisDir/*.fastq.gz;do
    echo $ii
    base=$(basename $ii)
    outFile=work/virusFind/${base%_trim.fastq.gz}
		if [ ! -f "${outFile}_match.fastq.gz" ];then
			date
      ~/installs/suffixc/suffixc $virusFile $ii -o $outFile -m 1 -l 12
      date
    else
			echo "Already processed. Skipping."
    fi
  done
done

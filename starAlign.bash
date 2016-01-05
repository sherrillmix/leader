#!/bin/bash 

#make errors end script
set -e 

for ii in data/*/*.fastq.gz;do
  echo $ii
  base=$(basename $ii)
  outFile=work/align/${base%_trim.fastq.gz}
  finalFile=$outFile.bam
  if [ ! -f "$finalFile" ];then
    date
    ./removeShort.py $ii >work/unzipped.fastq
    echo "Unzipped"
    date
    echo "Aligning"
    if [[ $ii =~ SIV ]];then
      echo "Aligning to monkey"
      index=~/installs/star/mac/
    else
      echo "Aligning to human"
      index=~/installs/star/index/
    fi
    ~/installs/star/bin/Linux_x86_64/STAR --genomeDir $index  --runThreadN 24 --readFilesIn work/unzipped.fastq --outFileNamePrefix $outFile.
    echo "Sorting"
    samtools view -bS "$outFile.Aligned.out.sam" > work/tmp.bam
    samtools sort -m 5000000000 work/tmp.bam $outFile #use ~5G of RAM
    samtools index $finalFile
    rm "$outFile.Aligned.out.sam"
    rm work/unzipped.fastq
    rm work/tmp.bam
    echo "Done"
  else
    echo "Already processed. Skipping."
  fi
done


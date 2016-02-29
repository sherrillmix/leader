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
    #from https://github.com/sherrillmix/dnapy/blob/master/dnapy/removeshort.py
    removeShort $ii -d 100000 -l 15>work/${base}_host_unzipped.fastq
    echo "Unzipped"
    date
    echo "Aligning"
    if [[ $ii =~ SIV ]];then
      echo "Aligning to monkey"
      index=index/mac/
    else
      echo "Aligning to human"
      index=index/index/
    fi
    ~/installs/star/bin/Linux_x86_64/STAR --genomeDir $index  --runThreadN 24 --readFilesIn work/${base}_host_unzipped.fastq --outFileNamePrefix $outFile.
    echo "Sorting"
    samtools view -bS "$outFile.Aligned.out.sam" > work/${base}_host_tmp.bam
    samtools sort -m 9000000000 work/${base}_host_tmp.bam $outFile 
    samtools index $finalFile
    rm "$outFile.Aligned.out.sam"
    rm work/${base}_host_unzipped.fastq
    rm work/${base}_host_tmp.bam
    echo "Done"
  else
    echo "Already processed. Skipping."
  fi
done


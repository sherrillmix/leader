#!/bin/bash 

#make errors end script
set -e 

for thisDir in data/*;do
	echo $thisDir
	virusFile=$thisDir/*.fasta
	echo $virusFile
	virusDir=work/$(basename $virusFile)
	virusDir=${virusDit%.fasta}
	echo $virusDir
	#if [ ! -e 
	#~/installs/star/bin/Linux_x86_64/STAR   --runMode genomeGenerate   --runThreadN 24   --genomeDir ./ --sjdbGTFfile knownGene.gtf  --genomeFastaFiles *.fa

	for ii in $thisDir/*.fastq.gz;do
		continue
		echo $ii
		base=$(basename $ii)
		outFile=work/align/${base%_trim.fastq.gz}_virus
		finalFile=$outFile.bam
		if [ ! -f "$finalFile" ];then
			date
			./removeShort.py $ii >work/unzipped.fastq
			echo "Unzipped"
			date
			echo "Aligning"
			~/installs/star/bin/Linux_x86_64/STAR --genomeDir ~/installs/star/index/  --runThreadN 24 --readFilesIn work/unzipped.fastq --outFileNamePrefix $outFile.
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
done


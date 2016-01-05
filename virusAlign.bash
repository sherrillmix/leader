#!/bin/bash 

#make errors end script
set -e 

for thisDir in data/*;do
	echo $thisDir
	virusFile=$thisDir/*.fasta
	echo $virusFile
	virusDir=work/virusGenomes/$(basename $virusFile)
	virusDir=${virusDir%.fasta}
	if [ ! -e $virusDir ];then
		echo "Working on $virusDir"
		mkdir -p $virusDir
		#map HIV89.6 junctions virus? try just mapping total first
		#~/installs/star/bin/Linux_x86_64/STAR   --runMode genomeGenerate   --runThreadN 24   --genomeDir ./ --sjdbGTFfile knownGene.gtf  --genomeFastaFiles *.fa --genomeSAindexNbases 4
		~/installs/star/bin/Linux_x86_64/STAR   --runMode genomeGenerate   --runThreadN 24   --genomeDir $virusDir --genomeFastaFiles $virusFile --genomeSAindexNbases 4 --sjdbFileChrStartEnd $thisDir/sjdbFile.tsv
	else
		echo "$virusDir genome already built"
	fi
	for ii in $thisDir/*.fastq.gz;do
		echo $ii
		base=$(basename $ii)
		outFile=work/virusAlign/${base%_trim.fastq.gz}_virus
		finalFile=$outFile.bam
		if [ ! -f "$finalFile" ];then
			date
			./removeShort.py $ii >work/unzipped_virus.fastq
			echo "Unzipped"
			date
			echo "Aligning"
			~/installs/star/bin/Linux_x86_64/STAR --genomeDir $virusDir  --runThreadN 20 --readFilesIn work/unzipped_virus.fastq --outFileNamePrefix $outFile.
			echo "Sorting"
			samtools view -bS "$outFile.Aligned.out.sam" > work/tmp_virus.bam
			samtools sort -m 5000000000 work/tmp_virus.bam $outFile #use ~5G of RAM
			samtools index $finalFile
			rm "$outFile.Aligned.out.sam"
			rm work/unzipped_virus.fastq
			rm work/tmp_virus.bam
			echo "Done"
		else
			echo "Already processed. Skipping."
		fi
	done
done


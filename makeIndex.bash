#!/bin/bash
set -e
cd index/index
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/ .
zcat ../knownGene.gtf.gz > knownGene.gtf  #downloaded from ucsc table browser
for ii in *.fa.gz;do
  zcat $ii > ${ii%.gz}
done
../bin/Linux_x86_64/STAR   --runMode genomeGenerate   --runThreadN 24   --genomeDir ./ --sjdbGTFfile knownGene.gtf  --genomeFastaFiles *.fa
rm *.fa
cd ../..

cd index/mac
rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/rheMac3/bigZips/rheMac3.fa.gz .
zcat ../rheMac3_refSeq.gtf.gz > refSeq.gtf #downloaded from ucsc table browser
zcat rheMac3.fa.gz > rheMac3.fa
../bin/Linux_x86_64/STAR   --runMode genomeGenerate   --runThreadN 12   --genomeDir ./ --sjdbGTFfile refSeq.gtf --limitGenomeGenerateRAM 91000000000  --genomeFastaFiles *.fa 
rm *.fa


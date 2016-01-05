# leader

## findSpliceSites.R 
Take annotated splice donors and acceptors, find on HIV/SIV genome, make sure cannonical and output all appropriate pairwise combinations for input into STAR aligner
## findVirusSeqs.bash 
Use [suffixc](https://github.com/sherrillmix/suffixc) to find reads starting/ending with SIV/HIV sequences (potential filter to speed up downstream alignment or look for chimeric RNA)
## makeIndex.bash 
Make macque/human indices for STAR aligner
## removeShort.py 
Remove short/empty reads (which trip up STAR aligner)
## starAlign.bash
Align to human/macaque genome using STAR aligner
## virusAlign.bash
Align to SIV/HIV genome using STAR aligner

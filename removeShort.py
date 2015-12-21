#!/usr/bin/env python
import sys
import getopt
import gzip
import atexit
import Bio.SeqIO.QualityIO


class Usage(Exception):
    def __init__(self,msg):
        self.msg=msg


def usage():
    return 'Usage: python removeShort.py [-v] -o outPrefix myFile.fastq\n'


def closeFiles(openFiles):
    for openFile in openFiles:
        openFile.close()

def writeRead(readFile,read):
    readFile.write("@%s\n%s\n+%s\n%s\n" %(read[0],read[1],read[0],read[2]))

def main(argv=None):
    if argv is None:
        argv=sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:],'hvd:o:l:',['help','output=','dots=','minLength='])
        except getopt.error, msg:
            raise Usage(msg)
    except Usage, err:
        sys.stderr.write(err.msg+"\n")
        sys.stderr.write(usage())
        return 2

    outFile='out_'
    verbose=True
    dots=1e6
    minLength=15
    for opt, arg in opts:
        if opt == "-v":
            verbose=True
        elif opt in ("-l","--minLength"):
            dots=int(arg)
        elif opt in ("-d","--dots"):
            dots=int(arg)
        elif opt in ("-h", "--help"):
            sys.stderr.write(usage())
            sys.exit()
        else:
            assert False, "unhandled option"

    if len(args)!=1:
        sys.stderr.write('Please provide a single fastq or fastq.gz file\n')
        return 4
    try:
        if args[0][-2:]=='gz' or args[0][-4]=='gzip':
            fastq=gzip.open(args[0], "r")
        else:
            fastq=open(args[0],"r")
    except IOError:
        sys.stderr.write("Problem opening file:"+args[0]+"\n")
        return 3
    openFiles=[fastq]
    atexit.register(closeFiles,openFiles)
    

    #assuming reads are sorted so that if present read2 would follow read1
    nGood=0
    nBad=0
    for currentRead in Bio.SeqIO.QualityIO.FastqGeneralIterator(fastq):
        if len(currentRead[1])>=minLength:
            nGood+=1
            writeRead(sys.stdout,currentRead)
        else:
            nBad+=1
        if verbose and (nGood+nBad) % dots == 0:
            sys.stderr.write('.')
    if verbose:
        sys.stderr.write("\nGood reads: "+str(nGood)+" Bad reads: "+str(nBad)+"\n")


if __name__ == '__main__':
    sys.exit(main())

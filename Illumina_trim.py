# this python script will read in forward reads, reverse reads and a tagList
# it splits file into interleaved fasta file sets with forward read followed by reverse read
# splits occur by tag in tagList.txt
# for paired reads to be interleaved into output, they must both match tag exactly
# If neither tags match, the reads (forward and reverse) are omitted completely
# It includes an option to keep reads that have mismatching barcodes ont eh forward and reverse read (saveMismatches)

import sys
import os
import argparse
import time

def reverseComplement(str):
    revcompDict={}
    revcompDict['A']='T'
    revcompDict['C']='G'
    revcompDict['G']='C'
    revcompDict['T']='A'
    revcompDict['N']='N'
    ret=""
    for character in list(str[::-1]):
        ret+=revcompDict[character]
    return(ret)

#arg defaults
minQC=15
phredOffset=33
logFileName="Paired_MiSeq_trim_QC.log"
minLength=50
fivePrimeBuffer=8


parser = argparse.ArgumentParser(description='Options for readSorterTrimmer.')
parser.add_argument('--forward','-f', dest='forward',help='REQUIRED. The name of the forward paired end read file',required=True)
parser.add_argument('--reverse','-r', dest='reverse',help='REQUIRED. The name of the reverse paired end read file',required=True)
parser.add_argument('--output','-o', dest='output',help='REQUIRED. The name of the output file prefix',required=True)
parser.add_argument('--log','-l', dest='logFileName', default=logFileName, help='outputLog {default:'+logFileName+'}')
parser.add_argument('--minLength','-m', dest='minLength', default=minLength, type=int,help='minimum length requirement for both reads {default:'+str(minLength)+'}')
parser.add_argument('--minQC','-q', dest='minQC', default=minQC, type=int,help='minimum QC allowed {default:'+str(minQC)+'}')
parser.add_argument('--phred', '-p',dest='phredOffset', default=phredOffset, type=int,help='Phred offset used in Illumina fastq files.{default:'+str(phredOffset)+'}')
parser.add_argument('--fivePrimeBuffer', '-b',dest='fivePrimeBuffer', default=phredOffset, type=int,help='Length of sequence at five prime end of read that is masked from QC trim.{default:'+str(fivePrimeBuffer)+'}')

args=parser.parse_args()
pairedEnd_file_name_forward=args.forward
pairedEnd_file_name_reverse=args.reverse
outputFastqName_forward=args.output+"_1.fastq"
outputFastqName_reverse=args.output+"_2.fastq"
minQC=args.minQC
minQLength=args.minLength
phredOffset=args.phredOffset
logFileName=args.logFileName
fivePrimeBuffer=args.fivePrimeBuffer

#open logfile for writing
logOutput=open(logFileName,"w")
fastqOutput_1=open(outputFastqName_forward,"w")
fastqOutput_2=open(outputFastqName_reverse,"w")

#open paired_end_files and barcode file simultaneously and read through them
readsProcessed=0
rawPairedReadCount=0
rawPairedReadBases=0
finalPairedReadCount=0
finalPairedReadBases=0
print("Processing reads...")
startTime=time.time()
with open(pairedEnd_file_name_forward) as forwardReadFile, open(pairedEnd_file_name_reverse) as reverseReadFile: 
    for taxa1,fasta1,taxa1_,qual1,taxa2,fasta2,taxa2_,qual2 in zip(forwardReadFile,forwardReadFile,forwardReadFile,forwardReadFile,reverseReadFile,reverseReadFile,reverseReadFile,reverseReadFile):
        taxa1=taxa1.rstrip()
        taxa2=taxa2.rstrip()
        fasta1=fasta1.rstrip()
        fasta2=fasta2.rstrip()
        qual1=qual1.rstrip()
        qual2=qual2.rstrip()
        readsProcessed+=1
        rawPairedReadCount+=1
        rawPairedReadBases+=len(fasta1)+len(fasta2)
#check QC and trim at first instance of a quality score lower than minQC
        if ord(qual1[fivePrimeBuffer]) < phredOffset+minQC:
            fasta1=""
        else:
            qual1=qual1[:fivePrimeBuffer]+"".join((c if ord(c) >= phredOffset+minQC else ' ' for c in qual1[fivePrimeBuffer:])).split()[0]
            fasta1=fasta1[0:len(qual1)]
        if ord(qual2[fivePrimeBuffer]) < phredOffset+minQC:
            fasta2=""
        else:
            qual2=qual2[:fivePrimeBuffer]+"".join((c if ord(c) >= phredOffset+minQC else ' ' for c in qual2[fivePrimeBuffer:])).split()[0]
            fasta2=fasta2[0:len(qual2)]
#write files
        if len(fasta1)>=minLength and len(fasta2)>=minLength:
            finalPairedReadCount+=1
            finalPairedReadBases+=len(fasta1)+len(fasta2)
            fastqOutput_1.write(taxa1+"\n")
            fastqOutput_1.write(fasta1+"\n")
            fastqOutput_1.write("+\n")
            fastqOutput_1.write(qual1+"\n")
            fastqOutput_2.write(taxa2+"\n")
            fastqOutput_2.write(fasta2+"\n")
            fastqOutput_2.write("+\n")
            fastqOutput_2.write(qual2+"\n")
print("completed processing reads.\n"+str(readsProcessed)+" reads processed in "+"{:1.1f}".format(time.time()-startTime)+" sec.")

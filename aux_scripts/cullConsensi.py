"""
Created on Tue Jan  7 10:45:32 2020

@author: balemj
"""

import os
import sys

from Bio import Seq, SeqIO, SeqRecord


def getSeq(inPath):
    sequence = list(SeqIO.parse(inPath, 'fasta'))[0]
    
    return(sequence)
#/def getSeq
    

def main(args):
    inPath = args[0]
    outPath = args[1]
    dirs = [fol for fol in os.listdir(inPath) if
            os.path.isdir(os.path.join(inPath, fol))]
    
    newNames = ["index_" + x + "_shiverBuild_consensus_MinCov_15_30.fasta" for x in dirs]
    
    seqsExist = [os.path.join(inPath, y, x) for x, y  in 
                  zip(newNames, dirs) if
                  os.path.exists(os.path.join(inPath, y, x))]
    
    seqs_unedit = [getSeq(x) for x in seqsExist ]
    
    seqs_edit = [SeqRecord.SeqRecord(
            Seq.Seq(str(x.seq).replace('?', '-')),
            id=x.id,
            description=''
            ) for x in seqs_unedit]
    seqs_final = [x for x in seqs_edit if len(x.seq) > 250]

    with open(outPath, 'w') as outHandle:
        SeqIO.write(seqs_final, outHandle, "fasta")

#/def main
        
        
if __name__ == "__main__":
    main(sys.argv[1:])
#/end if
    
#/end cullConsensi.py

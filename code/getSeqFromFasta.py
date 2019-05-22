#!/usr/bin/python
#encoding=utf-8

import argparse
from pyfasta import Fasta
import pandas as pd

def fuc(row):
    row = row['seq'][0]
    return row

def getSequence(genome):
    genome=Fasta(genome)
    RAD_seq = pd.read_csv('../data/input_data/peak.csv')
    result = map(lambda i:[genome.sequence({'chr':RAD_seq['chrom'][i],'start':RAD_seq['start'][i],'stop':RAD_seq['end'][i]})],range(len(RAD_seq)))
    RAD_seq['seq'] = result
    RAD_seq['seq'] = RAD_seq.apply(fuc,axis=1)
    RAD_seq.to_csv('../data/input_data/RAD_seq.csv',index=False)
    print 'getSequence is over,RAD_seq.csv is bulit!'

if __name__=='__main__':
    genome = '../model/hg19.fa'
    getSequence(genome)
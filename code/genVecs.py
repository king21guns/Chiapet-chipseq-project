#!/usr/bin/python
#coding=utf-8

from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
import logging
import pandas as pd
import numpy as np

def DNA2Sentence(Dna, K):
    sentence = ''
    length = len(Dna)
    for i in xrange(length-K+1):
        if 'N' in Dna[i:i+K]:break
        sentence += Dna[i:i+K] + ' '
    #delete extra space
    sentence = sentence[0:len(sentence) - 1]
    return sentence

def genCorpus(K):
    RAD_seq = pd.read_csv('../data/input_data/RAD_seq.csv')
    for index,row in RAD_seq.iterrows():
        Dna = row['seq'].upper()
        sentence = DNA2Sentence(Dna,K)
        with open('../model/corpus.txt','a') as fr:
            fr.write(sentence+'\n')
    print 'corpus.txt is bulit'

def getWord_model():
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
    sentence = LineSentence('../model/corpus.txt',max_sentence_length = 10000)
    print "Start Training Word2Vec model..."

    # Set values for various parameters
    num_features = 100    # Word vector dimensionality
    min_word_count = 1      # Minimum word count
    num_workers = 20         # Number of threads to run in parallel
    context = 20            # Context window size
    downsampling = 1e-3     # Downsample setting for frequent words 

    # Initialize and train the model
    print "Training Word2Vec model..."
    word_model = Word2Vec(sentence, workers=num_workers,\
                    size=num_features, min_count=min_word_count, \
                    window=context, sample=downsampling, seed=1,iter = 50)
    word_model.save('../model/DNAtoVec.model')
    print 'DNAtoVec.model is bulit,Training is over!'

if __name__ == '__main__':
    genCorpus(6)
    getWord_model()
#!/usr/bin/python
#coding=utf-8

import pandas as pd
import numpy as np
import warnings
import sys
import os

def gethead(head1,head2):
    head = []
    head.append('label')
    for name in head1:
        name = 'C1_'+name
        head.append(name)

    for name in head2:
        if name != 'chrom':
            name = 'C2_' + name
        head.append(name)
    return head

def getRAD1(RAD,index):
    temp = RAD.iloc[index,0:]
    nonpredictors = ['chrom','used']
    predictors_df = temp.drop(nonpredictors,axis = 1)
    head = list(predictors_df.columns.values)

    return np.asarray(predictors_df),head

def getRAD2(RAD,index):
    temp = RAD.iloc[index,0:]
    nonpredictors = ['used']
    predictors_df = temp.drop(nonpredictors,axis = 1)
    head = list(predictors_df.columns.values)
    return np.asarray(predictors_df),head

def gen_neg():
    warnings.filterwarnings("ignore")
    RAD = pd.read_csv('../data/input_data/RAD_after_pos.csv')
    Positive_data = pd.read_csv('../data/input_data/pos_sample.csv')

    Positive_data = Positive_data[Positive_data['label'] == 1]
    Positive_data.index = xrange(len(Positive_data))

    distance = Positive_data["C2_start"] - Positive_data["C1_start"]

    distance = np.sort(distance,axis = None)

    #covering 95% percent of positive samples. Dump the 'extreme' situation
    a = int(len(Positive_data)*0.025)
    #a = 0
    min_range = distance[a]
    max_range = distance[len(distance)-a - 1]
    
    print 'min_range:%d' %(min_range)
    print 'max_range:%d' %(max_range)
    
    total = len(RAD)
    RAD.index = xrange(total)
    
    count = 0
    
    indexlist1=[]
    indexlist2=[]
    label = []

    for i in xrange(total):
        if (i % 100 == 0) and (i != 0):
                print "Processing: %d of %d\r" %(i,total),
                sys.stdout.flush()

        chrom = RAD['chrom'][i]
        start = RAD['start'][i]

        potential_pair = RAD[(RAD['chrom'] == chrom) & (RAD['start'] <= start + max_range) & (RAD['start'] >= start + min_range)] ###根据loop长短取负例

        if len(potential_pair) != 0:
            for j in xrange(len(potential_pair)):
                indexlist1.append(i)
                indexlist2.append(potential_pair.index[j])
                if (RAD['used'][i] == 1) & (RAD['used'][potential_pair.index[j]] == 1):
                    label.append(0)
                elif (RAD['used'][i] != 1) & (RAD['used'][potential_pair.index[j]] != 1):
                    label.append(3)
                elif RAD['used'][i] == 1:
                    label.append(4)
                else:
                    label.append(2)

                count +=1
    print
    print "Negative Data Source: " + str(count)

    indexlist1 = np.asarray(indexlist1)
    indexlist2 = np.asarray(indexlist2)
    label = np.asarray(label)
    label = label.reshape((len(label),1))
    a,head1 =  getRAD1(RAD,indexlist1)
    b,head2 = getRAD2(RAD,indexlist2)
    
    arrays = np.concatenate((label,a,b,np.abs(indexlist2-indexlist1).reshape((len(label),1))),axis = 1)
    header = gethead(head1,head2)
    header.append('num_between')
    table= pd.DataFrame(arrays,columns=header)
    table.to_csv('../data/input_data/neg_sample.csv', index = False)

    indexlist1 = indexlist1.reshape((len(indexlist1),1))
    indexlist2 = indexlist2.reshape((len(indexlist2),1))
    indexs = np.concatenate((label,indexlist1,indexlist2),axis = 1)
    indexs = pd.DataFrame(indexs,columns=['label','index1','index2'])
    indexs.to_csv('../data/input_data/neg_index.csv', index = False)
    
if __name__ == '__main__':
    gen_neg()
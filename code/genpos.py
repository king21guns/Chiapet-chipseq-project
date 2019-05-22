#!/usr/bin/python
#coding=utf-8

import numpy as np
import pandas as pd
import sys
import os
import warnings

chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

threshold = 0
#####################################################################################################################
def RAD21_active():
    loop = pd.read_csv('../data/input_data/loop.csv')

    chrom = np.concatenate((loop['c1_chr'],loop['c2_chr']),axis = 0)
    start = np.concatenate((loop['c1_start'],loop['c2_start']),axis = 0)
    end = np.concatenate((loop['c1_end'],loop['c2_end']),axis = 0)
    start -= threshold
    end += threshold
    peak = {'chrom':chrom,'start':start,'end':end}
    peak = pd.DataFrame(peak)
    peak['covered'] = 0

    RAD = pd.read_csv('../data/input_data/peak.csv')
    RAD['active'] = 0
    #RAD = pd.read_table("../Temp/%s/RAD.csv" %(cell),sep = ",")
    total = len(RAD)
    counter = 0

    for chr in chrom_list:
        active = []
        peak_value = []
        temp_RAD = RAD[RAD['chrom'] == chr]
        temp_peak = peak[peak['chrom'] == chr]

        for i in xrange(len(temp_RAD)):
            if (counter % 100 == 0) and (i != 0):
                print 'Generating Active segments: %d of %d\r' %(counter,total),
                sys.stdout.flush()
            index = temp_RAD.index[i]
            start = RAD['start'][index]
            end = RAD['end'][index]
            summit = (start + end)/2
            a = temp_peak[(temp_peak['start'] <= summit) & (temp_peak['end'] >= summit)]

            if(len(a) != 0):
                active.append(1)
                peak['covered'][a.index] = 1
            else:
                active.append(0)
            counter += 1
        RAD['active'][temp_RAD.index] = active
    print 'active chip-seq segments: %d of %d\r' %(len(RAD[RAD['active']==1]),len(RAD))
    print 'covered chiA-pet: %d of %d\r' %(len(peak[peak['covered']==1]),len(peak))
    peak.to_csv('../data/input_data/RAD21_chiA-pet_cover.csv',index=False)
    RAD.to_csv('../data/input_data/RAD21_chip-seq_active.csv',index=False)
#############################################################################################################################

def getcases(loop,i,RAD,origin_RAD):
    threshold = 0
    cases1 = "" #####why
    i = loop.index[i]
    while (len(cases1) == 0) and (threshold <= 10):
        chrom = loop['c1_chr'][i]
        start = loop['c1_start'][i]- threshold
        end = loop['c1_end'][i]+ threshold

        cases1 = RAD[(RAD['chrom'] == chrom) & (RAD['start'] >= start - 18) & (RAD['end'] <= end + 18)]
    
        if len(cases1) != 0:
            origin_RAD['used'][cases1.index] = 1
        threshold +=50
    
    threshold = 0
    cases2 = ""
    while (len(cases2) == 0) and (threshold <= 10):

        chrom = loop['c2_chr'][i]
        start = loop['c2_start'][i]- threshold
        end = loop['c2_end'][i]+ threshold

        cases2 = RAD[(RAD['chrom'] == chrom) & (RAD['start'] >= start - 18) & (RAD['end'] <= end + 18)]
        
        if len(cases2) != 0:
            origin_RAD['used'][cases2.index] = 1
        threshold +=50
        
    return cases1,cases2

def gethead(head1,head2):
    head = []
    head.append('label')
    for name in head1:
        name = 'C1_'+name
        head.append(name)

    for name in head2:
        if name != "chrom":
            name = 'C2_' + name
        head.append(name)

    print head
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

def concat_seq_histone():
    RAD_histone = pd.read_csv('../data/input_data/RAD_epigenetic.csv')
    RAD_seq = pd.read_csv('../data/input_data/RAD_seq.csv')
    RAD_seq = RAD_seq.drop(RAD_seq.columns[0:3],axis=1)
    RAD = pd.concat([RAD_histone,RAD_seq],axis=1)
    RAD.to_csv('../data/input_data/RAD.csv',index=False)

def gen_pos():
    warnings.filterwarnings("ignore")
    loop = pd.read_csv('../data/input_data/loop.csv')
    RAD = pd.read_csv('../data/input_data/RAD.csv')

    pos,unlab,upside,downside,unmap= 0,0,0,0,0
    indexlist1=[]
    indexlist2=[]
    label = []

    RAD['used'] = 0

    count = 0
    for chr in chrom_list:
        temp_RAD = RAD[RAD['chrom'] == chr]
        temp_loop = loop[loop['c1_chr'] == chr]

        for i in xrange(len(temp_loop)):
            if (count % 100 == 0) and (i != 0):
                print "Mapping: %d of %d\r" %(count,len(loop)),
                sys.stdout.flush()
            count += 1

            cases1,cases2 = getcases(temp_loop,i,temp_RAD,RAD)
            length1,length2 = len(cases1),len(cases2)

            if length1 == 0 and length2 == 0:
                unmap += 1
                continue
            elif length1 == 0:
                upside +=1
                continue
            elif length2 == 0:
                downside +=1
                continue
            elif length1 == 1 and length2 == 1:
                label.append(1)
                indexlist1.append(cases1.index[0])
                indexlist2.append(cases2.index[0])
                pos += 1
            else:
                for j in xrange(length1):
                    for k in xrange(length2):
                        label.append(-1)
                        indexlist1.append(cases1.index[j])
                        indexlist2.append(cases2.index[k])
                        unlab += 1

    print
    print "unmappable: %d  no upside: %d no downside : %d positive: %d potential: %d\n" %(unmap,upside,downside,pos,unlab)
    
    indexlist1 = np.asarray(indexlist1)
    indexlist2 = np.asarray(indexlist2)
    
    label = np.asarray(label)
    label = label.reshape((len(label),1))
    
    a,head1 =  getRAD1(RAD,indexlist1)
    b,head2 = getRAD2(RAD,indexlist2)

    indexlist1 = indexlist1.reshape((len(indexlist1),1))
    indexlist2 = indexlist2.reshape((len(indexlist2),1))
    indexs = np.concatenate((label,indexlist1,indexlist2),axis = 1)
    indexs = pd.DataFrame(indexs,columns=['label','index1','index2'])
    indexs.to_csv('../data/input_data/pos_index.csv',index = False)
    arrays = np.concatenate((label,a,b,np.abs(indexlist2-indexlist1)),axis = 1)
    header = gethead(head1,head2)
    header.append('num_between')
    table= pd.DataFrame(arrays,columns=header)
    table.to_csv('../data/input_data/pos_sample.csv',index = False) 

    RAD.to_csv('../data/input_data/RAD_after_pos.csv',index = False)

if __name__ == '__main__':
    concat_seq_histone()
    gen_pos()
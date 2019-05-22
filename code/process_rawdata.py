#!/usr/bin/python
#coding=utf-8

import pandas as pd
import sys

def overlap_RAD21(file_chiApet,file_chipsep):
    input_all = pd.read_csv(file_chipsep,sep='\t',header=None)
    input_all = input_all.drop(input_all.columns[3:],axis=1)
    col_name = ['chrom','start','end']
    input_all.columns = col_name
    input_all.sort_values(by=['chrom','start'],inplace=True)
    input_all.index = range(len(input_all))
    input_all.to_csv('../data/input_data/peak.csv',index=False)
    print 'peak.csv is built!'
    input_all['mid'] = (input_all['start']+input_all['end'])/2

    pos = pd.read_excel(file_chiApet)
    pos = pos.drop(pos.columns[6:],axis=1)
    col_name = ['c1_chr','c1_start','c1_end','c2_chr','c2_start','c2_end']
    pos.columns = col_name
    pos.sort_values(by=['c1_chr','c1_start'],inplace=True)
    pos.index = range(len(pos))
    pos.to_csv('../data/input_data/loop.csv',index=False)
    print 'loop.csv is built!'
    pos['c1_mid'] = (pos['c1_start']+pos['c1_end'])/2
    pos['c2_mid'] = (pos['c2_start']+pos['c2_end'])/2

    dis_c1 = []
    dis_c2 = []
    total = len(pos)
    for index,row in pos.iterrows():
        if index%1000==0:
            print 'dealing with: %d of %d\r' %(index,total),
            sys.stdout.flush()
        temp = abs(input_all[input_all['chrom']==row['c1_chr']]['mid']-row['c1_mid'])
        dis_c1.append(min(temp))
        temp = abs(input_all[input_all['chrom']==row['c2_chr']]['mid']-row['c2_mid'])
        dis_c2.append(min(temp))        

    pos['c1_dis'] = dis_c1
    pos['c2_dis'] = dis_c2
    pos.to_csv('../data/input_data/RAD21_overlap_pair.csv',index=False)
    counter = len(pos[(pos['c1_dis']<250)|(pos['c2_dis']<250)])
    #total = len(pos)
    print 'overlap rows(either dis < 1000): %d of %d\r' %(counter,total)
    #sys.stdout.flush()

if __name__ == '__main__':
    overlap_RAD21('../data/raw_data/RAD21_loop_chiApet.xlsx','../data/raw_data/RAD21_peak_chipsep.bed')


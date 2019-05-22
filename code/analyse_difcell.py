#!/usr/bin/python
#coding=utf-8

import numpy as np
import pandas as pd
import sys
from matplotlib_venn import venn2, venn2_circles,venn3, venn3_circles
from matplotlib import pyplot as plt
import re
import os
import warnings

chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

threshold = 0

def pos_xlsx(file_name):
    pos = pd.read_excel('../other_cell/raw_data/%s/chiA-pet.xlsx' %(file_name))
    pos = pos.drop(pos.columns[6:],axis=1)
    col_name = ['c1_chr','c1_start','c1_end','c2_chr','c2_start','c2_end']
    pos.columns = col_name
    pos.sort_values(by=['c1_chr','c1_start'],inplace=True)
    pos.index = range(len(pos))
    pos.to_csv('../other_cell/processed/%s_chiA-pet.csv' %(file_name),index=False)
    print 'chiA-pet.csv is built! num:%d' %(len(pos))
    return pos

def pos_bed(file_name):
    pos = pd.read_csv('../other_cell/raw_data/bed_bed/%s/chiA-pet.bed' %(file_name),sep='\t',header=None)
    temp = pos[3].values
    pos = []
    for i in temp:
        pos.append(re.split(r'[\:\.\-\,]+',i))
    pos = pd.DataFrame(pos)
    pos = pos.drop(pos.columns[6:],axis=1) 
    col_name = ['c1_chr','c1_start','c1_end','c2_chr','c2_start','c2_end']
    pos.columns = col_name
    pos = pos[pos['c1_chr']==pos['c2_chr']]
    pos.index = range(len(pos))
    #pos.apply(pd.to_numeric, errors='ignore')
    pos[['c1_start','c1_end','c2_start','c2_end']] = pos[['c1_start','c1_end','c2_start','c2_end']].astype(int)
    pos.sort_values(by=['c1_chr','c1_start'],inplace=True)
    pos.to_csv('../other_cell/processed/%s_chiA-pet.csv' %(file_name),index=False)
    print 'chiA-pet.csv is built! num:%d' %(len(pos))
    return pos

def pos_fastq(file_name):
    file_list = os.listdir('../other_cell/raw_data/%s/chiA-pet-fastq/' %(file_name))
    chiApet = []
    for chiApet_name in file_list:
        pos = pd.read_csv('../other_cell/raw_data/%s/chiA-pet-fastq/%s' %(file_name,chiApet_name),sep='\t',header=None) 
        pos = pos.drop(pos.columns[6:],axis=1)
        col_name = ['c1_chr','c1_start','c1_end','c2_chr','c2_start','c2_end']
        pos.columns = col_name
        pos = pos[pos['c1_chr']==pos['c2_chr']]
        pos.index = range(len(pos))
        pos[['c1_start','c1_end','c2_start','c2_end']] = pos[['c1_start','c1_end','c2_start','c2_end']].astype(int)
        chiApet.append(pos)

    pos = chiApet[0]
    for i in range(len(chiApet)-1):
        pos = pd.concat([pos,chiApet[i+1]],axis=0)
    pos.sort_values(by=['c1_chr','c1_start'],inplace=True)
    pos.to_csv('../other_cell/processed/%s_chiA-pet.csv' %(file_name),index=False)
    print 'chiA-pet.csv is built! num:%d' %(len(pos))
    return pos

def cal_dis(file_name):
    input_all = pd.read_csv('../other_cell/raw_data/%s/chip-seq.bed' %(file_name),sep='\t',header=None)
    input_all = input_all.drop(input_all.columns[3:],axis=1)
    col_name = ['chrom','start','end']
    input_all.columns = col_name
    input_all.sort_values(by=['chrom','start'],inplace=True)
    input_all.index = range(len(input_all))
    input_all.to_csv('../other_cell/processed/%s_chip-seq.csv' %(file_name),index=False)
    print 'chip-seq.csv is built!'
    input_all['mid'] = (input_all['start']+input_all['end'])/2

    #chiA-pet is .xlsx
    #pos = pos_xlsx(file_name)

    #chiA-pet is .bed
    #pos = pos_bed(file_name)

    #chiA-pet is some .fastq
    pos = pos_fastq(file_name)

    pos['c1_mid'] = (pos['c1_start']+pos['c1_end'])/2
    pos['c2_mid'] = (pos['c2_start']+pos['c2_end'])/2

    pos['c1_length'] = pos['c1_end']-pos['c1_start']
    pos['c2_length'] = pos['c2_end']-pos['c2_start']

    dis_c1 = []
    dis_c2 = []
    total = len(pos)
    for index,row in pos.iterrows():
        if index%1000==0:
            print 'dealing with %s: %d of %d\r' %(file_name,index,total),
            sys.stdout.flush()
        temp = abs(input_all[input_all['chrom']==row['c1_chr']]['mid']-row['c1_mid'])
        if len(temp)>1:
            dis_c1.append(min(temp))
        else:
            dis_c1.append(1e+5)
        temp = abs(input_all[input_all['chrom']==row['c2_chr']]['mid']-row['c2_mid'])
        if len(temp)>1:
            dis_c2.append(min(temp))
        else:
            dis_c2.append(1e+5)       

    pos['c1_dis'] = dis_c1
    pos['c2_dis'] = dis_c2

    pos.to_csv('../other_cell/processed/%s_overlap_pair.csv' %(file_name),index=False)
    print '%s is done' %(file_name)
    #total = len(pos)
    #print 'overlap rows(either dis < 1000): %d of %d\r' %(counter,total)
    #sys.stdout.flush()
    return pos,input_all,file_name

def calculate(pos,input_all,file_name):
    c1_mean_len = int(pos['c1_length'].mean())
    c2_mean_len = int(pos['c2_length'].mean())
    dis_standard = (c1_mean_len+c2_mean_len)/2
    print 'dis_standard:%d bp' %(dis_standard)

    counter_single = len(pos[(pos['c1_dis']<dis_standard)|(pos['c2_dis']<dis_standard)])
    single = [counter_single,len(pos)-counter_single]
    single_ratio = float(len(pos)-counter_single)/counter_single

    counter_double = len(pos[(pos['c1_dis']<dis_standard)&(pos['c2_dis']<dis_standard)])
    double = [counter_double,len(pos)-counter_double]
    
    pos = pos.drop(['c1_length','c2_length','c1_mid','c2_mid'],axis=1)
    pos = pos[['c1_chr','c1_start','c1_end','c1_dis','c2_chr','c2_start','c2_end','c2_dis']]
    a = pos.drop(pos.columns[4:],axis=1)
    b = pos.drop(pos.columns[0:4],axis=1)
    b.columns=['c1_chr','c1_start','c1_end','c1_dis']
    c = pd.concat([a,b],axis=0)
    d = c.drop_duplicates(['c1_chr','c1_start'],keep='first')

    counter_no = len(d[d['c1_dis']<dis_standard])
    no_repeat_single = [counter_no,len(d)-counter_no]
    print 'calculate:%s single:%d of %d,double:%d of %d,no_repeat_single:%d of %d' %(file_name,counter_single,len(pos),counter_double,len(pos),counter_no,len(d))
    return single,double,no_repeat_single,single_ratio

def active(loop,input_all,file_name,single_ratio):
    chrom = np.concatenate((loop['c1_chr'],loop['c2_chr']),axis = 0)
    start = np.concatenate((loop['c1_start'],loop['c2_start']),axis = 0)
    end = np.concatenate((loop['c1_end'],loop['c2_end']),axis = 0)
    start -= threshold
    end += threshold
    peak = {'chrom':chrom,'start':start,'end':end}
    peak = pd.DataFrame(peak)
    peak['covered'] = 0

    RAD = input_all
    RAD['active'] = 0
    #CTCF = pd.read_table("../Temp/%s/CTCF.csv" %(cell),sep = ",")
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
    print 'active chip-seq segments of %s: %d of %d\r' %(file_name,len(RAD[RAD['active']==1]),len(RAD))
    print 'covered chiA-pet of %s: %d of %d\r' %(file_name,len(peak[peak['covered']==1]),len(peak))

    return [int(len(RAD[RAD['active']==1])*single_ratio),len(RAD)-len(RAD[RAD['active']==1]),len(RAD[RAD['active']==1])]

def venn_pie_plot(single,double,no_repeat_single,overlap,file_name):
    title_name = file_name + ' venn&pie diagram'
    labels_1 = ['mapped(%d)'%(single[0]),'unmapped(%d)'%(single[1])]
    labels_2 = ['mapped(%d)'%(double[0]),'unmapped(%d)'%(double[1])]
    labels_3 = ['mapped(%d)'%(no_repeat_single[0]),'unmapped(%d)'%(no_repeat_single[1])]
    fig = plt.figure(figsize=(12,12))
    plt.subplot(2,2,1) 
    plt.pie(single,labels=labels_1,autopct='%1.2f%%')
    plt.title('single mapped pie chart(chiA-pet)',bbox={'facecolor':'0.8', 'pad':5})
    plt.subplot(2,2,2) #两行两列,这是第二个图
    plt.pie(double,labels=labels_2,autopct='%1.2f%%')
    plt.title('double mapped pie chart(chiA-pet)',bbox={'facecolor':'0.8', 'pad':5})
    plt.subplot(2,2,3) 
    plt.pie(no_repeat_single,labels=labels_3,autopct='%1.2f%%')
    plt.title('no_repeat_single mapped pie chart(chiA-pet)',bbox={'facecolor':'0.8', 'pad':5})
    plt.subplot(2,2,4)
    venn2(subsets=overlap,set_labels=('chiA-pet','chip-seq'),set_colors=('r','g'))
    plt.title("chiA-pet&chip-seq venn diagram",bbox={'facecolor':'0.8', 'pad':5})
    fig.suptitle(title_name, fontsize=16)
    #plt.show()
    plt.savefig('../other_cell/pic/%s_venn&pie.jpg'%(title_name))
    plt.close()

if __name__ == '__main__':
    warnings.filterwarnings('ignore')
    file_name = 'HepG2-RAD21'
    pos,input_all,file_name = cal_dis(file_name)
    single,double,no_repeat_single,single_ratio = calculate(pos,input_all,file_name)
    overlap = active(pos,input_all,file_name,single_ratio)
    venn_pie_plot(single,double,no_repeat_single,overlap,file_name)
    '''
    file_list = '../other_cell/raw_data/bed_bed'
    file_list = os.listdir(file_list)
    for file_name in file_list:
        print '%s is dealing with' %(file_name)
        pos,input_all,file_name = cal_dis(file_name)
        single,double,no_repeat_single,single_ratio = calculate(pos,input_all,file_name)
        overlap = active(pos,input_all,file_name,single_ratio)
        venn_pie_plot(single,double,no_repeat_single,overlap,file_name)
    '''
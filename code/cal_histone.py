#!/usr/bin/python
#coding=utf-8

import os,sys
import pandas as pd
import numpy as np
import multiprocessing
from functools import partial
import time

chrom_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", \
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", \
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]

def gen_histone_matrix(RAD_fileList,histone_fileList):
    RAD = pd.read_csv(RAD_fileList)

    histone = os.listdir(histone_fileList)
    part = ['ups_','ins_','downs_']
    histone_columns = []
    for i in part:
        for j in histone:
            histone_columns.append(i+j)

    RAD = pd.concat([RAD,pd.DataFrame(np.zeros([RAD.shape[0], len(histone_columns)]), columns=histone_columns)],axis=1)
    print 'matrix is finished'
    return RAD,histone,histone_fileList

def count_histone_parallel(RAD,stream_num,histone_fileList,histone):
    #print 'couting is begining'
    #os.mkdir('temp')
    start_time = time.time()
    fileList = os.listdir(histone_fileList+'/'+histone)
    #counter = 0
    for j in fileList:
        file = pd.read_csv(histone_fileList+'/'+histone+'/'+j,sep='\t',header = None)
        header = ['chrom','start','end','name','no','+-']
        file.columns = header
        file['real_start'] = range(len(file))
        file.loc[file['+-']=='+','real_start']=file[file['+-']=='+']['start']
        file.loc[file['+-']=='-','real_start']=file[file['+-']=='-']['end']
        file = file.drop(['start','end','name','no','+-'],axis=1)
        now_time = time.time()
        print(now_time-start_time)

        for index, row in RAD.iterrows():
            '''
            line = line.strip().split('\t')
            if len(line)<3:break
            '''
            chrom = row['chrom']
            start = row['start']
            end = row['end']
            #if chrom == 'chrY':break
            #counter+=1
            
            if index%20000==0:
                RAD.to_csv('../result/cal_histone/temp/%s_%s_%d_RAD.csv' %(histone,j,index),index=False)
            
            #start,end = int(row['start']),int(row['end'])
            if index % 100 == 0:
                print 'DNA %d of %s in %s\r' %(index, j, histone)
                #print 'DNA %d of %s in %s\r' %(index, j, histone),
                now_time = time.time()
                print(now_time-start_time)
                #sys.stdout.flush()

            counter = file[(file['chrom']==chrom)&(file['real_start']>=start-stream_num)&(file['real_start']<=start)].index
            RAD.loc[index,'ups_'+histone]+=len(counter)

            counter = file[(file['chrom']==chrom)&(file['real_start']>=start)&(file['real_start']<=end)].index
            RAD.loc[index,'ins_'+histone]+=len(counter)

            counter = file[(file['chrom']==chrom)&(file['real_start']>=end)&(file['real_start']<=end+stream_num)].index
            RAD.loc[index,'downs_'+histone]+=len(counter)

    proc = os.getpid()
    print('{0} by process id: {1}'.format(histone, proc))
    if len(fileList)>0:
        RAD.to_csv('../result/cal_histone/final/%s_RAD.csv'%(histone),index=False)
    print '%s is finished' %(histone)
    return RAD

def count_histone_parallel_2(RAD,stream_num,histone_fileList,histone):
    #print 'couting is begining'
    #os.mkdir('temp')
    start_time = time.time()
    fileList = os.listdir(histone_fileList+'/'+histone)
    #counter = 0
    for j in fileList:
        file = pd.read_csv(histone_fileList+'/'+histone+'/'+j,sep='\t',header = None)
        header = ['chrom','start','end','name','no','+-']
        file.columns = header
        file['real_start'] = range(len(file))
        file.loc[file['+-']=='+','real_start']=file[file['+-']=='+']['start']
        file.loc[file['+-']=='-','real_start']=file[file['+-']=='-']['end']
        file = file.drop(['start','end','name','no','+-'],axis=1)
        now_time = time.time()
        print(now_time-start_time)

        total = len(RAD)
        count = 0
        for chr in chrom_list:
            ups_value = []
            ins_value = []
            downs_value = []
            temp_RAD = RAD[RAD["chrom"] == chr]
            temp_file = file[file["chrom"] == chr]

            for i in xrange(len(temp_RAD)):
                if (count % 100 == 0) and (i != 0):
                    print 'Dealing with %s: %d of %d\r' %(j,count,total),
                    sys.stdout.flush()
                index = temp_RAD.index[i]
                start = RAD['start'][index]
                end = RAD['end'][index]
                a= temp_file[(temp_file['real_start']>=start-stream_num)&(temp_file['real_start']<=start)]
                if(len(a) != 0):
                    ups_value.append(len(a))
                else:
                    ups_value.append(0)
                b= temp_file[(temp_file['real_start']>=start)&(temp_file['real_start']<=end)]
                if(len(b) != 0):
                    ins_value.append(len(b))
                else:
                    ins_value.append(0)
                c= temp_file[(temp_file['real_start']>=end)&(temp_file['real_start']<=end+stream_num)]
                if(len(c) != 0):
                    downs_value.append(len(c))
                else:
                    downs_value.append(0)
                count += 1
            RAD['ups_'+histone][temp_RAD.index] += ups_value
            RAD['ins_'+histone][temp_RAD.index] += ins_value
            RAD['downs_'+histone][temp_RAD.index] += downs_value

    proc = os.getpid()
    print('{0} by process id: {1}'.format(histone, proc))
    if len(fileList)>0:
        RAD.to_csv('../result/cal_histone/final/%s_RAD.csv'%(histone),index=False)
    print '%s is finished' %(histone)
    return RAD

def addToALL(result):
    RAD_notadd = result[0].drop(result[0].columns[3:],axis=1)
    RAD_data = result[0].drop(result[0].columns[0:3],axis=1)
    for i in range(len(result)-1):
        temp = result[i+1].drop(result[i+1].columns[0:3],axis=1)
        RAD_data += temp
    RAD = pd.concat([RAD_notadd,RAD_data],axis=1)
    #RAD.to_csv('../data/input_data/RAD_histone.csv',index=False)
    print 'all is finished'
    return RAD

def caL_DNase(RAD,stream_num):
    RAD['ups_DNase'] = 0
    RAD['DNase'] = 0
    RAD['downs_DNase'] = 0
    fileList = os.listdir('../data/raw_data/DNase')
    for j in fileList:
        file = pd.read_csv('../data/raw_data/DNase/'+j,sep='\t',header = None)
        file = file.drop(file.columns[3:],axis=1)
        header = ['chrom','start','end']
        file.columns = header


        total = len(RAD)
        count = 0
        for chr in chrom_list:
            DNase_ups_value = []
            DNase_value = []
            DNase_downs_value = []

            temp_RAD = RAD[RAD['chrom'] == chr]
            temp_file = file[file['chrom'] == chr]

            for i in xrange(len(temp_RAD)):
                if (count % 100 == 0) and (i != 0):
                    print 'Dealing with %s: %d of %d\r' %(j,count,total),
                    sys.stdout.flush()
                index = temp_RAD.index[i]
                start = RAD['start'][index]
                end = RAD['end'][index]
                a= temp_file[(temp_file['start']>=start-stream_num)&(temp_file['start']<=start)]
                if(len(a) != 0):
                    DNase_ups_value.append(len(a))
                else:
                    DNase_ups_value.append(0)
                b= temp_file[(temp_file['start']>=start)&(temp_file['start']<=end)]
                if(len(b) != 0):
                    DNase_value.append(len(b))
                else:
                    DNase_value.append(0)
                c= temp_file[(temp_file['start']>=end)&(temp_file['start']<=end+stream_num)]
                if(len(c) != 0):
                    DNase_downs_value.append(len(c))
                else:
                    DNase_downs_value.append(0)
                count += 1
            RAD['ups_DNase'][temp_RAD.index] += DNase_ups_value
            RAD['DNase'][temp_RAD.index] += DNase_value
            RAD['downs_DNase'][temp_RAD.index] += DNase_downs_value
    RAD.to_csv('../data/input_data/RAD_epigenetic.csv',index=False)
'''
def SaveTocsv(raw_pos_data,raw_neg_data,TotalMappedReads,output_path):
    raw_pos_data.to_csv(output_path+'raw_pos_data'+'.csv',sep='\t')
    raw_neg_data.to_csv(output_path+'raw_neg_data'+'.csv',sep='\t')
    TotalMappedReads.to_csv(output_path+'TotalMappedReads'+'.csv',sep='\t')
'''
################################################################################################################################
def count_histone_serial(RAD,stream_num,histone_fileList,histone_list):
    #print 'couting is begining'
    #os.mkdir('temp')
    start_time = time.time()
    for histone in histone_list:
        fileList = os.listdir(histone_fileList+'/'+histone)
        #counter = 0
        for j in fileList:
            file = pd.read_csv(histone_fileList+'/'+histone+'/'+j,sep='\t',header = None)
            header = ['chrom','start','end','name','no','+-']
            file.columns = header
            file['real_start'] = range(len(file))
            file.loc[file['+-']=='+','real_start']=file[file['+-']=='+']['start']
            file.loc[file['+-']=='-','real_start']=file[file['+-']=='-']['end']
            file = file.drop(['start','end','name','no','+-'],axis=1)
            now_time = time.time()
            print(now_time-start_time)

            for index, row in RAD.iterrows():
                '''
                line = line.strip().split('\t')
                if len(line)<3:break
                '''
                chrom = row['chrom']
                start = row['start']
                end = row['end']
                #if chrom == 'chrY':break
                #counter+=1
                
                if index%20000==0:
                    RAD.to_csv('../result/cal_histone/temp/%s_%s_%d_RAD.csv' %(histone,j,index),index=False)
                
                #start,end = int(row['start']),int(row['end'])
                if index % 100 == 0:
                    print 'DNA %d of %s in %s\r' %(index, j, histone)
                    #print 'DNA %d of %s in %s\r' %(index, j, histone),
                    now_time = time.time()
                    print(now_time-start_time)
                    #sys.stdout.flush()

                counter = file[(file['chrom']==chrom)&(file['real_start']>=start-stream_num)&(file['real_start']<=start)].index
                RAD.loc[index,'ups_'+histone]+=len(counter)

                counter = file[(file['chrom']==chrom)&(file['real_start']>=start)&(file['real_start']<=end)].index
                RAD.loc[index,'ins_'+histone]+=len(counter)

                counter = file[(file['chrom']==chrom)&(file['real_start']>=end)&(file['real_start']<=end+stream_num)].index
                RAD.loc[index,'downs_'+histone]+=len(counter)

        #proc = os.getpid()
        #print('{0} by process id: {1}'.format(histone, proc))
        if len(fileList)>0:
            RAD.to_csv('../result/cal_histone/final/%s_RAD.csv'%(histone),index=False)
        print '%s is finished' %(histone)
        return RAD

################################################################################################################################
if __name__ == '__main__':
    RAD,histone,histone_fileList=gen_histone_matrix('../data/input_data/peak.csv','../data/raw_data/histone')
    # multiprocesss
    # parrarel
    
    #cores = multiprocessing.cpu_count()

    cores = 8
    pool = multiprocessing.Pool(processes=cores)
    stream_num = 1000
    func = partial(count_histone_parallel_2,RAD,stream_num,histone_fileList)
    result = pool.map(func,histone)
    pool.close()
    pool.join()
    RAD = addToALL(result)
    caL_DNase(RAD,stream_num)
    # serial
    #RAD = count_histone_serial(RAD,1000,histone_fileList,histone)



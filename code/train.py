#!usr/bin/python
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
from gensim.models import Word2Vec
import genVecs
import xgboost as xgb
from xgboost import plot_importance
from matplotlib import pyplot as plt
from sklearn import metrics
from sklearn.cross_validation import cross_val_score,StratifiedKFold,train_test_split
from sklearn import preprocessing
from sklearn.metrics import f1_score,accuracy_score,precision_score,recall_score,roc_auc_score
from sklearn.externals import joblib
from sklearn.model_selection import GridSearchCV,ParameterGrid
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
import warnings
import sys
from puAdapter import PUAdapter

def gen_test(LabelData):
    data = pd.read_csv(LabelData)
    pos = data[data['label']==1]
    pos_test = pos.sample(n=1121,axis=0) 
    pos_train = pos.drop(pos_test.index,axis=0)
    neg = data[data['label']!=1]
    neg_test = neg.sample(n=1121,axis=0) 
    neg_train = neg.drop(neg_test.index,axis=0)

    train = pd.concat([pos_train,neg_train],axis=0)
    train.index = range(len(train))
    train.to_csv('../result/train/train_all.csv',index=False)
    test = pd.concat([pos_test,neg_test],axis=0)
    test.index = range(len(test))
    test.to_csv('../result/train/test_all.csv',index=False)
    return train,test

##############################################################################################################
def SeqToVec(DNAlist1,DNAlist2,num_features):
    model = Word2Vec.load('../model/DNAtoVec.model')
    DNAFeatureVecs = np.zeros((len(DNAlist1),2*num_features), dtype="float32")
    counter = 0   
    for DNA in DNAlist1:
        if counter%100==0:
            print "Processing DNAlist1: %d of %d\r" %(counter,len(DNAlist1)),
            sys.stdout.flush()
        if len(DNA)>1:
            DNAFeatureVecs[counter][0:num_features] = np.mean(model[DNA],axis = 0)
        counter+=1

    counter = 0
    for DNA in DNAlist2:
        if counter%100==0:
            print "Processing DNAlist2: %d of %d\r" %(counter,len(DNAlist2)),
            sys.stdout.flush()
        if len(DNA)>1:
            DNAFeatureVecs[counter][num_features:2*num_features] = np.mean(model[DNA],axis = 0)
        counter+=1

    seqVecData = pd.DataFrame(DNAFeatureVecs)
    return seqVecData

def getDNA_split(DNAdata,word):
    DNAlist1 = []
    DNAlist2 = []
    counter = 0
    for DNA in DNAdata['C1_seq']:
        if counter%1000==0:
            print "Spliting C1_seq: %d of %d\r" %(counter,len(DNAdata)),
            sys.stdout.flush()       
        DNA = str(DNA).upper()
        DNAlist1.append(genVecs.DNA2Sentence(DNA,word).split(' '))
        counter += 1

    counter = 0
    for DNA in DNAdata['C2_seq']:
        if counter%1000==0:
            print "Spliting C1_seq: %d of %d\r" %(counter,len(DNAdata)),
            sys.stdout.flush()   
        DNA = str(DNA).upper()
        DNAlist2.append(genVecs.DNA2Sentence(DNA,word).split(' '))
        counter += 1

    return DNAlist1,DNAlist2

def getSeqFeature(seq_data,word=6,num_features=100):
    print 'preprocessing seqFeature!'
    DNAlist1,DNAlist2 = getDNA_split(seq_data,word)
    seqVecData = SeqToVec(DNAlist1,DNAlist2,num_features)
    print 'getSeqFeature is Ok!'
    return seqVecData

###############################################################################################################
def getLabelData(LabelData):
    data = pd.read_csv(LabelData)
    data['label'][data[data['label']!=1].index]=0
    seq_data = data[['C1_seq','C2_seq']]
    seqVecData = getSeqFeature(seq_data)
    data = data.drop(['chrom','C1_start','C1_end','C1_seq','C2_start','C2_end','C2_seq','num_between','index_back'],axis=1)# 'index_back' in balance
    data = pd.concat([data,seqVecData],axis=1)
    print(data.columns)
    data.to_csv('../result/train/test_datavecs_all.csv',index=False)
    return data

def train(data,LabelData):
    y = data['label']
    dataVecs = data.drop(['label'],axis=1)
    # only histone
    # dataVecs = dataVecs.drop(dataVecs.columns[48:],axis=1)
    # only DNAtoVec
    #dataVecs = dataVecs.drop(dataVecs.columns[0:48],axis=1)
    X = preprocessing.StandardScaler().fit_transform(dataVecs)

    print 'trainging'
    max_depth= 6
    learning_rate= 0.05
    n_estimators = 500
    fold_num = 5
    scale_pos_weight = 0.1
    subsample=0.8
    colsample_bytree = 0.8
    cv = StratifiedKFold(y = y, n_folds = fold_num, shuffle = True, random_state = 0)  


    for train_index, test_index in cv:
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        #joblib.dump(forest,'../model/xgb.dat')
        #forest = joblib.load('../model/xgb.dat')
        forest = xgb.XGBClassifier(max_depth=max_depth,learning_rate=learning_rate, n_estimators = n_estimators,nthread=30) #,scale_pos_weight=scale_pos_weight
        forest.fit(X_train,y_train)

        #traning set
        y_pred = forest.predict(X_train)
        acc = accuracy_score(y_train, y_pred)
        f1 = f1_score(y_train, y_pred, average='micro')
        pre = precision_score(y_train, y_pred, average='micro')
        recall = recall_score(y_train, y_pred, average='micro')
        auc = roc_auc_score(y_train, y_pred)

        fpr, tpr, thresholds = metrics.roc_curve(y_train, y_pred, pos_label=1)
        print(metrics.auc(fpr, tpr))
        print 'acc:%f,f1:%f,pre:%f,recall:%f,auc:%f,auc(pos_label=1):%f\n' %(acc,f1,pre,recall,auc,metrics.auc(fpr, tpr))

        with open('../result/train/train.txt','a') as fr:
            fr.write('%s\n'%(LabelData))
            fr.write('DNAtoVec and histone\n')
            fr.write('max_depth:%d,learning_rate:%f,fold_num:%d,n_estimators:%d,scale_pos_weight:%f\n' %(max_depth,learning_rate,fold_num,n_estimators,scale_pos_weight))
            fr.write('acc:%f,f1:%f,pre:%f,recall:%f,auc:%f\n' %(acc,f1,pre,recall,auc))
        
        #test set
        y_pred = forest.predict(X_test)
        acc = accuracy_score(y_test, y_pred)
        f1 = f1_score(y_test, y_pred, average='micro')
        pre = precision_score(y_test, y_pred, average='micro')
        recall = recall_score(y_test, y_pred, average='micro')
        auc = roc_auc_score(y_test, y_pred)

        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred, pos_label=1)
        print(metrics.auc(fpr, tpr))
        print 'acc:%f,f1:%f,pre:%f,recall:%f,auc:%f,auc(pos_label=1):%f\n' %(acc,f1,pre,recall,auc,metrics.auc(fpr, tpr))
        
        '''
        with open('../result/train/record.txt','a') as fr:
            fr.write('%s\n'%(LabelData))
            fr.write('DNAtoVec and histone\n')
            fr.write('max_depth:%d,learning_rate:%f,fold_num:%d,n_estimators:%d,scale_pos_weight:%f\n' %(max_depth,learning_rate,fold_num,n_estimators,scale_pos_weight))
            fr.write('acc:%f,f1:%f,pre:%f,recall:%f,auc:%f\n' %(acc,f1,pre,recall,auc))
        '''

    #scores= cross_val_score(forest, dataVecs, label,cv = cv,n_jobs = 8)
    #print(scores)

def gridSearchCV(data):
    y_train = data['label'].values
    X_train = data.drop(['label'],axis=1)
    X_train = preprocessing.StandardScaler().fit_transform(X_train)

    cv_params = {'n_estimators': [1100]}
    other_params = {'learning_rate': 0.08,'n_estimators': 1100,'max_depth': 7, 'min_child_weight': 1, 'seed': 0,
                    'subsample': 0.9, 'colsample_bytree': 0.9, 'gamma': 0, 'reg_alpha': 0, 'reg_lambda': 1}
    forest = xgb.XGBClassifier(**other_params)
    #cv_params = {'n_estimators':[300,500,700,900]}
    #other_params = {'n_estimators':500,'min_samples_split':100,'min_samples_leaf':20,'max_depth':8,'max_features':'sqrt'}
    #forest = RandomForestClassifier(**other_params)
    #cv_params = {'kernel': ['rbf','linear','poly']}
    #forest = SVC()
    
    optimized_forest = GridSearchCV(estimator=forest, param_grid=cv_params, scoring='roc_auc', cv=5, verbose=1, n_jobs=8)
    optimized_forest.fit(X_train, y_train)
    evalute_result = optimized_forest.grid_scores_
    print('每轮迭代运行结果:{0}'.format(evalute_result))
    print('参数的最佳取值：{0}'.format(optimized_forest.best_params_))
    print('最佳模型得分:{0}'.format(optimized_forest.best_score_))

def test(train_data,test_data):
    y_train = train_data['label'].values
    X_train = train_data.drop(['label'],axis=1)
    X_train = preprocessing.StandardScaler().fit_transform(X_train)

    y_test = test_data['label'].values
    X_test = test_data.drop(['label'],axis=1)
    X_test = preprocessing.StandardScaler().fit_transform(X_test)

    other_params = {'learning_rate': 0.08,'n_estimators': 1100,'max_depth': 7, 'min_child_weight': 1, 'seed': 0,
                    'subsample': 0.9, 'colsample_bytree': 0.9, 'gamma': 0, 'reg_alpha': 0, 'reg_lambda': 1}
    forest = xgb.XGBClassifier(**other_params)
    #forest = PUAdapter(forest)
    forest.fit(X_train,y_train)

    #traning set
    y_pred = forest.predict(X_train)
    acc = accuracy_score(y_train, y_pred)
    f1 = f1_score(y_train, y_pred, average='micro')
    pre = precision_score(y_train, y_pred, average='micro')
    recall = recall_score(y_train, y_pred, average='micro')
    auc = roc_auc_score(y_train, y_pred)

    fpr, tpr, thresholds = metrics.roc_curve(y_train, y_pred, pos_label=1)
    print(metrics.auc(fpr, tpr))
    print 'acc:%f,f1:%f,pre:%f,recall:%f,auc:%f,auc(pos_label=1):%f\n' %(acc,f1,pre,recall,auc,metrics.auc(fpr, tpr))

    y_pred = forest.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred, average='micro')
    pre = precision_score(y_test, y_pred, average='micro')
    recall = recall_score(y_test, y_pred, average='micro')
    auc = roc_auc_score(y_test, y_pred)

    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred, pos_label=1)
    print(metrics.auc(fpr, tpr))
    print 'acc:%f,f1:%f,pre:%f,recall:%f,auc:%f,auc(pos_label=1):%f\n' %(acc,f1,pre,recall,auc,metrics.auc(fpr, tpr))
    
    
    with open('../result/test/record.txt','a') as fr:
        for k,v in other_params.items():
            fr.write(k+str(v)+'\t')
        fr.write('\n')
        fr.write('acc:%f,f1:%f,pre:%f,recall:%f,auc:%f\n' %(acc,f1,pre,recall,auc))

if __name__ == '__main__':
    word = 6
    num_features = 100
    LabelData = '../data/input_data/LabelData_neg3_balance.csv'
    train_data = '../result/train/train_datavecs_all.csv'
    test_data = '../result/train/test_all.csv'
    warnings.filterwarnings('ignore')
    #gen_test(LabelData)
    train_data = pd.read_csv('../result/train/train_datavecs_all.csv')
    #data = getLabelData(test_data)
    #train(data,LabelData)
    #gridSearchCV(train_data)
    test_data = pd.read_csv('../result/train/test_datavecs_all.csv')
    test(train_data,test_data)
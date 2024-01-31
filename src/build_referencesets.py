import sys
import math
import random
import numpy as np
import pandas as pd
import scipy.stats as stats

from scipy.special import betaln
from scipy.optimize import minimize


def binomial_likelihood(mean,data1,data2,exons):
    lm = np.log(mean)
    lm1 = np.log(1.0-mean)
    return -sum([data1[i]*lm + data2[i]*lm1 for i in exons])


def neg_log_likelihood(params,data1,data2,exons):
    a, b = params
    ll = 0
    c = betaln(a,b)
    for i in exons:
        ll += betaln(data1[i]+a, data2[i]+b)-c
    return -ll

def filter_lmexons(allexons_df,file1=None):
    if file1 == None: 
        return allexons_df
    filtered_exons = pd.read_csv(file1,sep='\t',header=None)
    exon_list = filtered_exons[filtered_exons.columns[3]].tolist()
    exon_index = {}

    for index in exon_list: 
        exon_index[index] = 1
    rows_to_remove = []

    for i in range(allexons_df.shape[0]):
        if allexons_df.iloc[i,3] not in exon_index and allexons_df.iloc[i,4] not in exon_index:            
            rows_to_remove.append(i)
    
    mask = ~allexons_df.index.isin(rows_to_remove)
    filtered_df=allexons_df[mask]

    return filtered_df
        

## df is counts data table for all exons from ExomeDepth
def sample_readsums(count_matrix_file,exons_for_beta=10000,max_ref=15):
    allexons_df = pd.read_csv(count_matrix_file,sep='\t')
    
    for i in range(allexons_df.shape[1]):
        if 'bam' in allexons_df.columns[i]: 
            first = i
            break
    print(allexons_df.columns[0:10],allexons_df.shape,first,allexons_df.shape[1]-first,file=sys.stderr)

    filtered_df = filter_lmexons(allexons_df)
    print(allexons_df.shape,'filtered down to',filtered_df.shape,file=sys.stderr)
        
    exons = random.sample([i for i in range(filtered_df.shape[0])],min(exons_for_beta,filtered_df.shape[0]))
    df1 = filtered_df.drop(columns=filtered_df.columns[0:first])
    
    corrmat = df1.corr()
    for i in range(df1.shape[1]): 
        print('index',i,df1.columns[i],df1[df1.columns[i]].mean(axis=0),','.join([str(round(a,3)) for a in corrmat.iloc[i,:].tolist()]))

    maxpairs =min(max_ref,allexons_df.shape[1])

    for i in range(df1.shape[1]):
        vec1 = df1.iloc[:,i].tolist()
        corr_list = list(enumerate(corrmat.iloc[i,:].tolist()))
        corr_list.sort(reverse=True,key=lambda x:x[1])
        #print(i,corr_list[1:1+maxpairs],file=sys.stderr)

        for j,corr in corr_list[1:1+maxpairs]:
            vec2 = df1.iloc[:,j].tolist()
            s1 = sum(vec1); s2 = sum(vec2)
            p = s1/(s1+s2)

            mean = p
            best_results = [0,-1,-1]

            left = 1
            right = 1000000
            iters = 0

            while iters < 100:
                mid = (left+right)/2
                ll_left = neg_log_likelihood((mean*left,left-mean*left),vec1,vec2,exons)
                ll_mid = neg_log_likelihood((mean*mid,mid-mean*mid),vec1,vec2,exons)
                ll_right = neg_log_likelihood((mean*right,right-mean*right),vec1,vec2,exons)

                if abs(ll_mid-ll_left)/ll_mid < 1e-7 and iters > 2: 
                    best_results = [ll_mid,mean*mid,mid-mid*mean]
                    break
                #print(ll_left,ll_right,ll_mid,left,right,mid)
                if ll_right < ll_left: 
                    left = mid
                else: 
                    right = mid
                iters += 1

            binll = binomial_likelihood(p,vec1,vec2,exons)
            #print('binomial',binll,simll,binll-simll)
            betabinll = best_results[0]
            delta = binll-betabinll

            if delta < 30: 
                print('best',best_results[0],best_results[1],best_results[2],i,j,corr,'fit',binll,betabinll,delta,mean)
            else: 
                print('BB',best_results[0],best_results[1],best_results[2],i,j,corr,'fit',binll,betabinll,delta,mean)

sample_readsums(sys.argv[1])


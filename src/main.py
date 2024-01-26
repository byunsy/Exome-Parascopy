import pandas as pd
import numpy as np
import os,sys,gzip
import scipy.optimize as opt
#import scipy.misc as misc
#from scipy.stats import spearmanr,kendalltau
import math
import random
import networkx as nx
from exomecounts import *
from optimize_functions import *

import argparse

MAX_PENALTY=1000

def parse_args():
    # Parse arguments from command-line
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="read count file for indiviidual gene")
    parser.add_argument("-o", "--outfile", required=True,help="path to output file")
    parser.add_argument("-p", "--params", required=True,      help="path to parameters file")
    parser.add_argument("--truth", required=False,       help="path to file with true copy number values",default=None)
    parser.add_argument("-c", "--refcn", default=4,         help="reference copy number")
    parser.add_argument("--max", default=1000,         help="maximum value of penalty scale parameter")
    args_parsed = parser.parse_args()
    return args_parsed

    
def local_updates(x,exome_data,update=True,DEBUG=False,order='optimized'):
    iterations=0

    L = list(range(exome_data.n))
    if order =='optimized' and exome_data.order != None: L = exome_data.order

    total_updates=0
    current_ll = fopt(x,exome_data)
    initial_ll = current_ll
    while iterations < 20:
        updates=0
        if order =='randomize': random.shuffle(L) ## every iteration..
        for i in L:
            current = x[i]
            ll_list=[]
            for c in range(1,10):
                x[i] = c
                ll_list.append((fopt(x,exome_data)-current_ll,c))
            x[i] = current
            ll_list.sort()
            #print('iter',iterations,i,ll_list[0],ll_list[1],'debug')
            if ll_list[0][1] != current and update: 
                x[i] = ll_list[0][1]
                updates +=1
                current_ll = fopt(x,exome_data)
                if DEBUG: print('updating',i,'iter',iterations,ll_list[0],current_ll,'from',current,exome_data.trueCN[i])
            elif not update:
                print(i,ll_list[0],current_ll,'from',current,exome_data.trueCN[i])
        iterations +=1 
        total_updates += updates
        if updates==0: break
    if DEBUG: print('CN local updates',initial_ll,current_ll,-1*current_ll+initial_ll,'iter',iterations,'updated',total_updates)
    return initial_ll-current_ll,iterations,total_updates,current_ll


def penalized_likelihood(self,result_start,step=20,control=None):
    ll_nopenalty=result_start.fun
    x0 = [self.refCN]*self.n
    for i in range(self.n): x0[i] = result_start.x[i]

    prev = result_start.fun
    pcnvec = result_start.x*(self.refCN/np.mean(result_start.x))
    flag=0
    for ld in range(20,MAX_PENALTY,step):
        self.Lambda = ld
        result = opt.minimize(fopt,x0,bounds=self.cn_bounds,tol=1e-8,args=(self,),method='L-BFGS-B') # Nelder-Mead
        #rankcc = spearmanr(np.around(result.x,2),self.trueCN)

        meancn = np.mean(result.x[0:self.n]); 
        scalefactor= self.refCN/meancn
        cnvec = result.x*scalefactor
        varcn = np.var(result.x[0:self.n]); 
        sum_delta = sum([abs(cnvec[i]-pcnvec[i]) for i in range(self.n)])
        pcnvec=cnvec

        roundedvec = [round(a*scalefactor,2) for a in sorted(result.x)]
        diffvalues=1
        for i in range(self.n-1):
            if roundedvec[i+1]-roundedvec[i] >0.2: diffvalues +=1
        print(roundedvec)
        #print('diff',diffvalues,roundedvec,self.n,diffvalues*math.log(self.n),file=sys.stderr)
        self.Lambda=0; ll = fopt(result.x[0:self.n],self)
        #print('ll penalty',result.fun,ld,ll-ll_nopenalty,result.fun-ll,diffvalues,'meancn',meancn,varcn,file=sys.stderr)
        print('ll penalty',result.fun,ld,ll-ll_nopenalty,result.fun-ll,diffvalues,'meancn',meancn,varcn,file=sys.stdout)
        if abs(prev-result.fun)/(step*result.fun) < 1e-7: 
            print('break loop',abs(prev-result.fun)/(step*result.fun))
            flag +=1
            if flag ==2: break
        prev =result.fun
        for i in range(self.n): x0[i] = result.x[i]
    return result,ld


def convert_integerCN(self,result):
    counts ={}
    meancn = np.mean(result.x[0:self.n]); 
    scalefactor= self.refCN/meancn
    for i in range(self.n): result.x[i] *= scalefactor
    for a in result.x:
        try: counts[round(a,1)] +=1
        except KeyError: counts[round(a,1)] = 1
    clusters = sorted([(cn,count) for cn,count in counts.items()],reverse=False)
    bestpair = [0,0]
    for pair in clusters:
        if pair[1] > bestpair[1]: bestpair =pair
    print('mean',meancn,scalefactor,clusters,bestpair)

    for CN in [0,self.refCN-1,self.refCN,self.refCN+1]:
        scalefactor = CN/bestpair[0]
        newcn1 = [round(result.x[i]*scalefactor,3) for i in range(self.n)]
        newcn= [int(round(a,0)) for a in newcn1]
        if CN ==0: newcn = [self.refCN for i in range(self.n)]
        delta,a,b,current = local_updates(newcn,self,update=True,DEBUG=False,order='optimized')
        print(CN,'local',delta,a,b,current,fopt(newcn,self),newcn)
        print(CN,'error',self.n,sum([ math.ceil(abs(newcn[i]-self.trueCN[i])/10) for i in range(self.n)]))



def best_fractional(self,minCN=1.01,maxCN=10.01,control=None):
    self.cn_bounds = [(minCN,maxCN) for i in range(self.n)]

    self.graph = pairwise_graph(self)

    x0 = [self.refCN]*self.n
    result0 = opt.minimize(fopt,x0,bounds=self.cn_bounds,tol=1e-8,args=(self,),method='L-BFGS-B') # Nelder-Mead
    print('no-penalty',round(result0.fun,3),'\n','CN-vector',' '.join([str(a) for a in sorted(result0.x)]),'\n')
    ll_nopenalty= result0.fun
    
    self.compare = sorted([[result0.x[i],self.trueCN[i],i] for i in range(self.n)])
    for i in range(self.n):  print(self.compare[i],end=' ')
    print()
    self.order= [self.compare[i][2] for i in range(self.n)]

    #######
 
    if MAX_PENALTY> 0:    
        meancn = np.mean(result0.x)
        result0.x *= self.refCN/meancn
        """
        for i in range(self.n): 
            if result0.x[i] < 1: self.cn_bounds[i] = (0.01,2)
            else: self.cn_bounds[i] = (1,maxCN)
        """
        result,ldval = penalized_likelihood(self,result0)
    else:    result = result0

    
    meancn = np.mean(result.x[0:self.n]); varcn = np.var(result.x[0:self.n])
    scalefactor= self.refCN/meancn
    for i in range(self.n): result.x[i] *= scalefactor
    print('mean',meancn,'var',varcn,scalefactor)

    self.bestcnvec = np.around(result.x[0:self.n],3)
    print('contCN maxll',round(result.fun,3),'\n','CN-vector',' '.join([str(a) for a in sorted(self.bestcnvec)]),'\n')
    self.bestLL = result.fun
    print('contCN maxll',round(result.fun,3),self.n,file=sys.stderr)#,np.around(result.x/low,2))

    convert_integerCN(self,result)

    #diagnostic(result.x,self,self.trueCN)
    self.Lambda=0
    return result,result0


def analyze_component(data1,refCN,outfile=sys.stdout):
    print('\n\nanalyzing connected component',c,comp,data1.n)
    data1.refCN = refCN
    data1.cn = [data1.refCN]*data1.n

    result,result0 = best_fractional(data1)

    data1.cn = [refCN]*data1.n
    delta,a,b,current = local_updates(data1.cn,data1,update=True,DEBUG=False,order='optimized')

    # data1.trueCN[i] is zero if missing
    errors = sum([ min(data1.trueCN[i],math.ceil(abs(data1.cn[i]-data1.trueCN[i])/10)) for i in range(data1.n)])
    errors1 = sum([ min(data1.trueCN[i],math.ceil(abs(data1.cn[i]-data1.trueCN[i]-1)/10)) for i in range(data1.n)]) ## all samples +1
    print('stats',delta,current,fopt(data1.cn,data1),'errors',min(errors,errors1),data1.n)
    print('best',''.join([str(a) for a in data1.cn]))
    print('trut',''.join([str(a) for a in data1.trueCN]))
    print('diff',''.join([str(abs(data1.cn[i]-data1.trueCN[i])) for i in range(data1.n)]))
    errors2=min(errors,errors1)

    print('#component',c,'size',data1.n,'likelihood',result.fun,file=outfile)
    for i in range(data1.n):  print('final-CN',c,data1.n,data1.samplenames[i],data1.cn[i],result.x[i],result0.x[i],data1.trueCN[i],file=outfile)

    fval0 = fopt(data1.cn,data1); 
    fval1 = fopt(data1.trueCN,data1); 
    print('----------testing trueCN likelihood ------------',fval0,fval1,'\n')

    return errors2

################################################################################################################

args = parse_args()

MAX_PENALTY=int(args.max)
refCN = int(args.refcn)
data =ExomeData()
data.get_parameters(args.params)
data.gene_counts(args.input)
if args.truth != None: data.read_trueCN(args.truth)
else: data.trueCN = None
outfile=open(args.outfile,'w')

stats = [0,0,0]
c=0
for comp in data.comp_list:
    data_cc = data.construct_subset(c)
    if comp <= 2: 
        print('small component',file=sys.stderr)
        print('#component',c,'size',data1.n,'likelihood',-1,file=outfile)
        for i in range(data1.n):  print('final-CN',c,data1.n,data1.samplenames[i],'-1',file=outfile)
        stats[2] += comp
        continue
    errors = analyze_component(data_cc,refCN,outfile=outfile)
    stats[0] += errors; stats[1] += data_cc.n
    c+=1

print('final accuracy stats',stats[0],stats[2],stats[1])
outfile.close()



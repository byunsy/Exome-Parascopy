import os
import sys
import math
import random
import numpy as np
import pandas as pd
import networkx as nx
import scipy.optimize as opt

class ExomeData:
    calls = 0

    def __init__(self,n=0):
        self.samples_index = {}
        self.samplenames = [] ## list
        self.n = n ## number of samples
        self.means = []
        self.counts = []
        self.ecounts = []
        self.counts_scaled = []
        self.reference_sets = [] 
        self.ref_sums = []
        self.components = {}
        self.comp_list = []
        self.betamatrix = {}
        self.correlations = []
        self.graph = {}
        self.Lambda = 0
        self.cn_bounds = []
        self.cn = None
        self.fracCN = None
        self.trueCN = None
        self.minCN = 0
        self.maxCN = 10
        self.order = None
        self.bestLL =- 100000
        self.bestLL_discrete =- 1
        self.compare = None
        self.bestcnvec = None
        self.map = []
        self.prior = {}

    def construct_subset(self,comp): ## create new ExomeData object for a connected component
        subdata = ExomeData()
        subdata.trueCN = []
        subdata.betamatrix = {}
        subdata.correlations = []
        subdata.samplenames=[]
        subdata.prior = self.prior
        subrows = []
        temp_map = [0]*self.n
        
        for i in range(self.n):
            if i not in self.components or self.components[i] != comp: 
                continue
            temp_map[i] = subdata.n
            subrows.append(i)
            subdata.samplenames.append(self.samplenames[i])
            #print('map',i,subdata.n,end=' ')
            subdata.n += 1
        if len(self.correlations) > 0 : 
            for i in range(subdata.n): 
                subdata.correlations.append([ self.correlations[subrows[i]][j] for j in subrows ])
        ###print('subdata')
        for key,value in self.betamatrix.items():
            newkey = (temp_map[key[0]],temp_map[key[1]])
            subdata.betamatrix[newkey] = value

        for i in range(self.n):
            if i not in self.components or self.components[i] != comp: 
                continue
            subdata.means.append(self.means[i])
            subdata.counts.append(self.counts[i])
            subdata.ecounts.append(self.ecounts[i])
            subdata.reference_sets.append([temp_map[s] for s in self.reference_sets[i]])
            subdata.trueCN.append(self.trueCN[i])
        subdata.map = [0]*subdata.n
        
        for i in range(subdata.n): 
            subdata.map[i] = i
        subdata.bestLL_discrete =- 1
        subdata.ref_sums = [sum([subdata.counts[j] for j in subdata.reference_sets[i]]) for i in range(subdata.n)]
        
        return subdata
    
    def connected_comp(self):
        G = nx.Graph()
        for i in range(self.n):
            G.add_edges_from([(i,j) for j in self.reference_sets[i]])
        components = list(nx.connected_components(G))
        i = 0
        #for i in range(self.n): self.components[i] = -1
        print("\nCONNECTED COMPONENTS (undirected)")
        for comp in components:             
            list_nodes = sorted(comp)
            self.comp_list.append(len(list_nodes))
            for v in list_nodes: 
                self.components[v] = i
            print(f'conn-comp-{i}: {list_nodes}')
            i += 1
    
    def connected_comp_directed(self,directed=True):
        G = nx.DiGraph()
        for i in range(self.n):
            G.add_edges_from([(j,i) for j in self.reference_sets[i]])
        components = list(nx.strongly_connected_components(G))
        i = 0
        #for i in range(self.n): self.components[i] = -1
        print("\nCONNECTED COMPONENTS (directed)")
        for comp in components:             
            list_nodes = sorted(comp)
            self.comp_list.append(len(list_nodes))
            for v in list_nodes: 
                self.components[v] = i
            print(f'conn-comp-{i}: {list_nodes}')
            i += 1
    
    def get_parameters(self, pfile, prior_file=None, min_sum=20000): ## read the betafit.out file
        self.means = []
        self.n = 0
        with open(pfile, 'r') as file1:
            for line in file1:
                if line.startswith('index'): 
                    v = line.strip().split()
                    self.samples_index[v[2]] = int(v[1])
                    simple_id = v[2].split('.')[0] ## HG00100.mapped.ILLUMINA.exome.bam -> HG00110
                    self.samples_index[(simple_id,0)] = int(v[1])
                    self.samplenames.append(simple_id)
                    self.n += 1
                    self.means.append(float(v[3]))
                    if len(v) > 4: 
                        self.correlations.append([float(a) for a in v[4].split(',')])
                        #print(self.correlations[-1])
                elif line.startswith('best') or line.startswith('BB'):
                    v = line.strip().split()
                    s1 = int(v[4])
                    s2 = int(v[5])
                    alpha = float(v[2])
                    beta = float(v[3])
                    corr = float(v[6])
                    if alpha + beta > min_sum: 
                        self.betamatrix[(s1,s2)] = (alpha,beta)

        for i in range(self.n): 
            self.reference_sets.append([])

        for key,value in self.betamatrix.items(): 
            self.reference_sets[key[0]].append(key[1])
        # for i in range(self.n): self.reference_sets[i] = self.reference_sets[i][0:4] ## sets reference-set size to 4

        self.connected_comp()

        # Get prior probs if provided
        if prior_file:
            print("\nPRIOR PROBABILITY")
            print("- Found prior probabilities:", prior_file)
            prior = pd.read_csv(prior_file, sep="\t")
            self.prior = prior.set_index('cn')['prob'].to_dict()

    def gene_counts(self,cfile):
        df = pd.read_csv(cfile,sep='\t')
        self.counts = [0]*self.n
        self.ecounts = [[] for i in range(self.n)]
        ## pandas ignores the first column with the 1,2,3... indexes | no problem
        #print(df.columns[4:10])
        #indices = [i for i, item in enumerate(df.columns) if 'bam' in item]

        matched = 0
        missing = []
        ## assumption that fifth column is the bam file names
        if 'bam' in df.columns[0] or 'bam' in df.columns[3]:
            print('Missing exon columns',file=sys.stderr)
            sys.exit()
        for i in range(4,df.shape[1]): 
            sample = df.columns[i].strip('\"')
            #print(sample,sample.split('.')[0])
            sample_index = -1
            try:
                sample_index = self.samples_index[sample]
            except KeyError: 
                try: 
                    sample_index = self.samples_index[(sample.split('.')[0],0)]
                except KeyError: 
                    pass
            if sample_index != -1:
                count_sum = df[df.columns[i]].sum(axis=0)
                #print(sample,count_sum,sample_index,df[df.columns[i]].tolist())
                self.counts[sample_index] = count_sum
                self.ecounts[sample_index] = df[df.columns[i]].tolist()
                matched += 1
            else:
                missing.append(sample)
        print("\nGENE COUNTS")
        print('- Matched:', matched)
        print('- Missing:', len(missing))
        ###print(self.counts)
        if matched < self.n: 
            print('Number of samples with count data is less than the number in the statistics file', matched, self.n, file=sys.stderr)
            print('Check input files', file=sys.stderr)
            sys.exit()
       
    def read_trueCN(self,filename): 
        missing = 0
        try: 
            self.trueCN = [0]*self.n
            with open(filename, 'r') as file1:
                for line in file1: 
                    sample = line.split()[0]
                    CN = int(line.split()[1])
                    try: 
                        sample_index = self.samples_index[(sample,0)]
                        self.trueCN[sample_index] = CN
                    except KeyError:
                        missing += 1
            print("\nTRUE-CN")
            print('- TrueCN matched:', len(self.trueCN))
            #print(self.trueCN)
            print('- TrueCN missing:', missing)
        except FileNotFoundError: self.trueCN = None



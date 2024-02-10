import sys
import math

def fopt(x, exome_data):
    exome_data.calls += 1
    total_ll = 0.0
    vec = [x[i]*exome_data.means[i] for i in range(exome_data.n)]
    penalty = 0

    for i in range(exome_data.n):
        try: 
            edgelist = sorted(exome_data.graph[i])
            pterm = sum([ value*(x[i]-x[j])*(x[i]-x[j]) for (value,j) in edgelist[0:15]])
            #pterm =sum([ value*math.sqrt( (x[i]-x[j])*(x[i]-x[j]) ) for (value,j) in edgelist[0:15]])
            penalty += pterm
        except KeyError: 
            pass
     
    for i in range(exome_data.n):
        try: 
            p1 = vec[i] + sum([vec[j] for j in exome_data.reference_sets[i]])
            a = math.log(p1)
            p = vec[i] / p1
            prior_prob = 1
            if exome_data.prior:
                prior_prob = exome_data.prior.get(round(vec[i]), 1e-8) ### testing prior
            ll = math.log(p)*exome_data.counts[i] + math.log(1.0-p)*exome_data.ref_sums[i] + math.log(prior_prob) ### testing prior
        except ValueError: 
            ll = 0

        total_ll += ll

    #if exome_data.n < len(x): return -1*total_ll+ x[exome_data.n]*penalty
    return -1*total_ll + exome_data.Lambda*math.sqrt(penalty) #math.sqrt(penalty)

    #return -1*total_ll + exome_data.Lambda*penalty
    #print('LL', total_ll, penalty)


def pairwise_graph(data, thresh=5): ## calculate likelihoods for each edge
    graph = {}
    edges = []

    for i in range(data.n): 
        edges.append(0)

    for i in range(data.n):
        for j in data.reference_sets[i]:
            p = data.means[i] / (data.means[i] + data.means[j])
            best_triple = [-1e8, -1e8, -1e8]
            for cn1 in range(1, data.refCN+3):
                for cn2 in range(1, data.refCN+3):
                    p = cn1*data.means[i] / (cn1*data.means[i] + cn2*data.means[j])
                    ll = math.log(p)*data.counts[i] + math.log(1.0-p)*data.counts[j]
                    if cn1 == cn2 and ll > best_triple[0]: 
                        best_triple[0] = ll
                    elif cn1 < cn2 and ll > best_triple[1]: 
                        best_triple[1] = ll
                    elif cn1 > cn2 and ll > best_triple[2]: 
                        best_triple[2] = ll
            d1 = best_triple[0] - best_triple[1]
            d2 = best_triple[0] - best_triple[2]

            if d1 >= thresh and d2 >= thresh:     
                #print('edge',i,j,data.correlations[i][j])
                try: 
                    graph[i].append((min(d1,d2),j))
                except KeyError: 
                    graph[i] = [(min(d1,d2),j)]
                edges[i] += 1
            else: 
                d1 = best_triple[2] - best_triple[0]
                d2 = best_triple[2] - best_triple[1]
                #if d1 >= thresh*2 or d2 >= thresh*2: graph[(i,j)] = (min(d1,d2),'gt')
        #print('node',i,'CN',data.trueCN[i],data.reference_sets[i],edges)
    ###print('noedges', edges.count(0), edges)
    return graph


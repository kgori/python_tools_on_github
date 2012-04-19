#!/usr/bin/env python

def pairings(l):
    import numpy as np
    m = []
    for i in range(len(l)):
        m.append([])
        for j in range(len(l)):
            if i==j: m[i].append(0)
            elif l[i]==l[j]: m[i].append(1)
            elif l[i]!=l[j]: m[i].append(-1)
    return np.array(m)

def comp2(t,f,ttp=1,tp=1,tn=1,fn=-1,fp=-1):
    import numpy as np
    m = []
    for row in range(len(t)):
        m.append([])
        for col in range(len(t)):
            if f[row][col]==t[row][col]:
                val = t[row][col]
                if val==0: m[row].append(ttp) #trivial tp - self always pairs with self
                elif val==1: m[row].append(tp) #tp - pairing exists in t and in f
                else: m[row].append(tn) #tn - pairing doesn't exist in t or in f
            else:
                val = t[row][col]
                if val == 1: m[row].append(fn) #fn - pairing exists in t, not in f
                elif val == -1: m[row].append(fp) #fp - pairing pairing exists in f, not in t
    return np.array(m)

def w(t,f):
    f = pairings(f)
    t = pairings(t)
    tp = comp2(t,f,0,1,0,0,0) # N11
    tn = comp2(t,f,0,0,1,0,0) # N00
    fn = comp2(t,f,0,0,0,1,0) # N10
    fp = comp2(t,f,0,0,0,0,1) # N01
    l = [0.5*sum(sum(x)) for x in [tp,tn,fn,fp]]
    return l

def make_clustering(l):
    clusters = list(set(l))
    result = []
    for c in clusters:
        cluster_k = []
        for i in range(len(l)):
            if c==l[i]:
                cluster_k.append(i)
        result.append(set(cluster_k))
    return result

def confusion_matrix(t,f): # where t and f are lists describing two clusterings of the same data
    import numpy as np
    clustering_t = make_clustering(t)
    clustering_f = make_clustering(f)
    result = []
    for i in range(len(clustering_t)):
        row = []
        for j in range(len(clustering_f)):
            row.append(len(clustering_t[i].intersection( clustering_f[j])))
        result.append(row)
    return np.array(result)

def probability_distribution(l):
    cl = make_clustering(l)
    total = float(len(l))
    return [ len(x)/total for x in cl ]

def entropy(l):
    from math import log
    pd = probability_distribution(l)
    return abs(sum([x*log(x,2) for x in pd]))

def joint_distribution(t,f):
    cm = confusion_matrix(t,f)
    total = float(len(t))
    return cm/total

def mutual_information2(t,f):
    """ shame this one doesn't always work, because it looks neater
        P.S. This one doesn't always work (but in those cases, jd = jd.T fixes it"""
    from math import log
    s = 0 
    jd = joint_distribution(t,f)
    pdt = probability_distribution(t)
    pdf = probability_distribution(f)
    lt=len(pdt)
    lf=len(pdf)
    for i in range(lt):
        for j in range(lf):
            if jd[i][j]==0: continue # because 0 * log(0) = 0 (lim x->0: xlog(x)->0)
            s+=jd[i][j]*log(jd[i][j]/(pdt[i]*pdf[j]),2)
            print i,j,jd[i][j]/(pdt[i]*pdf[j])
    return s

def mutual_information(t,f):
    # Note: Proof that 0*log(0)=0, by L'hopital's rule on limits:
    # lim (x->c) f(x)/g(x) = lim(x->c) f'(x)/g'(x)
    # Rewrite xlog(x) as log(x) / 1/x
    # if f(x) = log(x) and g(x) = 1/x, then xlog(x) = f(x)/g(x)
    # Also f'(x) = 1/x and g'(x) = -1/x^2 
    # In the limit (x->0): xlog(x) = 1/x / -1/x^2 = -x, and -x -> 0
    from math import log
    ct = make_clustering(t)
    cf = make_clustering(f)
    lt = len(ct)
    lf = len(cf)
    m = len(t)
    s = 0
    for i in range(lt):
        for j in range(lf):
            intersect = float(len(ct[i].intersection(cf[j])))
            if intersect == 0: 
                continue # because 0 * log(0) = 0 
            else: 
                s += (intersect/m)*log(m*intersect/(len(ct[i])*len(cf[j])),2)
    return s

def variation_of_information(t,f): 
    return entropy(t)+entropy(f)-2*mutual_information(t,f)

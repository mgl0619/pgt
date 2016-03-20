# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 22:38:50 2014

@author: gangliu
"""

std_aa_mass = {
    'G': 57.02146,
    'A': 71.03711,
    'S': 87.03203,
    'P': 97.05276,
    'V': 99.06841,
    'T': 101.04768,
    'C': 103.00919,
    'L': 113.08406,
    'I': 113.08406,
    'N': 114.04293,
    'D': 115.02694,
    'Q': 128.05858,
    'K': 128.09496,
    'E': 129.04259,
    'M': 131.04049,
    'H': 137.05891,
    'F': 147.06841,
    'R': 156.10111,
    'Y': 163.06333,
    'W': 186.07931,
    }
"""A dictionary with monoisotopic masses of the twenty standard
amino acid residues.
"""

glycanresidue_mass = {
    'h': 162.0528235,
    'n': 203.0793724,
    's': 291.0954164,
    'g': 307.0903311,
    'f': 146.0579089,
    'x': 132.0422588,
    'z': 79.95681459,
    'p': 79.96633093,
    'u': 176.0320881,
    'k': 250.0688675,
    'q': 292.07289454
    }
    
"""A dictionary with monoisotopic masses of the common monosaccharide 
residues.
"""

def staffordRandFixedSum(n, u, nsets):
    import numpy

    #deal with n=1 case
    if n == 1:
        return numpy.tile(numpy.array([u]),[nsets,1])

    k = numpy.floor(u)
    s = u
    step = 1 if k < (k-n+1) else -1
    s1 = s - numpy.arange( k, (k-n+1)+step, step )
    step = 1 if (k+n) < (k-n+1) else -1
    s2 = numpy.arange( (k+n), (k+1)+step, step ) - s

    tiny = numpy.finfo(float).tiny
    huge = numpy.finfo(float).max

    w = numpy.zeros((n, n+1))
    w[0,1] = huge
    t = numpy.zeros((n-1,n))

    for i in numpy.arange(2, (n+1)):
        tmp1 = w[i-2, numpy.arange(1,(i+1))] * s1[numpy.arange(0,i)]/float(i)
        tmp2 = w[i-2, numpy.arange(0,i)] * s2[numpy.arange((n-i),n)]/float(i)
        w[i-1, numpy.arange(1,(i+1))] = tmp1 + tmp2;
        tmp3 = w[i-1, numpy.arange(1,(i+1))] + tiny;
        tmp4 = numpy.array( (s2[numpy.arange((n-i),n)] > s1[numpy.arange(0,i)]) )
        t[i-2, numpy.arange(0,i)] = (tmp2 / tmp3) * tmp4 + (1 - tmp1/tmp3) * (numpy.logical_not(tmp4))

    m = nsets
    x = numpy.zeros((n,m))
    rt = numpy.random.uniform(size=(n-1,m)) #rand simplex type
    rs = numpy.random.uniform(size=(n-1,m)) #rand position in simplex
    s = numpy.repeat(s, m);
    j = numpy.repeat(int(k+1), m);
    sm = numpy.repeat(0, m);
    pr = numpy.repeat(1, m);

    for i in numpy.arange(n-1,0,-1): #iterate through dimensions
        e = ( rt[(n-i)-1,...] <= t[i-1,j-1] ) #decide which direction to move in this dimension (1 or 0)
        sx = rs[(n-i)-1,...] ** (1/float(i)) #next simplex coord
        sm = sm + (1-sx) * pr * s/float(i+1)
        pr = sx * pr
        x[(n-i)-1,...] = sm + pr * e
        s = s - e
        j = j - e #change transition table column if required

    x[n-1,...] = sm + pr * s

    #iterated in fixed dimension order but needs to be randomised
    #permute x row order within each column
    for i in xrange(0,m):
        x[...,i] = x[numpy.random.permutation(n),i]

    return numpy.transpose(x);
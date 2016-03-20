# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 22:38:50 2014

@author: gangliu
"""
class Glycan:
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
    
    glycanresidue_formula = {
    'h': {'C':6,'H':10,'O':5},
    'n': {'C':8,'H':13,'O':5,'N':1},
    's': {'C':11,'H':17,'O':8,'N':1},
    'g': {'C':11,'H':17,'O':9,'N':1},
    'f': {'C':6,'H':10,'O':4},
    'x': {'C':5,'H':8,'O':4},
    'z': {'S':1,'O':3},
    'p': {'P':1,'O':3,'H':1},
    'u': {'C':6,'H':10,'O':5},
    'k': {'C':6,'H':8,'O':6},
    'q': {'C':10,'H':16,'O':6,'N':2,'S':1}
    }
    
    def __init__(self,sgpstr):
        self.sgpstr = sgpstr
    def printmass(self,monolet):
        print Glycan.glycanresidue_formula(monolet)
print Glycan.glycanresidue_mass['h']   
print Glycan.glycanresidue_formula['q']   

        
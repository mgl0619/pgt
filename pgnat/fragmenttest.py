# -*- coding: utf-8 -*-
"""
Created on Wed Oct 15 22:30:11 2014

@author: gangliu
"""

#! /usr/bin/python

import unittest
from fragment import splitGlyPep

class splitTest(unittest.TestCase):
    
    def setUp(self):
        print "\nstart test"
    
    def tearDown(self):
        print "\nend test"

    def test1(self): 
        sgptest1 = splitGlyPep('GYM<o>KNCT<s>')
        self.assertEqual(sgptest1.peptide,'GYMKNCT')
        nonglycanptm=['o', 's']
        self.assertEqual(sgptest1.nonglycanptm,nonglycanptm)
        glycanptm=[]; nonglycanpos=[3, 7]; glycanpos=[]
        self.assertEqual(sgptest1.glycanptm,glycanptm)
        self.assertEqual(sgptest1.nonglycanpos,nonglycanpos)
        self.assertEqual(sgptest1.glycanpos,glycanpos)
    
    def test2(self):
        sgptest2 = splitGlyPep('GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT{n{h{s}}{n{h{s}{f}}}}R')
        glycanptm = ['{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}','{n{h{s}}{n{h{s}{f}}}}']
        glycanpos = [4,6]
        nonglycanpos = []
        nonglycanptm = []
        self.assertEqual(sgptest2.peptide,'GYLNCTR')
        self.assertEqual(sgptest2.nonglycanptm,nonglycanptm)
        self.assertEqual(sgptest2.glycanptm,glycanptm)
        self.assertEqual(sgptest2.nonglycanpos,nonglycanpos)
        self.assertEqual(sgptest2.glycanpos,glycanpos)
    # def test3(self):        

if __name__ == "__main__":
    unittest.main()    
        
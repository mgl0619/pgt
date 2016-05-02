# -*- coding: utf-8 -*-
"""
Created on Sun May 01 12:57:59 2016

@author: gang
"""

class MSScan():
    def __init__(self):
        self.id = 0
        self.scanNum = 0
        self.retentionTime = 0.0
        self.lowMz = 0.0
        self.highMz = 0.0
        self.peaksCount = 0
        self.polarity = 0
        self.scanType = 0
        self.basePeakMz = 0.0
        self.basePeakIntensity = 0.0
        self.totIonCurrent = 0
        self.precursorMz = 0
        self.precursorScanNum = 0
        self.precursorCharge = 0 
        self.msLevel = 0        
        self.msIntensityList = []
        self.actMethod = []
        self.precision = []        

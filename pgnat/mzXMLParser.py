# -*- coding: utf-8 -*-
"""
Created on Sat Apr. 16 10:05:42 2016
Modified on Apr.24 2016 11:00 2016

1) Use lxml.etree iterparse to parse mzXML fast_iter method learned from 
http://www.ibm.com/developerworks/library/x-hiperfparse/
2) Added methods/class to retrieve Instrument Info and Data Processing Info to mzXMLParser

@author: gang
"""

# -------------------------------------------------------------------------
#     This program is free software; you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation; either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.

#     Complete text of GNU GPL can be found in the file LICENSE.TXT in the
#     main directory of the program.
# -------------------------------------------------------------------------

from lxml import etree
import zlib
import base64
import numpy as np
import sys
import MSScan

__dataProcessingInfo_fieldNames = []

def _decode_base64_data_array(source, dtype, is_compressed):
    """Read a base64-encoded binary array.

    Parameters
    ----------
    source : str
        A binary array encoded with base64.
    dtype : str
        The type of the array in numpy dtype notation.
    is_compressed : bool
        If True then the array will be decompressed with zlib.

    Returns
    -------
    out : numpy.array
    """

    decoded_source = base64.b64decode(source.encode('ascii'))
    if is_compressed:
        decoded_source = zlib.decompress(decoded_source)
    output = np.frombuffer(decoded_source, dtype=dtype)
    return output
    
def mzxmlParser_xml(infile):
    return etree.parse(infile)
    
class mzXMLParser():
     __instrumentMSInfo_fieldNames = ["msManufacturer","msModel","msIonisation",\
                               "msDetector","msMassAnalyzer","msResolution"]
     __instrumentSWInfo_fieldNames = ["softwareVersion","softwareName","softwareType"]
     
     __instrumentOPInfo_fieldNames = ["operatorFirstName","operatorLastName","operatorPhone"]

     __dataProcessingSWInfo_fieldNames = ["softwareVersion","softwareType","softwareName"]
     
     __dataProcessingDataInfo_fieldNames = ["centroided","deisotoped","intensityCutoff","spotIntegration"]     
     
     #infile: input file name including path
     _inFile = None
    
     # instrument info
     _instrumentInfo = {}
    
      # data processing info
     _dataProcessingInfo = {}
     
      # scan data
     _scanData = None
     
     @property    
     def inFile(self):
         """ InFile Property"""
         return self._inFile
        
     @inFile.setter
     def x(self,inFilePath):
         self._inFile = inFilePath
    
     @property    
     def instrumentInfo(self):
         """ Instrument Info Property"""
         return self._instrumentInfo
        
     @property    
     def dataProcessingInfo(self):
         """ Data Processing Info. Property"""
         return self._dataProcessingInfo
         
     def _fast_iter(self,context,func,isSingleItem=True,eventStrList=['start'],*args):
        try:
            listOfXMLInfo = []
            for event, elem in context:
                for eventstr in eventStrList:
                    if event == eventstr:
                        xmlInfo = func(elem, *args)
                        #print xmlInfoNew
                        if isSingleItem:
                            if xmlInfo:
                                del context
                                return xmlInfo
                            else:
                                elem.clear()
                                while elem.getprevious() is not None:
                                    del elem.getparent()[0]
                        else:
                            listOfXMLInfo.append(xmlInfo)
                            print listOfXMLInfo
            del context
            return listOfXMLInfo
        except:
             print 'syntax error'
             e = sys.exc_info()[0]
             print e
             return {}    

     def _scanInfo(elem, tagName):
        xmlInfo = {};
        if((etree.QName(elem.tag)).localname==tagName):
           # retrieve properties of the tag
           eleattrib = dict(elem.attrib)
           if eleattrib:
               for key,value in eleattrib.items():
                   xmlInfo[key] = value
           xmlInfo['text'] = elem.text
           for child in elem:
               #print 'attrib=', child.attrib
               childName = (etree.QName(child.tag)).localname
               xmlInfo[childName] = dict(child.attrib)
               xmlInfo[(childName+'text')] = child.text
        return xmlInfo
    
     def _setInstrumentMSInfo(self,msKey,instrumentInfo):
         if instrumentInfo.get(msKey,None):
             self._instrumentInfo[msKey]=(instrumentInfo.get(msKey)).get('value',None)
         else:
             self._instrumentInfo[msKey]=None
             
     def _setInstrumentSWInfo(self,msKey,instrumentInfo):
         if instrumentInfo.get(u'software',None):
             self._instrumentInfo[msKey]=(instrumentInfo.get('software')).get(msKey[8:].lower(),None)
         else:
             self._instrumentInfo[msKey]= None
             
     def _setInstrumentOPInfo(self,msKey,instrumentInfo):
         if instrumentInfo.get(u'operator',None):
             self._instrumentInfo[msKey]=(instrumentInfo.get('operator')).get(msKey[8:],None)
         else:
             self._instrumentInfo[msKey]= None        
             
     def _setDataProcessingMSInfo(self,msKey,dataProcessingInfo):
         if dataProcessingInfo.get('software',None):
             print msKey[8:]
             self._dataProcessingInfo[msKey]=(dataProcessingInfo.get('software')).get(msKey[8:].lower(),None)
         else:
             self._dataProcessingInfo[msKey]= None
             
     def _setDataProcessingDataInfo(self,msKey,dataProcessingInfo):
         self._dataProcessingInfo[msKey]=dataProcessingInfo.get(msKey,None)
             
     def _scanInstrumentInfo(self):
         instrumentInfo = self._fast_iter(etree.iterparse(self.inFile,events = \
                                 ("start","end")),self._scanInfo,True,['start'],"msInstrument")
        
         for msKey in self.__instrumentMSInfo_fieldNames:
             self._setInstrumentMSInfo(msKey,instrumentInfo)
        
         for msKey in self.__instrumentSWInfo_fieldNames:
             self._setInstrumentSWInfo(msKey,instrumentInfo)
             
         for msKey in self.__instrumentOPInfo_fieldNames:
             self._setInstrumentSWInfo(msKey,instrumentInfo)
             
     
     def _scanDataProcessingInfo(self):
         dataProcessingInfo = self._fast_iter(etree.iterparse(self.inFile,events = \
                                 ("start","end")),self._scanInfo,True,['start'],"dataProcessing")
         print dataProcessingInfo                        
                                 
         for msKey in self.__dataProcessingDataInfo_fieldNames:
             self._setDataProcessingDataInfo(msKey,dataProcessingInfo)
                        
         for msKey in self.__dataProcessingSWInfo_fieldNames:
             self._setDataProcessingMSInfo(msKey,dataProcessingInfo)
    
     def __init__(self,path=None):
         if path:
             self._inFile = path
             self._scanInstrumentInfo()
             self._scanDataProcessingInfo()
             print self.instrumentInfo
             print self.dataProcessingInfo
    
     def readScanData(self,scanNum):
         context = etree.iterparse(self.inFile,events = ("start","end"))
         for event, elem in context:
             if event == 'start':
                  if((etree.QName(elem.tag)).localname=='scan'):
                     # retrieve properties of the tag
                     eleattrib = dict(elem.attrib)
                     if eleattrib['num']==str(scanNum):
                         msScan = MSScan
                         msScan.scanNum = scanNum
                         msScan.retentionTime = eleattrib['retentionTime']
                         msScan.lowMz = eleattrib['lowMz']
                         msScan.highMz = eleattrib['highMz']
                         msScan.polarity = eleattrib['polarity']
                         msScan.basePeakMz = eleattrib['basePeakMz']
                         msScan.basePeakIntensity = eleattrib['basePeakIntensity']
                         msScan.totIonCurrent = eleattrib['totIonCurrent']
                         msScan.msLevel = eleattrib['msLevel']
                            
                         for child in elem:
                            childName = (etree.QName(child.tag)).localname
                            if childName == 'peaks' :
                                eleattrib = dict(child.attrib)
                                msScan.precision = eleattrib['precision']
                                msScan.msIntensity = child.text
                         del context
                         return msScan
                     else:
                         elem.clear()
                         while elem.getprevious() is not None:
                             del elem.getparent()[0]    
                  else:
                      elem.clear()
                      while elem.getprevious() is not None:
                          del elem.getparent()[0]         
         return None
         
def main():
    infile   = 'test.mzXML'
    testmzXMLParser = mzXMLParser(infile)
    scan2 = testmzXMLParser.readScanData(6586)
    print scan2.msIntensity
    print 'end'
       
#    def retrieveScanData(self,scannum):
if  __name__ == '__main__':
    main()
          
    
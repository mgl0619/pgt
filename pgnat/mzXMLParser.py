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
    
def fast_iter(context,func,*args):
    xmlInfo = {}
    try:
        for event, elem in context:            
            xmlInfoNew = func(elem, *args)
            xmlInfo.update(xmlInfoNew)
            elem.clear()
            while elem.getprevious() is not None:
                del elem.getparent()[0]
        del context
    except:
         print 'syntax error'
         e = sys.exc_info()[0]
         print e
         return None
    return xmlInfo

def scanInfo(elem, tagName):
    xmlInfo = {};
    if((etree.QName(elem.tag)).localname==tagName):
       for child in elem:
           #print 'attrib=', child.attrib
           xmlInfo[(etree.QName(child.tag)).localname] = dict(child.attrib)
    return xmlInfo    
                
class mzXMLParser():
     __instrumentInfo_fieldNames = ["msManufacturer","msModel","msIonisation",\
                               "msDetector","msMassAnalyzer","softwareVersion",\
                               "softwareName","softwareType"]

     __dataProcessing_fieldNames = ["softwareVersion","softwareType","softwareName"]     
     
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
             
     def _setDataProcessingInfo(self,msKey,dataProcessingInfo):
         if dataProcessingInfo.get('software',None):
             self._dataProcessingInfo[msKey]=(dataProcessingInfo.get('software')).get(msKey[8:].lower(),None)
         else:
             self._dataProcessingInfo[msKey]= None        
    
        
     def _scanInstrumentInfo(self):
         instrumentInfo = fast_iter(etree.iterparse(infile,events = \
                                 ("start","end")),scanInfo,"msInstrument")
         for msKey in self.__instrumentInfo_fieldNames[0:4]:
             self._setInstrumentMSInfo(msKey,instrumentInfo)
        
         for msKey in self.__instrumentInfo_fieldNames[5:8]:
             self._setInstrumentSWInfo(msKey,instrumentInfo)
     
     def _scanDataProcessingInfo(self):
         dataProcessingInfo = fast_iter(etree.iterparse(infile,events = \
                                 ("start","end")),scanInfo,"dataProcessing")
         for msKey in self.__dataProcessing_fieldNames:
             self._setDataProcessingInfo(msKey,dataProcessingInfo)
    
     def __init__(self,path=None):
         if path:
             self._inFile = path
             self._scanInstrumentInfo()
             self._scanDataProcessingInfo()
             print self.instrumentInfo
             print self.dataProcessingInfo
    
    
#    def retrieveScanData(self,scannum):
if  __name__ == '__main__':
    emptyparser = mzXMLParser()
    infile   = 'test.mzXML'
    testmzXMLParser = mzXMLParser(infile)
          
    
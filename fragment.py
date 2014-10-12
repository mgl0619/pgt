# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 20:00:59 2014

@author: gangliu

"""
import re
import collections

def splitGlyPep(sgp):
    """Split glycopeptide into peptide, glycans and non-glycan modification
    return string
    
    @Syntax: 
        [peptide,glycan,nonglycanptm]=splitGlyPep(sgp)
    
    @Params:
        sgp: glycopeptide string, string type
        
    @Retuns:
        peptide: 
        glycan:
        nonglycanptm:
    
    @See also
    """
    
    # set up regular expression pattern for peptide, non-glycan ptm, glycan
    pep_regex = '[ARNDCQEGHILKMFPSTWYV]';
    nonglycan_regex = '(?<=<)[\+\-\da-z]+(?=>)';
    glycan_regex = '{[\+\-\d.a-z{}]+}(?=[A-Z])|(?<=[A-Z]){[\+\-\d.a-z{}]+}';
    pep_pattern = re.compile(pep_regex,re.UNICODE);
    nonglycan_pattern = re.compile(nonglycan_regex,re.UNICODE)
    glycan_pattern = re.compile(glycan_regex,re.UNICODE);
    
    # find peptide sequence and their position in SGP
    peptide = ''.join(pep_pattern.findall(sgp));
    print peptide
    #    aapos_in_pep = [];
    #    for match in pep_pattern.finditer(sgp):
    #      startpop = match.start();
    #      aaposinsgp.append(startpop);
    
    
    # find non-glycan ptm enclosed between angle brackets 
    nonglycanptm = nonglycan_pattern.findall(sgp)
    # print nonglycanptm
    
    nonglycanpos_in_sgp = [];
    for match in nonglycan_pattern.finditer(sgp):
       startpos            = match.start();
#       print startpos
       nonglycanptmpos     = len(pep_pattern.findall(sgp[0:startpos-1]));
       nonglycanpos_in_sgp.append(nonglycanptmpos);
#       print nonglycanpos_in_sgp
       
    
    # find glycan structure enclosed between curly brackets
    if(len(peptide)==0):
        glycanptm = sgp
        glycanpos_in_sgp = []
    else:
        glycanptm = glycan_pattern.findall(sgp);
        print glycanptm
        glycanpos_in_sgp = []
        for match in glycan_pattern.finditer(sgp):
           startpop         = match.start()
           glycanptmpos     = len(pep_pattern.findall(sgp[0:startpop]))
           glycanpos_in_sgp.append(glycanptmpos)      
    
    
    SGPelements = collections.namedtuple('SGPelements',['peptide',
       'nonglycanptm','glycanptm','nonglycanpos','glycanpos'])
    return SGPelements(peptide,nonglycanptm,glycanptm,nonglycanpos_in_sgp,
    glycanpos_in_sgp)
    
testsgp1 =  'GYM<o>KNCT<s>'
sgptest1 = splitGlyPep(testsgp1) 

testsgp2 = 'GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT{n{h{s}}{n{h{s}{f}}}}R'
sgptest2 = splitGlyPep(testsgp2) 

testsgp3 = '{n{h{s}}{n{h{s}{f}}}}'
sgptest3 = splitGlyPep(testsgp3) 

testsgp4 = 'GYLN{n{n{h{h{h{162.1}}}{h{h{h}}{h{h}}}}}}CT<+96>{n{h{s}}{n{h{s}{f}}}}R'
sgptest4 = splitGlyPep(testsgp4) 

testsgp5 = 'GYLN{n{n{h{h{h{h}}}{h{h{h}}{h{h}}}}}}CT<s>{n{h{s}}{n{h{s}{f}}}}R'
sgptest5 = splitGlyPep(testsgp5)     
    
    


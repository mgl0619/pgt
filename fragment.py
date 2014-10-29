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
        sgpelements = splitGlyPep(sgp)
    
    @Params:
        sgp: glycopeptide string, string type
        
    @Retuns:
        sgpelements.peptide: peptide element
        sgpelements.glycan: glycan element
        sgpelements.nonglycanptm: non-glycan modification element
    
    @See also joinGlyMat
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
    #    print peptide
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
        # print glycanptm
        glycanpos_in_sgp = []
        for match in glycan_pattern.finditer(sgp):
           startpop         = match.start()
           glycanptmpos     = len(pep_pattern.findall(sgp[0:startpop]))
           glycanpos_in_sgp.append(glycanptmpos)      
    
    
    SGPelements = collections.namedtuple('SGPelements',['peptide',
       'nonglycanptm','glycanptm','nonglycanpos','glycanpos'])
    return SGPelements(peptide,nonglycanptm,glycanptm,nonglycanpos_in_sgp,
    glycanpos_in_sgp)
    
def joinGlyPep(pep,gly,glyPos,nonglyMod,nonglyModPos):
    """Reconstruct the glycopeptide
    
    @Syntax: 
        sgp = joinGlyPep(peptide,glycan,nonglycanptm,glycanpos,modpos)
    
    @Params:
        peptide: peptide backbone
        glycan: glycan modification
        nonglycanptm: non-glycan modification
        
    @Retuns:
        sgp: glycopeptide structure in SGP format
    
    @See also splitGlyPep
    """ 
        
    if(len(pep)==0): # glycan only
        sgp =  gly;
    else:
        sgp = '';
        
        for i in range(len(pep)):
           print i  
           sgp = sgp + pep[i]
           
           ismodadded = 0;
           for k in range(len(nonglyModPos)) :
               singlemodpos = nonglyModPos[k]
               print singlemodpos
               if(i==singlemodpos):
                   sgp = sgp + '<'+nonglyMod[k]+'>'                 
                   print sgp
                   ismodadded = 1
                   break
           if(ismodadded==1):
               # nonglyMod.remove(nonglyMod[k])
               del nonglyMod[k]
               del nonglyModPos[k]
               #nonglyModPos.remove(singlemodpos)            
           
           isglyadded = 0
           for j in range(len(glyPos)):
               singleglypos = glyPos[j]
              # print singleglypos
               if(i==singleglypos):
                   sgp = sgp + gly[j]
                   isglyadded = 1                   
                   print sgp
                   break
               
           if(isglyadded==1):
               del gly[j]
               del glyPos[j]
               #gly.remove(gly[j])           
               #glyPos.remove(singleglypos)               
    return(sgp)
    
def swapAA(pepseq,swapoption='random'):
    """swap amino acid  
    
    @Syntax: 
        newpepseq = swapAA(pepseq,option)
    
    @Params:
        pepseq: peptide sequence
        option: swap option.
        
    @Retuns:
        newpepseq: new peptide sequence
    
    @Examples:
     >     
    
    @See also splitGlyPep
    """ 
    
    if(swapoption=='random'):
        import random
        pepseq = ''.join(random.sample(pepseq,len(pepseq)));
    elif(swapoption=='reverse'):
        pepseq = pepseq[::-1]
    else:
        print "unknown optioin"
        raise
    return(pepseq)

def decoysgp(sgp,decoypepopt='random',decoyglyopt=''):
    """swap amino acid  
    
    @Syntax: 
        newpepseq = swapAA(pepseq,option)
    
    @Params:
        pepseq: peptide sequence
        option: swap option.
        
    @Retuns:
        newpepseq: new peptide sequence
    
    @Examples:
     >     
    
    @See also splitGlyPep
    """ 
    
    glycan_regex = '(?<={)[a-z]';
    glycan_pattern = re.compile(glycan_regex,re.UNICODE);    
    
    sgpelements = splitGlyPep(sgp,decoypepopt)
    sgppep      = swapAA(sgpelements.peptide)
    
    maxmwchange = 40
    numglycans  = len(sgptest1.glycanptm)
    for i in range(numglycans):
        glystring   = sgptest1.glycanptm[i]
        glyresidues = glycan_pattern.findall(glystring)
        for j in glyresidues:
            glyresiduesmw(j)= glycanMW(glyresidue[j])
            
#maxmwchange = 50;
#for i = 1: numglycans
#    glystring     = glyMat(1,i).struct; 
#    glyresidues   = regexp(glystring,'(?<={)[a-z]','match');
#    glyresiduesmw = zeros(glyMat(1,i).len,1);
#    for j = 1 : glyMat(1,i).len
#       glyresiduesmw(j)=Glycan.glycanMSMap(glyresidues{j});
#    end
#    mwchange_residues  = randfixedsum(glyMat(1,i).len,1,0,-1,1)*maxmwchange;
#    mwchange_residues  = mwchange_residues + glyresiduesmw;  
#    
#    replacementstring  = cellstr(num2str(mwchange_residues));
#    replacementstring  = strtrim(replacementstring);
#    glyMat(1,i).struct = regexprep(glystring,'(?<={)[a-z]',replacementstring,'once');
#end
    
            
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
    
pep          = 'GYMKNCT';
gly          = ['{n{h{s}}}']
glyPos       = [4]
nonglyMod    = ['o','s']
nonglyModPos = [2,6]    
sgptest      = joinGlyPep(pep,gly,glyPos,nonglyMod,nonglyModPos)
aaseq2       = swapAA(pep,'random')

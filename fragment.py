# -*- coding: utf-8 -*-
"""
Created on Thu Oct 09 20:00:59 2014

@author: gangliu

"""
import re
import collections
import massconstant

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
          # print i  
           sgp = sgp + pep[i]
           
           ismodadded = 0;
           for k in range(len(nonglyModPos)) :
               singlemodpos = nonglyModPos[k]
             #  print singlemodpos
               if(i==singlemodpos):
                   sgp = sgp + '<'+nonglyMod[k]+'>'                 
                   #print sgp
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
                  # print sgp
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
        print "unknown option"
        raise
    return(pepseq)

def decoysgp(sgp,decoypepopt='random'):
    """swap amino acid  
    
    @Syntax: 
        newsgp = decoysgp(sgp,decoypepopt,decoyglyopt)
    
    @Params:
        sgp: glycopeptide sequence
        decoypepopt: decoy option for peptide 
        decoyglyopt: decoy option for glycan
        
    @Retuns:
        newsgp: new peptide sequence
    
    @Examples:
    
      > SmallGlyPep='M<o>{n{h{s}}}GHKL<o>FLM{n{h{x}}}L';
        decoyopt='flip';
        decoysgp(SmallGlyPep, decoyopt)     
    
    @See also splitGlyPep    
    """ 
    import random
    import numpy as np
    
    glycan_pattern = re.compile('(?<={)[a-z]',re.UNICODE);
    sgpelements    = splitGlyPep(sgp)
    numaa = len(sgpelements.peptide)
    
    # amino acid decoy    
    if(decoypepopt=='random'):
        newindex       = random.sample(range(numaa),numaa)        
        newpepseq      = ''.join(sgpelements.peptide[i] for i in newindex);
    elif(decoypepopt=='reverse'):
        newpepseq         = sgpelements.peptide[::-1]
        newindex          = range(numaa-1,-1,-1)
    else:
        print "unknown option"
        raise
    
    # update new positions for glycan and non-glycan ptm
    newnonglycanpos_in_sgp =[]
    newglycanpos_in_sgp = []
    numglycan = len(sgpelements.glycanptm)
    numnonglycanptm = len(sgpelements.nonglycanptm)
    for i in range(numaa):
        for j in range(numnonglycanptm):
            if(newindex[i]==sgpelements.nonglycanpos[j]):
                newnonglycanpos_in_sgp.append(i)
        
        for j in range(numglycan):
            if(newindex[i]==sgpelements.glycanpos[j]):
                newglycanpos_in_sgp.append(i)

    # glycan decoy
    maxmwchange = 40
    numglycans  = len(sgpelements.glycanptm)
    newglycanptm =[]
    
    for i in range(numglycans):
        replacementstring = []
        glyresiduesmw = []
        glystring   = sgpelements.glycanptm[i]
        glyresidues = glycan_pattern.findall(glystring)
        nummono = len(glyresidues)
        for j in range(nummono):
            glyresiduesmw.append(massconstant.glycanresidue_mass[glyresidues[j]])
        glyresiduesmw_change  = (massconstant.staffordRandFixedSum(nummono,1,1)-1./10)*maxmwchange
        glyresiduesmwarray    = np.array(glyresiduesmw,ndmin=2)
        glyresiduesnewmw      = glyresiduesmw_change + glyresiduesmwarray

        for j in range(nummono):
            replacementstring.append(str(glyresiduesnewmw[0,j]))
            glystring=re.sub('(?<={)[a-z]',replacementstring[j],glystring,1)
            
        newglycanptm.append(glystring)
    newsgp=joinGlyPep(newpepseq,newglycanptm,newglycanpos_in_sgp,sgpelements.nonglycanptm,newnonglycanpos_in_sgp)             
    return(newsgp)   
  
            
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
newsgp      = decoysgp(testsgp2)

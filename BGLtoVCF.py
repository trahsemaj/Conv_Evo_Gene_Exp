#! /usr/bin/env python 

import sys,getopt,module
'''

BGLtoVCF.py [options] unphased.vcf phased.bgl > phased.vcf
input: an unphased vcf file, and the phased bgl version of the file
output: a phased vcf file, using information from the bgl file

[options]:
    -i [LIST] individuals in the bgl file to write to the vcf, usually 2,5,8
    -v [LIST] indexes of individuals in the original vcf file 0,1,2,5,6

'''

##BGLtoVCF.py -i 2 -v 2 PAXB_VTP_RNA.vcf JCH003A_trio_tempxp.bgl > test_bgl.vcf

optlist,args = getopt.getopt(sys.argv[1:],'i:v:')

VCFFILE = args[0]
BGLFILE = args[1]
INDIV = []
INDICES = []

for a in optlist:
    if a[0] == '-i':
        INDIV  = map(lambda x: int(x), a[1].split(','))
    if a[0] == '-v':
        INDICES  = map(lambda x: int(x), a[1].split(','))



##takes a filepath as input, as well as individuals to grab
##returns a dict in the form BGLDICT[pos] = ['A,A','A,T','A,T']
def loadBgl(bfile,inds):
    
    bgldict = {}
    
    f = open(bfile,'r')
    
    for line in f:
        if not line.startswith('M'):
            continue
        
        line = line.strip()
        
        split = line.split()
        
        pos = split[1].split(',')[1]
        
        rawgenos = split[2:]
        selectgenos = []
        ##for each individual, add the phased genotypes to the list
        for i in inds:
            temp = []
            temp.append(rawgenos[i * 2])
            temp.append(rawgenos[(i * 2) + 1]) 
            selectgenos.append(temp)
        
        bgldict[pos] = selectgenos
        
    return bgldict


def phaseVCF(invcf, bgldict,indices):
    
    f = open(invcf, 'r')
    
    
    for line in f:
        line = line.strip()

        if line.startswith('##'):
            print line
            continue
        if line.startswith('#CHROM'):
            split = line.split('\t')
            print '\t'.join(split[:9])
            continue
        split = line.split('\t')
        
        
        chr = split[0]
        pos = split[1]
        ref = split[3]
        alt = split[4]
        stem = split[:9]
        
        ##filter out the indices not used in the phased bcf file (which usually only has offspring)
        rawunphased = split[9:]
        unphased = []
        for i in indices:
            unphased.append(rawunphased[i])
        key = chr + ':' + pos
        
        ##marker dict, returns 0 or 1 depending on the state of the marker
        markers = {}
        markers[ref] = '0'
        markers[alt] = '1'
        
        if not key in bgldict:
            continue
        
        ##for phasing DP - read counts
        AD = []
        fields = []
        ##print >> sys.stderr,unphased
        for gt in unphased:
            ##print >> sys.stderr,gt
            ##if len(gt.split(':')) < 2:
            if gt.startswith('.'):
                AD.append(['.','.'])
                fields.append([''])
                continue

            read_counts = gt.split(':')[1].split(',')
            AD.append(read_counts)
            ##fields which aren't GT or AD
            fields.append(gt.split(':')[2:])
        ##print >> sys.stderr,"AD " +str(AD)
        phased = bgldict[key]
        ##print  >> sys.stderr,phased,fields
        ##for each individual write the genotype
        for i,ind in enumerate(phased):
            
            rawReads = AD[i]
            ##print >> sys.stderr,"rawReads " +str(rawReads)
            curFields = fields[i][:]
            if ind[0] in markers and ind[1] in markers:
                genotype = markers[ind[0]] + '|' + markers[ind[1]]
                supReads = rawReads[0] + ',' + rawReads[1]
            ##I don't know when this would happen, but setting it to ref seems safest
            else:
                genotype = '0|0'
                supReads = ','.join(rawReads)
            ###insert the sup reads into fields at the 0th index - this could be done better but didn't want to re-write everything
           ## print >> sys.stderr,fields
            curFields.insert(0,supReads)
            ##print >> sys.stderr,fields
            curFields.insert(0,genotype)
            ##print >> sys.stderr,fields
            toappend = ':'.join(curFields)
            stem.append(toappend)
        ##print >> sys.stderr,stem
        toprint = '\t'.join(stem)
        print toprint    
        


BGLDICT = loadBgl(BGLFILE, INDIV)
phaseVCF(VCFFILE,BGLDICT,INDICES)
    
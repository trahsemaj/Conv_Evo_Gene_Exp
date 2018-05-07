#! /usr/bin/env python

import sys


'''
VCFtoBeagle.py in.vcf [TRIOS]
[TRIOS] = STR, in the from 0,2,3/4,2,3/5,2,3 with each trio seperated by / and each collum in the vcf seperated by ,
'''


VCF = sys.argv[1]
RAWTRIOS = sys.argv[2]
TRIOS = []
temp = RAWTRIOS.split('/')
for t in temp:
    TRIOS.append(map(lambda x: int(x),t.split(',')))
    

header = ['I','id']
num_ind = 0
for t in TRIOS:
    for i in t:
        header.append(str(num_ind))
        header.append(str(num_ind))
        num_ind +=1
print '\t'.join(header)


f = open(VCF,'r')
prev = 0
chr = ''
count = 0
for line in f:
    
    if line.startswith('#'):
        continue
    
    line = line.strip()
    split = line.split()
    
    ##ref/alt allele, reference with markers[genotype]
    markers = split[3:5]
    rawGenos = split[9:]
    chr = split[0]
    pos = int(split[1])
    '''if pos <= prev:
        continue'''
    prev = pos
    toprint = ['M',str(count) + ',' + chr + ':' + str(pos)]
    count +=1
    for t in TRIOS:
        for ind in t:
            ##append '?' for missing data
            if rawGenos[ind].startswith('.'):
                toprint.append('?')
                toprint.append('?')
                continue
            ##extracts the genotypes and turns them into ints
            genos = map(lambda x: int(x) , rawGenos[ind].split(':')[0].split('/'))
            ##ignore multiallelic
            if genos[0] > 1 or genos[1] > 1:
                continue
            toprint.append(markers[genos[0]])
            toprint.append(markers[genos[1]])
            
        
    
    ##if toprint is full length, no missing or multiallelic were encountered, so continue
    if len(toprint) == (num_ind * 2) + 2:
        print '\t'.join(toprint)
    
    
    
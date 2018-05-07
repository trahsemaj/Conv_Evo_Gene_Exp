#! /usr/bin/env python
## check_mendelian_trio.py
##AUTHOR: Sharon R. Browning, sguy@uw.edu 
##print "usage: python check_mendelian_trio.py trios.bgl new_trios.bgl missing_code"

import sys
triosfile = open(sys.argv[1])
outfile = open(sys.argv[2],"w")
missing_code = sys.argv[3]

def isConsistent(a1,a2,c,missing_code):
    return c == missing_code or a1 == missing_code or a2 == missing_code or c == a1 or c == a2
count = 0
total = 0
for line in triosfile.xreadlines():
  elements = line.split()
  if elements[0] == "M":
    alleles = elements[2:]
    rs = elements[1]
    for j in range(0,len(alleles),6):
        b1 = isConsistent(alleles[j],alleles[j+1],alleles[j+4],missing_code) and isConsistent(alleles[j+2],alleles[j+3],alleles[j+5],missing_code)
        b2 = isConsistent(alleles[j],alleles[j+1],alleles[j+5],missing_code) and isConsistent(alleles[j+2],alleles[j+3],alleles[j+4],missing_code)
        total +=1
        if not (b1 or b2):
          ##print "mendelian inconsistency at",rs,"trio",j/6+1,":",alleles[j:(j+6)],"(setting all to missing)"
          count +=1
          for i in range(6):
            alleles[j+i] = missing_code
    print >> outfile, "M",rs,
    for a in alleles:
        print >> outfile,a,
    print >> outfile
  else:
      print >> outfile, line,
    
print total,count
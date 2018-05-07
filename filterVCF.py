#!/usr/bin/python
import sys,getopt

'''
VCF_filter.py [options] IN.vcf
Generic VCF filter, writes to STDOUT.  Easier to write this than use VCFTools (POS)
    -d [INT] The min depth needed to keep a site [0]
    -g [INT] The min genotype quality needed to output a given variant
    -b [FLOAT] The minimum inbalance at HET variants, given as a ratio. E.G. -b 3 will filter out 3:1 and 1:3 REF:ALT
    -i [INT] The index of the individual whose fields are being considered [0]
'''
##minimum depth needed to output line
MIN_DEPTH = 0
##minimum genotype quality
MIN_GQ = 0
##minimum inbalance of het SNPs
MIN_INBAL = 0
INDEX = 0
args = sys.argv[1:]    
optlist, args = getopt.getopt(args,'d:g:b:i:')
    
    
##Add filters as needed    
for a in optlist:
    if a[0] == '-d':
        MIN_DEPTH = int(a[1])
    if a[0] == '-g':
        MIN_GQ = int(a[1])
    if a[0] == '-b':
        MIN_INBAL = float(a[1])+1e-20
    if a[0] == '-i':
        INDEX = int(a[1])
VCF = args[0]

##input: a fpath pointing to a vcf file, and the index of an individual, with filtering options
##output: the filtered vcf printed to stdout
##filtering rules applied in order (missing,depth,quality,inbalance)
def filterVCF(vcf,index,min_depth=0,min_gq=0,min_inbal=0):
    
    f = open(vcf,'r')
    total_variants = 0
    missing = 0
    gq_filtered = 0
    depth_filtered = 0
    inbalanced_filtered = 0
    passed_filters = 0
    for line in f:
        line = line.strip()
        
        if line.startswith('#'):
            print line
            continue
        
        split = line.split('\t')
        
        ##skip malformed lines
        if len(split) < 8:
            continue
        total_variants += 1
        ##all the genotype calls
        rawGenotypes = split[9:]
        ##the genotype of the ind in question
        indGeno = rawGenotypes[index]
        
        
        ##gt line GT:AD:DP:GQ:PL
        gt_info = indGeno.split(':')
        genotype = gt_info[0]
        if genotype.startswith('.'):
            missing+=1
            continue
        
        ##check if missing:
        AD = gt_info[1]
        if AD[0] =='.':
            missing+=1
            continue
        
        ##depth_filter
        cur_depth = int(gt_info[2])
        if cur_depth < min_depth:
            depth_filtered +=1
            continue
        
        
        ##gq filter
        cur_gq = int(gt_info[3])
        if cur_gq < min_gq:
            gq_filtered += 1
            continue
        
        ##balance - only look at het
        
        if genotype == '0/1':
            
            read_counts = map(lambda x: float(x)+1e-20, gt_info[1].split(','))
            
            ##if the read count ratio is inbalanced
            if (read_counts[0]/read_counts[1] > min_inbal) or (read_counts[0]/read_counts[1] < (1./min_inbal)):
                inbalanced_filtered +=1
                continue
      
    
        ##if all filters passed, print the line
        print line
        passed_filters += 1
        
        
    ##print stderr report
    print >> sys.stderr,'%d total variants checked' % total_variants
    print >> sys.stderr,'%d (%f) missing genotype calls' % (missing,float(missing)/total_variants)
    print >> sys.stderr,'%d (%f) passed all filters' % (passed_filters, float(passed_filters)/total_variants)
    print >> sys.stderr,'%d (%f) failed read_depth < %d' % (depth_filtered, float(depth_filtered)/total_variants, min_depth)
    print >> sys.stderr,'%d (%f) failed gq < %d' %(gq_filtered, float(gq_filtered)/total_variants,min_gq)
    print >> sys.stderr,'%d (%f) failed inbalance > %f' % (inbalanced_filtered, float(inbalanced_filtered)/total_variants,min_inbal)
    
    
filterVCF(VCF,INDEX,min_depth=MIN_DEPTH,min_gq=MIN_GQ,min_inbal=MIN_INBAL)


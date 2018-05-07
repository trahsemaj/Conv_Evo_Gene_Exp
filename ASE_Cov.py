'''
series of functions to look at covriance of ASE data
'''

import numpy as np
import module
import math
import matplotlib.pyplot as plt
import random
import scipy.spatial.distance
import scipy.stats
from operator import itemgetter

##preforms the MantelHaenszel repeated test of independance
##array is in the form 2x2xR, where R is the number of replicates. [ [[A,B],[C,D]],....]
##uses the continuity correction
def MantelHaenszel(array):

    numsum = 0
    denomsum = 0

    for table in array:

        a = float(table[0][0]) + 1e-20
        b = float(table[1][0]) + 1e-20
        c = float(table[0][1]) + 1e-20
        d = float(table[1][1]) + 1e-20

        
        n = a + b + c + d
        
        
        numsum += (a - ((a+b)*(a+c)/n ))

        denomsum += (a+b)*(a+c)*(b+d)*(c+d)/(n**3-n**2)
    ##with continuity
    chi2 = ((abs(numsum) - .5) ** 2 ) / denomsum
    ##chi2 = ((abs(numsum)) ** 2 ) / denomsum
    p =  stats.chisqprob(chi2,1)
    return chi2, p


##loads all replicates from a given list of SAM files
def loadSAMs(samlist,CUTOFF = 75):
    
    ##stores read counts in the form counts[ENSID] = [[F1,M1],[F2,M2],....]
    counts = {}
    
    for sam in samlist:

        samcount,samlength = countSAM(sam)
        
        samids = samcount.keys()
        print len(samcount.keys())
        for id in samids:
            
            
            if id not in counts:
                
                counts[id] = [samcount[id]]
            
            else:
                counts[id].append(samcount[id])
                
                
    ids = counts.keys()
    ##only keep the id if in all 5 replicates
    for id in ids:
        
        if len(counts[id]) < len(samlist):
            del counts[id]
            continue
        
        listsum = sum(module.flatten(counts[id]))
        
        if listsum < CUTOFF:
            del counts[id]
            continue
        
        ##print counts[id]
        
    return counts

##loads all replicates from a SAM list, needing at least CUTOFF reads in at least MINREPS replicates
def loose_loadSAMs(samlist,CUTOFF = 10,MINREPS=3):
    
    ##stores read counts in the form counts[ENSID] = [[F1,M1],[F2,M2],....]
    counts = {}
    
    for i,sam in enumerate(samlist):

        samcount,samlength = countSAM(sam)
        
        samids = samcount.keys()
        print len(samcount.keys())
        for ensid in samids:
            
            
            if ensid not in counts:
                ##initialize as no expression in all reps
                counts[ensid] = [[0,0] for j in range(len(samlist))]
            
            
	    counts[ensid][i] = samcount[ensid]
                
                
    ids = counts.keys()
    ##only keep the id if in all 5 replicates
    for ensid in ids:
        
	cur_counts = counts[ensid]
	##print cur_counts
	over_threshold = 0
	##check how many reps have more than CUTOFF reads
	for rep in cur_counts:
	    reads = sum(rep)
	    ##print rep,reads
	    if reads > CUTOFF:
		over_threshold += 1
	##print over_threshold
	##if less than MINREPS had more than CUTOFF reads, filter the gene
	if over_threshold < MINREPS:
	    del counts[ensid]

        
    return counts

##input: a SAM file with only uniquely mapping reads
##output: the counts of uniquely mapping reads for each gene
##loads all replicates from a given list of SAM files
##out.counts format 
##ENSID	F_length	M_length	F_Counts	M_counts
def countSAM(samfile):
    ##in the form dict[ensid]=[F_length,M_length]
    lengthDict = {}

    ##in the form dict[esnid]=[F_counts,M_counts]
    countDict = {}





    f = open(samfile,'r')

    for line in f:
        
        line = line.strip()
        
        ##if a header line
        if line.startswith('@SQ'):
            
            split = line.split('\t')
            
            id = split[1].split(':')[1].split('_')[0]
            allele = int(split[1].split(':')[1].split('_')[1]) - 1
            length = split[2].split(':')[1]
            
            
            ##if id isn't in the dictionary, initialize
            if id not in lengthDict:
                lengthDict[id] = ['-1','-1']
            
            lengthDict[id][allele] = length
                
        
        ##if not the header, its a read
        if not line.startswith('@'):
            
            split = line.split('\t')
            
            id = split[2].split('_')[0]
            allele = int(split[2].split('_')[1]) - 1
            
            ##if id isn't in the dictionary, initialize
            if id not in countDict:
                countDict[id] = [0,0]
                
            countDict[id][allele] += 1
        ##else:
            ##print line
    
    
    ids = countDict.keys()
    
    ids.sort()
    
    for id in ids:
        
        F_length = lengthDict[id][0]
        M_length = lengthDict[id][1]
        F_count = countDict[id][0]
        M_count = countDict[id][1]
        
        ##print '\t'.join(map(lambda x: str(x),[id,F_length,M_length,F_count,M_count]))
        
    return countDict, lengthDict

##input: a gtf file, (in ensGene format)
##output: a dict in the from startDict[name]=chrpos_start,chrpos_end(of exon1)
def loadGTFStart(ingtf):
    
    startDict = {}
    
    f = open(ingtf,'r')
    
    for line in f:
        line = line.strip()
        split = line.split('\t')
        
        if len(split) < 8:
            continue
        
        chr = split[0]
        start = split[3]
        end = split[4]
        
        ensid = split[8].split('"')[1]
        
        if ensid not in startDict:
            chrpos_start = chr + ':' + start
            chrpos_end = chr + ':' + end
            startDict[ensid] = (chrpos_start,chrpos_end)
            
    return startDict

##takes a dictionary of counts across ASE replicates
##returns a list of counts sorted by gene order
##if CHR != '', then only return genes on the given chromosome
def sort_CountDict(sam_count_dict,GTF_file='/data/James/BED/ensGene.gtf',CHR=''):

    
    sorted_ase_list = []

    gtfstarts = loadGTFStart(GTF_file)

    ensids = sam_count_dict.keys()
    
    ##keeps track of the amount of bases in all prev chromosomes - chrI has no offset, chrII offset is len(chrI), ect
    chr_offset = module.chrpostoDict()
    
    gnames = module.ENSidDict()
    
    print len(ensids)




    if CHR:
	filtered = []
	for ensid in ensids:
	    ##print gtfstarts[ensid]
	    start_chrpos = gtfstarts[ensid][0]
	    if start_chrpos.split(':')[0] == CHR:
		filtered.append(ensid)

	ensids = filtered

    print len(ensids)

    sorted_ids = sorted(ensids,key=lambda id: int(gtfstarts[id][0].split(':')[1]) + chr_offset[gtfstarts[id][0].split(':')[0]])

    for gene in sorted_ids:
	sorted_ase_list.append(sam_count_dict[gene])

    return sorted_ase_list

    
##Input: a dictionary of fresh/marine unique read counts
##Output: a dictionary of ASE values, log2(F_count/M_counts)
def counts_toASE(sam_count_dict,CUTOFF=10):

    ensIDs = sam_count_dict.keys()

    ase_dict = {}

    for gene in ensIDs:
	

	counts = sam_count_dict[gene]

	ase_vals = []
	total = 0
	for f,m in counts:
	    ase = math.log(float(f+1e-20)/(m+1e-20),2)
	    if f+m < CUTOFF:
		     continue
	    ase_vals.append(ase)
	    total += f+m



	ase_dict[gene]=ase_vals
    
    return ase_dict


def plotCorrMatrix(matrix,CMAP_TYPE='Reds'):

    cmap = plt.get_cmap(CMAP_TYPE)

    color_matrix = []

    print np.array(matrix)

    matrix = np.array(matrix)


    for row in list(matrix):
	newrow = []
	##if np.array(row).mean() > 1:
	    ##continue
	for i in row:		
	    ##newrow.append(cmap((float(i)+1)/2))
	    newrow.append(cmap((2-i)))
	color_matrix.append(newrow)


    plt.imshow(color_matrix, aspect='auto',interpolation='nearest')
    ##plt.show()

##takes two lists length m and n of ase values as input
##returns a mxn matrix of the euclidean distance between each element in the lists
def list_distance(list_a,list_b):
    
    distance_matrix = []


    for ase_vals_a in list_a:
	temprow = []

	for ase_vals_b in list_b:
	    ase_vals_a = np.array(ase_vals_a)
	    ase_vals_b = np.array(ase_vals_b)
	    dist = np.linalg.norm(ase_vals_a-ase_vals_b)
	    ##dist = math.fabs(ase_vals_a.mean()-ase_vals_b.mean())
	    ##dist = scipy.spatial.distance.mahalanobis(ase_vals_a,ase_vals_b,np.linalg.inv(np.cov(list_a).T))
	    temprow.append(dist)
	distance_matrix.append(temprow)

    return distance_matrix

##takes a 2D array (numpy array)
##returns a 1D array, with the rows collapsed to their mean values
def collapse_array(input_array):

    toret = []
    for i in input_array:
	i = np.array(i)
	##if math.fabs(i.mean()) > 3:
	    ##continue
	toret.append(i.mean())

    return np.array(toret)
	    
##takes a sorted collapsed ase list as input
##plots the correlation score with a given window size
def binned_mean(ase_list,window):
    
    N = len(ase_list)

    mean_blocks = []

    xs = np.arange(window,N-window+1)
    zeros = np.zeros(len(xs))

    for i in range(window,N-window+1):
	cur_mean = np.mean(ase_list[i-window:i+window])

	##print ase_list[i-window:i+window]
	mean_blocks.append(cur_mean)

    plt.plot(xs,mean_blocks)
    plt.plot(xs,zeros,ls='dashed',color='red')
    ##plt.show()

##takes a sorted collapsed ase list as input
##plots the correlation score with a given window size
def binned_correlation(ase_list,window,step):
    
    N = len(ase_list)

    mean_blocks = []

    print ase_list.mean()

    xs = np.arange(window,N-window+1-step,step)
    zeros = np.zeros(len(xs))

    for i in range(window,N-window+1-step,step):
	prev_vals = ase_list[i-window:i+window]
	cur_vals = ase_list[i-window+step:i+window+step]
	##print prev_vals,cur_vals
	cur_corr = np.correlate(prev_vals,cur_vals)

	##print ase_list[i-window:i+window]
	mean_blocks.append(cur_corr)

    plt.plot(xs,mean_blocks)
    plt.plot(xs,zeros,ls='dashed',color='red')
    ##plt.show()


##input: a sorted_ase array, and distance, 1+ the number of genes between two genes
##returns the mean expression divergence between genes that are distance apart
def expressionDiff_distance(sorted_ase, distance):

    expressiondiff = []

    for i in range(len(sorted_ase)):
	
	curr_exp = sorted_ase[i]
	
	curr_diff_mean = 0

	##only look to the right for the first few genes
	if i - distance < 0:
	    upper = sorted_ase[i+distance]
	    curr_diff_mean = math.fabs(curr_exp-upper)/(math.fabs(curr_exp) + math.fabs(upper) + 1e-20)
	elif i+distance >= len(sorted_ase):
	    lower = sorted_ase[i-distance]
	    curr_diff_mean = math.fabs(curr_exp-lower)/(math.fabs(curr_exp) + math.fabs(lower) + 1e-20)
	else:
	    upper = sorted_ase[i+distance]
	    lower = sorted_ase[i-distance]
	    
	    upper_diff =  math.fabs(curr_exp-upper)/(math.fabs(curr_exp) + math.fabs(upper) + 1e-20)
	    lower_diff =  math.fabs(curr_exp-lower)/(math.fabs(curr_exp) + math.fabs(lower) + 1e-20)
	    curr_diff_mean = (upper_diff)## + lower_diff)/2
	    ##min(upper_diff,lower_diff)
	expressiondiff.append(curr_diff_mean)

    mean_diff = sum(expressiondiff)/len(expressiondiff)
    return mean_diff

def writeAseBga(name,ase_dict, outfile):
    gtfstarts = loadGTFStart('/data/James/BED/ensGene_collapsed.gtf')
    w = open(outfile,'w')
    ##write the BGA header
    w.write('track type=bedGraph name=%s description=%s\n' % (name, name))

    for ensid in ase_dict:
        chrpos_start,chrpos_end = gtfstarts[ensid]
        chrome,start = chrpos_start.split(':')
        chrome,end = chrpos_end.split(':')
        ase = ase_dict[ensid]
        ##log2ase = math.log(ase/(1-ase),2)

        line = '\t'.join([chrome,start,end,str(ase)])
        w.write(line + '\n')

##input: a ASE.l file
##output: the dictionary form of the file
def loadASE_list(ase_file):

    toret = {}

    f = open(ase_file,'r')

    for line in f:
	line = line.strip()
	split = line.split('\t')
	ensid = split[0]
	ase = float(split[1])
	toret[ensid]=ase

    return toret

##takes a sorted ase list
##plots the up vs dwon matrix
def plot_updown(sorted_ase,cmap = 'seismic'):

    N = len(sorted_ase)
    
    updown_matrix = []
    
    colors = plt.get_cmap(cmap)
    
    for indI,i in enumerate(sorted_ase):
	temp = []
	for indJ,j in enumerate(sorted_ase):
	    ##print i,j
	    '''if (i > 0 and j > 0) or (i < 0 and j < 0):
		temp.append(1)
	    else:
		temp.append(0)
	    '''
	    if indJ >= indI:
		diff = (math.fabs(i-j)/((math.fabs(i) + math.fabs(j))+ 1e-20))
	    
		diff = (diff**.5)/2

		distance = 1-(diff)

		if i > 0:
		    distance = diff

		temp.append(colors(distance))
	    else:
		temp.append(colors(.5))
	updown_matrix.append(temp)
    ##print updown_matrix
    plt.imshow(updown_matrix, aspect='auto',interpolation='nearest',cmap='seismic')


def plot_ase_expDiff(sorted_ase,sorted_expDiff,cmap = 'seismic'):
    
    N = len(sorted_ase)
    
    updown_matrix = []
    
    colors = plt.get_cmap(cmap)
    
    for indI,i in enumerate(sorted_ase):
	temp = []
	for indJ,j in enumerate(sorted_ase):
	    ##print i,j
	    '''if (i > 0 and j > 0) or (i < 0 and j < 0):
		temp.append(1)
	    else:
		temp.append(0)
	    '''
	    if indJ >= indI:
		diff = (math.fabs(i-j)/((math.fabs(i) + math.fabs(j))+ 1e-20))
	    
		diff = (diff**.5)/2

		distance = 1-(diff)

		if i > 0:
		    distance = diff

		temp.append(colors(distance))
	    else:
		x = sorted_expDiff[indI]
		y = sorted_expDiff[indJ]
		diff = (math.fabs(x-y)/((math.fabs(x) + math.fabs(y))+ 1e-20))
	    
		diff = (diff**.5)/2

		distance = 1-(diff)
		if y > 0:
		    distance = diff
		temp.append(colors(distance))
	updown_matrix.append(temp)
    ##print updown_matrix
    plt.imshow(updown_matrix, aspect='auto',interpolation='nearest',cmap='seismic')



##takes a dict of read counts (pre-filtered)
##returns a dict of p values and significant p-values
def countStats(counts,fdr):
    
    ##dict of all pvalues
    pdict = {}
    
    ##dict of all significant p-values
    sigdict = {}
    
    ##dict of all esitmated ASE values
    asedict = {}
    
    
    ids = counts.keys()
    
    tests = len(ids)
    
    
    
    for id in ids:
        table = []
        
        readCounts = counts[id]
        
        F_counts = 0
        M_counts = 0
        
        ##for each individual make a 2x2 O vs E table
        for indiv in readCounts:
            
            F_counts += indiv[0]
            M_counts += indiv[1]
            
            expected = sum(indiv)/2
            
            table.append( [ indiv,[expected,expected]])
        
        asedict[id] = float(F_counts) / (F_counts + M_counts)


        chi2, p = MantelHaenszel(table)
        
        pdict[id] = p
        
        
    plist = sorted(pdict.items(), key=itemgetter(1))
    
    for i,pair in enumerate(plist):
        id,pval = pair
        
        ##print pval
        ##print (i+1)/float(tests) * fdr
        
        if pval > (i+1)/float(tests) * fdr:
            break
        
        ##print i,id,pval
        sigdict[id] = pval
    
    return pdict,sigdict,asedict



SAMBASE = '/data/James/AlleleSpecificExpression/VTP/CERC/JCH00{0}_true_haplopairs_STARAligned.out.sam'
##SAMBASE = '/data/James/AlleleSpecificExpression/VTP/JCH00{0}_haplopairs_STARAligned.out.sam'
SAMLIST = [SAMBASE.format('3D'),SAMBASE.format('3E'),SAMBASE.format('3F'),SAMBASE.format('3I'),SAMBASE.format('3J')]
##SAMLIST = [SAMBASE.format('3A'),SAMBASE.format('3B'),SAMBASE.format('3C'),SAMBASE.format('3G'),SAMBASE.format('3H')]
##SAMBASE = '/data/James/AlleleSpecificExpression/PAXB_BS/JCH00{0}_haplopairs_STARAligned.out.sam'
##SAMLIST = [SAMBASE.format('1D'),SAMBASE.format('1E'),SAMBASE.format('1F'),SAMBASE.format('3K'),SAMBASE.format('3L')]
SAM_COUNTS = loadSAMs(SAMLIST)
loose_SAM_COUNTS = loose_loadSAMs(SAMLIST)

SAM_PDICT,SAM_SIGDICT,SAM_ASEDICT = countStats(SAM_COUNTS,.05)
print len(SAM_COUNTS.keys())
print len(SAM_SIGDICT.keys())
##multiASE_heatmap([SAM_ASEDICT])
##writeAseBga('sim_noASE',SAM_ASEDICT,'/data/James/AlleleSpecificExpression/simVTP/simVTP_noASE.bga')

w = open('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l','w')
for ensid in SAM_SIGDICT.keys():
    ase_val = SAM_ASEDICT[ensid]
    w.write(ensid + '\t' + str(ase_val) + '\n')
    
w.flush()
w.close()

w = open('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l','w')
for ensid in SAM_ASEDICT.keys():
    ase_val = SAM_ASEDICT[ensid]
    w.write(ensid + '\t' + str(ase_val) + '\n')
    
w.flush()
w.close()

COUNT_DICT = counts_toASE(loose_SAM_COUNTS)
SORTED_ENSIDS = sorted(COUNT_DICT.keys())
SORTED_COUNT_ARRAY = []
for i in SORTED_ENSIDS:
    SORTED_COUNT_ARRAY.append(COUNT_DICT[i])
MEAN_ASE_ARRAY = collapse_array(SORTED_COUNT_ARRAY)
ASE_DICT = {}
for i,ensid in enumerate(SORTED_ENSIDS):
    ASE_DICT[ensid] = MEAN_ASE_ARRAY[i]

w = open('loose_CERC_VTP_ASE.l','w')
for i,ensid in enumerate(SORTED_ENSIDS):
    towrite = ensid + '\t' + str(ASE_DICT[ensid]) + '\n'
    w.write(towrite)
w.flush()
w.close()


'''
SAMBASE = '/data/James/AlleleSpecificExpression/VTP/CERC/JCH00{0}_haplopairs_STARAligned.out.sam'
##SAMBASE = '/data/James/AlleleSpecificExpression/VTP/JCH00{0}_haplopairs_STARAligned.out.sam'
SAMLIST = [SAMBASE.format('3D'),SAMBASE.format('3E'),SAMBASE.format('3F'),SAMBASE.format('3I'),SAMBASE.format('3J')]
##SAMLIST = [SAMBASE.format('3A'),SAMBASE.format('3B'),SAMBASE.format('3C'),SAMBASE.format('3G'),SAMBASE.format('3H')]
##SAMBASE = '/data/James/AlleleSpecificExpression/PAXB_BS/JCH00{0}_haplopairs_STARAligned.out.sam'
##SAMLIST = [SAMBASE.format('1D'),SAMBASE.format('1E'),SAMBASE.format('1F'),SAMBASE.format('3K'),SAMBASE.format('3L')]
SAM_COUNTS = loadSAMs(SAMLIST)
loose_SAM_COUNTS = loose_loadSAMs(SAMLIST)

SAM_PDICT,SAM_SIGDICT,SAM_ASEDICT = countStats(SAM_COUNTS,.05)
print len(SAM_COUNTS.keys())
print len(SAM_SIGDICT.keys())
##multiASE_heatmap([SAM_ASEDICT])
##writeAseBga('sim_noASE',SAM_ASEDICT,'/data/James/AlleleSpecificExpression/simVTP/simVTP_noASE.bga')

w = open('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l','w')
for ensid in SAM_SIGDICT.keys():
    ase_val = SAM_ASEDICT[ensid]
    w.write(ensid + '\t' + str(ase_val) + '\n')
    
w.flush()
w.close()

w = open('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l','w')
for ensid in SAM_ASEDICT.keys():
    ase_val = SAM_ASEDICT[ensid]
    w.write(ensid + '\t' + str(ase_val) + '\n')
    
w.flush()
w.close()
'''


'''
COUNT_DICT = counts_toASE(loose_SAM_COUNTS)
SORTED_ENSIDS = sorted(COUNT_DICT.keys())
SORTED_COUNT_ARRAY = []
for i in SORTED_ENSIDS:
    SORTED_COUNT_ARRAY.append(COUNT_DICT[i])
MEAN_ASE_ARRAY = collapse_array(SORTED_COUNT_ARRAY)
ASE_DICT = {}
for i,ensid in enumerate(SORTED_ENSIDS):
    ASE_DICT[ensid] = MEAN_ASE_ARRAY[i]

w = open('loose_PAXB_BS_ASE.l','w')
for i,ensid in enumerate(SORTED_ENSIDS):
    ratio = (2**ASE_DICT[ensid])/((2**ASE_DICT[ensid]) + 1)
    if ratio < .05 or ratio > .95:
	print ASE_DICT[ensid]
    towrite = ensid + '\t' + str(ratio) + '\n'
    w.write(towrite)
w.flush()
w.close()
writeAseBga('LOOSE_CERC_VTP_ASE',ASE_DICT,'loose_CERC_VTP_ASE.bga')
'''

expDiff_dict = loadASE_list('/data/James/Parental_RNA_seq/CERC_RABS_diff/PAXB_RABS_expDiff.l')

shared_keys = set(expDiff_dict.keys()).intersection(COUNT_DICT.keys())

for key in expDiff_dict.keys():
    if key not in shared_keys:
	del(expDiff_dict[key])

for key in COUNT_DICT.keys():
    if key not in shared_keys:
	del(COUNT_DICT[key])

chromes = module.stickleChrs()
for chrome in ['chrXXI']:
    SORTED_LIST = sort_CountDict(COUNT_DICT,CHR=chrome)
    collapsed_sorted_array = collapse_array(np.array(SORTED_LIST))
    SORTED_expDiff_LIST = sort_CountDict(expDiff_dict,CHR=chrome)
    plot_ase_expDiff(collapsed_sorted_array,SORTED_expDiff_LIST)
    
    plt.suptitle(chrome)
    plt.show()




COUNT_DICT = counts_toASE(loose_SAM_COUNTS)
BS_COUNT_DICT = counts_toASE(loose_BS_SAM_COUNTS)


shared_keys = set(BS_COUNT_DICT.keys()).intersection(COUNT_DICT.keys())

for key in BS_COUNT_DICT.keys():
    if key not in shared_keys:
	del(BS_COUNT_DICT[key])

for key in COUNT_DICT.keys():
    if key not in shared_keys:
	del(COUNT_DICT[key])

chromes = module.stickleChrs()
for chrome in ['chrIV','chrXXI']:
    SORTED_LIST = sort_CountDict(COUNT_DICT,CHR=chrome)
    collapsed_sorted_array = collapse_array(np.array(SORTED_LIST))
    SORTED_BS_LIST = sort_CountDict(BS_COUNT_DICT,CHR=chrome)
    collapsed_sorted_bs_array = collapse_array(np.array(SORTED_BS_LIST))
    plot_ase_expDiff(collapsed_sorted_array,collapsed_sorted_bs_array)
    
    plt.suptitle(chrome)
    plt.show()

RANDOM_LIST = SORTED_LIST[:]
random.shuffle(RANDOM_LIST)
##SORTED_LIST = RANDOM_LIST
random_sorted_array = collapse_array(np.array(RANDOM_LIST))
ASE_VALS = ASE_DICT.values()
random_sorted_array = random.sample(ASE_VALS,len(collapsed_sorted_array))
plt.subplot(122)
plot_updown(random_sorted_array)
plt.tight_layout()
plt.show()







'''

SAMBASE = '/data/James/AlleleSpecificExpression/VTP/JCH00{0}_haplopairs_STARAligned.out.sam'
SAMLIST = [SAMBASE.format('3A'),SAMBASE.format('3B'),SAMBASE.format('3C'),SAMBASE.format('3G'),SAMBASE.format('3H')]

SAM_COUNTS = loadSAMs(SAMLIST)
loose_SAM_COUNTS = loose_loadSAMs(SAMLIST)

SAMBASE = '/data/James/AlleleSpecificExpression/PAXB_BS/JCH00{0}_haplopairs_STARAligned.out.sam'
SAMLIST = [SAMBASE.format('1D'),SAMBASE.format('1E'),SAMBASE.format('1F'),SAMBASE.format('3K'),SAMBASE.format('3L')]
loose_BS_SAM_COUNTS = loose_loadSAMs(SAMLIST)


COUNT_DICT = counts_toASE(loose_SAM_COUNTS)
BS_COUNT_DICT = counts_toASE(loose_BS_SAM_COUNTS)
SORTED_ENSIDS = sorted(COUNT_DICT.keys())
SORTED_COUNT_ARRAY = []
for i in SORTED_ENSIDS:
    SORTED_COUNT_ARRAY.append(COUNT_DICT[i])
MEAN_ASE_ARRAY = collapse_array(SORTED_COUNT_ARRAY)
ASE_DICT = {}
for i,ensid in enumerate(SORTED_ENSIDS):
    ASE_DICT[ensid] = MEAN_ASE_ARRAY[i]

w = open('loose_VTP_ASE.l','w')
for i,ensid in enumerate(SORTED_ENSIDS):
    towrite = ensid + '\t' + str(ASE_DICT[ensid]) + '\n'
    w.write(towrite)
w.flush()
w.close()
writeAseBga('LOOSE_VTP_ASE',ASE_DICT,'LOOSE_VTP_ASE.BGA')







chromes = module.stickleChrs()
chromes.remove('chrUn')
total_bias = [0,0,0,0]
total_mean = 0
for index,chrome in enumerate(chromes):

 if chrome == 'chrM':
     continue

 SORTED_LIST = sort_CountDict(COUNT_DICT,CHR=chrome)
 RANDOM_LIST = SORTED_LIST[:]
 random.shuffle(RANDOM_LIST)
 ##SORTED_LIST = RANDOM_LIST
 collapsed_sorted_array = collapse_array(np.array(SORTED_LIST))
 random_sorted_array = collapse_array(np.array(RANDOM_LIST))

 diff_distance_real = []
 diff_distance_total = []
 print chrome
 for j in range(1,5):
     curr_diff = expressionDiff_distance(collapsed_sorted_array,j)
     diff_distance_real.append(curr_diff)
 for j in range(1,len(collapsed_sorted_array)/2):
     curr_diff = expressionDiff_distance(collapsed_sorted_array,j)
     diff_distance_total.append(curr_diff)

 plt.subplot(5,5,index)
 xs = np.arange(1,5)
 mean_diff = [sum(diff_distance_total)/len(diff_distance_total)] * len(xs)
 plt.plot(xs,diff_distance_real,color='blue')
 plt.plot(xs,mean_diff,color='black',ls='dashed')
 print scipy.stats.spearmanr(diff_distance_real,xs)
 print scipy.stats.pearsonr(diff_distance_real,xs)
 print sum(diff_distance_total)/len(diff_distance_total)
 print 'first 5: ' + str(sum(diff_distance_total[:5])/5)

 m,b = np.polyfit(xs,diff_distance_real,1)

 plt.plot(xs,m*xs+b,color='red',linestyle='dashed')
 plt.title(chrome)

 for i in range(4):
     ##print diff_distance_real
     total_bias[i] += ((diff_distance_real[i] * len(SORTED_LIST)/len(COUNT_DICT.keys())))
 total_mean += (sum(diff_distance_total)/len(diff_distance_total) * len(SORTED_LIST)/len(COUNT_DICT.keys()))

print total_bias, total_mean
plt.tight_layout()
plt.show()

    

 





SORTED_LIST = sort_CountDict(ASE_DICT,CHR='chrI')	
RANDOM_LIST = SORTED_LIST[:]
random.shuffle(RANDOM_LIST)
collapsed_sorted_array = collapse_array(np.array(SORTED_LIST))
random_sorted_array = collapse_array(np.array(RANDOM_LIST))
plt.subplot(221)
plt.ylim([-1,1])
binned_mean(collapsed_sorted_array,3)
plt.subplot(222)
##plt.ylim([-1.5,1.5])
binned_correlation(collapsed_sorted_array,3,1)
plt.subplot(223)
plt.ylim([-1,1])
binned_mean(random_sorted_array,3)
plt.subplot(224)
##plt.ylim([-1.5,1.5])
binned_correlation(random_sorted_array,3,1)
plt.show()


SAMBASE = '/data/James/AlleleSpecificExpression/VTP/JCH00{0}_haplopairs_STARAligned.out.sam'
SAMLIST = [SAMBASE.format('3A'),SAMBASE.format('3B'),SAMBASE.format('3C'),SAMBASE.format('3G'),SAMBASE.format('3H')]

SAM_COUNTS = loadSAMs(SAMLIST)


ASE_DICT = counts_toASE(SAM_COUNTS)
SORTED_LIST = sort_CountDict(ASE_DICT,CHR='chrXXI')
RANDOM_LIST = SORTED_LIST[:]


##CORR_MATRIX = np.corrcoef(SORTED_LIST)
CORR_MATRIX = list_distance(SORTED_LIST,SORTED_LIST)
random.shuffle(RANDOM_LIST)
##RAND_MATRIX = np.corrcoef(RANDOM_LIST)
RAND_MATRIX = list_distance(RANDOM_LIST,RANDOM_LIST)

collapsed_sorted_array = collapse_array(np.array(SORTED_LIST))
random_sorted_array = collapse_array(np.array(RANDOM_LIST))

plt.subplot(121)
plotCorrMatrix(CORR_MATRIX)
plt.subplot(122)
plotCorrMatrix(RAND_MATRIX)
plt.tight_layout()
plt.show()

'''




__author__ = 'james'

import numpy as np
from scipy import stats
import module
import math
from matplotlib import pyplot as plt
import matplotlib
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Arial')
matplotlib.rc('text', usetex='false')
matplotlib.rc('pdf', fonttype=42)

##load ASE data into dictionary in the form ase_dict[ensid] = [[F1,M1],...,[FN,Mn]]
def load_ase_dict(infile):

    ase_dict = {}

    f = open(infile)
    f.readline()

    for line in f:
        line = line.strip()
        split = line.split('\t')
        cur_id = split[0]
        ase_reps = map(lambda x: float(x), split[1:])

        toappend = []
        for i in range(0,len(ase_reps),2):
            toappend.append([ase_reps[i],ase_reps[i+1]])

        ase_dict[cur_id] =  toappend

    return ase_dict


##returns a dict in the form count_dict[ensid] = [F_counts],[M_count]
def M_F_counts_reps(in_dict):
    M_F_dict = {}

    for ensid in in_dict.keys():
        cur_ase = in_dict[ensid]

        n_reps = len(cur_ase)
        total_log_fc = 0
        cur_M = []
        cur_F = []
        for ase_rep in cur_ase:
            cur_F.append(ase_rep[0])
            cur_M.append(ase_rep[1])
            cur_fc = ase_rep[0]/(ase_rep[1] + 1e-30)
            cur_log_fc = math.log(cur_fc + 1e-30,2)
            total_log_fc += cur_log_fc
        mean_fc = (total_log_fc / n_reps)


        if mean_fc > 3 or mean_fc < -3:
            continue
        M_F_dict[ensid] = (np.array(cur_F),np.array(cur_M))

    return M_F_dict


##pooled binomial test for significance (AIR = allelic imbalance ratio (total(F)/total(M))
##returns a dict in the form sigDict[ensid] = (AIR,pval)
##if outfile, write to outfile
def getSig_ASE(in_dict,FDR=.01,COV_CUT = 20,OUTFILE=''):

    sigDict = {}

    stickleNames = module.ENSidDict()

    pval_list = []
    ensid_list = []
    AIB_list = []

    for ensid in in_dict.keys():

        cur_F,cur_M = in_dict[ensid]
        cur_N = cur_F + cur_M

        total_F = np.sum(cur_F)
        total_M = np.sum(cur_M)

        if (total_F + total_M) < COV_CUT:
            continue
        total_N = float(total_F + total_M)
        binom_p = stats.binom_test(total_F,n=total_N)

        '''if ensid == 'ENSGACT00000017212':
            print cur_F,cur_M
            print total_F,total_M
            print total_N
            print binom_p
        '''
        zscores = ((cur_F + .5) - .5*cur_N)/ ((.25 * cur_N)**.5)

        AIB = (cur_F/cur_N).mean()
        sum_zscores = np.sum(zscores)
        ##sum variance, assume uncorrelated (null model, no ASE)
        sum_var = np.sum(.25 * cur_N)

        ##tval,pval =  stats.ttest_1samp(zscores,0)

        pval_list.append(binom_p)
        ensid_list.append(ensid)
        AIB_list.append(AIB)






    pval_array = np.array(pval_list)
    ensid_array = np.array(ensid_list)
    AIB_array = np.array(AIB_list)
    cutoff = 0.

    ##FDR by B-H
    ensid_array = ensid_array[np.argsort(pval_array)]
    AIB_array = AIB_array[np.argsort(pval_array)]
    pval_array = pval_array[np.argsort(pval_array)]

    w = open('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_ASEpval.l','w')

    for ensid in ensid_array:
        w.write(stickleNames[ensid] + '\n')

    w.flush()
    w.close()



    fdr_list = []
    for k in range(len(pval_array)):
        cutoff += 1.
        cur_fdr = (cutoff/len(pval_array)) * FDR
        ##print cutoff,len(pval_array)
        ##print cur_fdr
        ##print pval_array[k]
        if cur_fdr < pval_array[k]:
            print cutoff
            break


    if OUTFILE:
        w = open(OUTFILE,'w')

        header = ['ensid']
        for i in range(len(in_dict.values()[0][0])):
            F_rep = 'F_%d' % i
            header.append(F_rep)

        for i in range(len(in_dict.values()[0][0])):
            M_rep = 'M_%d' % i
            header.append(M_rep)
        header.append('pval')
        header.append('AIB')
        w.write('\t'.join(header) + '\n')
        for i in range(int(cutoff)):

            cur_ensid = ensid_array[i]
            cur_name = stickleNames[cur_ensid]
            cur_F,cur_M = in_dict[cur_ensid]
            cur_pval = pval_array[i]
            cur_AIB = AIB_array[i]

            printlist = [cur_name] + list(cur_F) + list(cur_M) + [cur_pval,cur_AIB]
            printlist = map(lambda x: str(x), printlist)
            w.write('\t'.join(printlist) + '\n')
        w.flush()
        w.close()

    ##make the AIB,pval dict

    for i in range(int(cutoff)):

        cur_ensid = ensid_array[i]
        cur_F,cur_M = in_dict[cur_ensid]
        cur_pval = pval_array[i]
        cur_AIB = AIB_array[i]
        sigDict[cur_ensid] = (cur_AIB,cur_pval)

    return sigDict



##correlates two F_M read count dictionaries
def cor_2dicts(key_list,dict1,dict2,OUTBASE=''):

    stickleNames = module.ENSid_Human_Dict()

    aib1 = []
    aib2 = []
    samedir = 0
    diffdir = 0
    samedir_set = set()
    diffdir_set = set()

    for ensid in key_list:

        if type(dict1[ensid]) == tuple:
            F1,M1 = dict1[ensid]
            F2,M2 = dict2[ensid]

            ##cur_aib1 = float(F1.sum())/(F1.sum() + M1.sum())
            ##cur_aib2 = float(F2.sum())/(F2.sum() + M2.sum())
            cur_aib1 = math.log(float(F1.sum())/M1.sum(),2)
            cur_aib2 = math.log(float(F2.sum())/M2.sum(),2)
        else:
            cur_aib1 = dict1[ensid]
            cur_aib2 = dict2[ensid]
        ##print cur_aib1,cur_aib2
        aib1.append(cur_aib1)
        aib2.append(cur_aib2)

        ##if (cur_aib1 > .5 and cur_aib2 > .5) or (cur_aib1 < .5 and cur_aib2 < .5):
        if (cur_aib1 > 0 and cur_aib2 > 0) or (cur_aib1 < 0 and cur_aib2 < 0):
            ##print stickleNames[ensid]
            samedir += 1
            samedir_set.add(ensid)
        else:
            ##
            diffdir += 1
            diffdir_set.add(ensid)

    aib1 = np.array(aib1)
    aib2 = np.array(aib2)

    print stats.pearsonr(aib1,aib2)
    print samedir,diffdir
    return samedir_set,diffdir_set


##laods D stat for each ensid into a dictionary
##D_dict[ensid] = D_val
def load_D_bygene(inbed):

    f = open(inbed,'r')

    D_dict = {}

    for line in f:
        line = line.strip()
        split = line.split('\t')

        cur_ensid = split[3]
        D_stat = float(split[4])

        D_dict[cur_ensid] = D_stat


    return D_dict

##takes an array and a sample_name
##plots the smoothed histogram (density) of the data
def plotSmoothedHist(array,sample_name,x_range=[-2,8]):

    density = stats.gaussian_kde(array)

    xmin = x_range[0]
    xmax = x_range[1]
    print xmin,xmax
    xs = np.linspace(xmin,xmax,200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    plt.plot(xs,density(xs),label=sample_name)

##takes a list of gene_sets, and a list of names for each gene set
##plots the density of each
def plot_D_sets(gene_set_list,sample_name_list,D_VAL_FILE='/home/james/Dropbox/Miller/data/Introgression/ensGene_collapsed_TSS_Dstat_20kb.bed'):

    D_vals_all = load_D_bygene(D_VAL_FILE)
    array_list = []

    for gene_set,sample in zip(gene_set_list,sample_name_list):

        cur_d_list = []

        for gene in gene_set:
            if gene in D_vals_all:
                cur_d_list.append(D_vals_all[gene])


        ##cur_d_list = [D_vals_all[x] for x in gene_set]
        cur_d_array = np.array(cur_d_list)
        print cur_d_array.mean()
        array_list.append(cur_d_array)
        plotSmoothedHist(cur_d_array,sample,x_range=[-1,1])

    print stats.mannwhitneyu(array_list[0],array_list[1])
    ##print stats.mannwhitneyu(array_list[0],array_list[2])
    ##print stats.mannwhitneyu(array_list[1],array_list[2])
    plt.legend()
    ##plt.show()

##load start positions for each ensid from a given BED
##returns a dict in the form start_dict[ensid] = (chrome,pos) (0-based)
##returns the last coordinate for - strand genes
def load_starts(inbed='/home/james/Dropbox/Miller/data/BED/ensGene_collapsed.bed'):

    start_dict = {}

    f = open(inbed)

    for line in f:
        line = line.strip()
        split = line.split('\t')

        chrome = split[0]
        start = split[1]
        end = split[2]
        ensid = split[3]
        strand = split[5]

        if strand == '-':
            start_dict[ensid] = (chrome,end)
        else:
            start_dict[ensid] = (chrome,start)

    return start_dict

##takes an iterable with ensids
def write_starts_bed(id_list,outfile):

    starts = load_starts()

    w = open(outfile,'w')

    for ensid in id_list:

        chrome, cur_start = starts[ensid]

        cur_stop = int(cur_start) + 1

        w.write('\t'.join([chrome,str(cur_start),str(cur_stop),ensid]) + '\n')

    w.flush()
    w.close()


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
def distanceHist(filea,fileb,filec=''):

    a_dist = []
    b_dist = []

    a = open(filea,'r')

    for line in a:
        line = line.strip()
        split = line.split('\t')
        distance = int(split[-1])
        if distance < 0:
            continue

        a_dist.append(distance)

    b = open(fileb,'r')

    for line in b:
        line = line.strip()
        split = line.split('\t')

        distance = int(split[-1])

        if distance < 0:
            continue
        b_dist.append(distance)

    a_dist = np.array(a_dist)
    b_dist = np.array(b_dist)

    plotSmoothedHist(a_dist,'DUUD',x_range=[0,10000000])
    plotSmoothedHist(b_dist,'UUDD',x_range=[0,10000000])

    if filec:

        c_dist = []

        c = open(filec,'r')

        for line in c:
            line = line.strip()
            split = line.split('\t')

            distance = int(split[-1])

            if distance < 0:
                continue
            c_dist.append(distance)

        c_dist = np.array(c_dist)


        plotSmoothedHist(c_dist,'genes',x_range=[0,10000000])

    print a_dist.mean(),b_dist.mean(),c_dist.mean()
    print stats.f_oneway(a_dist,b_dist)
    print stats.mannwhitneyu(a_dist,b_dist)
    print stats.f_oneway(a_dist,c_dist)
    print stats.mannwhitneyu(a_dist,c_dist)
    print stats.f_oneway(c_dist,b_dist)
    print stats.mannwhitneyu(c_dist,b_dist)

    plt.legend()
    ##plt.show()

##returns the set of genes within the 1,11,21 inversions
def load_inv_genes(inbed):
    toret = set()

    f = open(inbed)

    for line in f:
        line = line.strip()
        split = line.split('\t')
        ensid = split[-1]

        toret.add(ensid)

    return toret


##loads a file into a dictionary
##if log > 0, log transform using the given base
def load_logfc(infile,LOG=0):

    toret = {}

    f = open(infile)

    for line in f:
        line = line.strip()
        split = line.split('\t')
        cur_id = split[0]
        cur_fc = float(split[1])
        if LOG > 0:
            if cur_fc < .01:
                continue
            #print cur_fc,LOG
            cur_fc = math.log(cur_fc,LOG)
        toret[cur_id] = cur_fc

    return toret


##takes two dicts - parental fc and hybrid AIB
##returns the calculated trans values
def calc_trans(parental,hybrid,outfile=''):


    shared_ids = set(hybrid.keys()).intersection(set(parental.keys()))

    ##stores abs value of cis and trans exp changes

    trans_vals = {}
    percent_cis_vals = {}

    opposing_cis = 0
    opposing_trans = 0
    supportive_cis = 0
    supportive_trans = 0

    for ensid in shared_ids:
        cis_trans = parental[ensid]
        cis = hybrid[ensid]

        trans = cis_trans - cis

        percent_cis = math.fabs(cis)/(math.fabs(cis) + math.fabs(trans))##math.fabs(cis)/cis_trans##
        ##percent_cis = cis/cis_trans


        trans_vals[ensid] = trans
        percent_cis_vals[ensid] = percent_cis

        if (cis > 0 and trans > 0) or (cis < 0 and trans < 0):
            if math.fabs(cis) > math.fabs(trans):
                supportive_cis += 1
            else:
                supportive_trans += 1
        else:
            if math.fabs(cis) > math.fabs(trans):
                opposing_cis += 1
            else:
                opposing_trans += 1

    print opposing_cis,opposing_trans,supportive_cis,supportive_trans



    return trans_vals,percent_cis_vals


def load_chr_size(GFILE='/home/james/Dropbox/Miller/data/ReferenceGenomes/stickleback.genome'):

    f = open(GFILE)

    genome_dict = {}

    for line in f:
        line = line.strip()
        split = line.split()
        chr = split[0]
        size = int(split[1])

        genome_dict[chr]=size

    return genome_dict



##writes a bed with the mean logfc for each gene
##size specifies the window
def write_logFC_bga(indict,outbed,SIZE=2000,NAME='',DESCRIPTION=''):

    TSS_starts = load_starts()

    w = open(outbed,'w')

    gsize = load_chr_size()

    if not DESCRIPTION:
        DESCRIPTION = NAME

    if NAME:
        w.write('track type=bedGraph name="%s" description="%s"' % (NAME,DESCRIPTION) + '\n')

    for gene in indict.keys():

        if gene not in TSS_starts:
            continue

        cur_lfc = indict[gene]


        cur_chr,cur_start = TSS_starts[gene]

        cur_size = gsize[cur_chr]


        cur_stop = str(int(cur_start) + SIZE)


        if int(cur_stop) >= cur_size:
            cur_stop = str(cur_size-1)


        w.write(cur_chr + '\t' + str(cur_start) + '\t' + cur_stop + '\t' + str(cur_lfc) + '\n')

    w.flush()
    w.close()






fpath_base = '/home/james/Dropbox/Miller/data/ASE/'
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE/'
'''
distanceHist(fpath_base+'PAXB_CERC_all_ase_starts_diffDir_strict.bed',fpath_base+'PAXB_CERC_all_ase_starts_sameDir_strict.bed',filec=fpath_base+'PAXB_CERC_sig_ase_starts_strict.bed')
plt.show()
distanceHist(fpath_base+'PAXB_CERC_all_ase_starts_diffDir_union.bed',fpath_base+'PAXB_CERC_all_ase_starts_sameDir_union.bed',filec=fpath_base+'PAXB_CERC_sig_ase_starts_union.bed')
plt.show()
'''



PAXB_VTP_expdiff = load_logfc(fpath_base + 'PAXB_RABS_logfc.l')
PAXB_VTP_expdiff_sig = load_logfc(fpath_base + 'PAXB_RABS_sig_logfc.l')
PAXB_VTP_ase = load_logfc(fpath_base + 'PAXB_VTP_logAI.l')
PAXB_VTP_sig_ase = load_logfc(fpath_base + 'PAXB_VTP_logAI_sig.l')

CERC_VTP_expdiff = load_logfc(fpath_base + 'CERC_RABS_logfc.l')
CERC_VTP_expdiff_sig = load_logfc(fpath_base + 'CERC_RABS_sig_logfc.l')
CERC_VTP_ase = load_logfc(fpath_base + 'CERC_VTP_logAI.l')
CERC_VTP_sig_ase = load_logfc(fpath_base + 'CERC_VTP_logAI_sig.l')

PAXB_BS_ase = load_logfc(fpath_base + 'PAXB_BS_logAI.l')
PAXB_BS_sig_ase = load_logfc(fpath_base + 'PAXB_BS_logAI_sig.l')

PAXB_trans,PAXB_percent_cis = calc_trans(PAXB_VTP_expdiff,PAXB_VTP_ase)
CERC_trans,CERC_percent_cis = calc_trans(CERC_VTP_expdiff,CERC_VTP_ase)


write_logFC_bga(PAXB_VTP_ase,fpath_base + 'PAXB_VTP_ASE.bga',NAME='PAXB_VTP_AIB')
write_logFC_bga(PAXB_BS_ase,fpath_base + 'PAXB_BS_ASE.bga',NAME='PAXB_BS_AIB')
write_logFC_bga(CERC_VTP_ase,fpath_base + 'CERC_VTP_ASE.bga',NAME='CERC_VTP_AIB')


write_logFC_bga(PAXB_VTP_expdiff,fpath_base + 'PAXB_VTP_expDiff.bga',NAME='PAXB_VTP_expDiff')
write_logFC_bga(CERC_VTP_expdiff,fpath_base + 'CERC_VTP_expDiff.bga',NAME='CERC_VTP_expDiff')

'''
##Get PAXB_VTP
PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_VTP_M_F_dict = M_F_counts_reps(PAXB_VTP_ase_dict)
PAXB_VTP_sig_ase_dict = getSig_ASE(PAXB_VTP_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_ASEpval.tsv')

##Get PAXB BS
PAXB_BS_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_AI_all.tsv')
PAXB_BS_M_F_dict = M_F_counts_reps(PAXB_BS_ase_dict)
PAXB_BS_sig_ase_dict = getSig_ASE(PAXB_BS_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_BS_ASEpval.tsv')

##Get CERC VTP
CERC_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_AI_all.tsv')
CERC_VTP_M_F_dict = M_F_counts_reps(CERC_VTP_ase_dict)
CERC_VTP_sig_ase_dict = getSig_ASE(CERC_VTP_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/CERC_VTP_ASEpval.tsv')

PAXB_VTP_sig_genes = set(PAXB_VTP_sig_ase_dict.keys())
PAXB_BS_sig_genes = set(PAXB_BS_sig_ase_dict.keys())
CERC_VTP_sig_genes = set(CERC_VTP_sig_ase_dict.keys())

PAXB_VTP_genes = set(PAXB_VTP_M_F_dict)
PAXB_BS_genes = set(PAXB_BS_M_F_dict)
CERC_VTP_genes = set(CERC_VTP_M_F_dict)

'''
PAXB_VTP_sig_genes = set(PAXB_VTP_sig_ase.keys())
PAXB_BS_sig_genes = set(PAXB_BS_sig_ase.keys())
CERC_VTP_sig_genes = set(CERC_VTP_sig_ase.keys())

PAXB_VTP_genes = set(PAXB_VTP_ase.keys())
PAXB_BS_genes = set(PAXB_BS_ase.keys())
CERC_VTP_genes = set(CERC_VTP_ase.keys())

##get shared ASE ids
shared_VTP_ids = PAXB_VTP_genes.intersection(CERC_VTP_genes)
shared_sig_VTP_ids = PAXB_VTP_sig_genes.intersection(CERC_VTP_sig_genes)
shared_all_ids = CERC_VTP_genes.intersection(PAXB_BS_genes).intersection(PAXB_VTP_genes)
shared_all_sig_ids = CERC_VTP_sig_genes.intersection(PAXB_BS_sig_genes).intersection(PAXB_VTP_sig_genes)

##get shared expDiff ids
shared_VTP_expdiff_ids = set(PAXB_VTP_expdiff.keys()).intersection(set(CERC_VTP_expdiff.keys()))
shared_VTP_expdiff_sig_ids = set(PAXB_VTP_expdiff_sig.keys()).intersection(set(CERC_VTP_expdiff_sig.keys()))
##sig trans are sig ExpDiff but not sigASE
shared_VTP_trans_ids = set(PAXB_trans.keys()).intersection(CERC_trans.keys())

##correlate all dicts, get sameDir ASE, Trans, and ExpDiff
VTPsamedir_sig,VTPdiffdir_sig = cor_2dicts(shared_sig_VTP_ids,PAXB_VTP_ase,CERC_VTP_ase)
VTPsamedir,VTPdiffdir = cor_2dicts(shared_VTP_ids,PAXB_VTP_ase,CERC_VTP_ase)


VTP_trans_samedir,VTP_trans_diffdir = cor_2dicts(shared_VTP_trans_ids,PAXB_trans,CERC_trans)

VTP_expdiff_samedir_sig,VTP_expdiff_diffdir_sig = cor_2dicts(shared_VTP_expdiff_sig_ids,PAXB_VTP_expdiff,CERC_VTP_expdiff)


write_starts_bed(shared_VTP_ids,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_all_ase_starts.bed')
write_starts_bed(shared_sig_VTP_ids,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_sig_ase_starts.bed')
write_starts_bed(VTPsamedir_sig,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_sig_ase_starts_sameDir.bed')
write_starts_bed(VTPdiffdir_sig,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_sig_ase_starts_diffDir.bed')
write_starts_bed(VTPsamedir,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_all_ase_starts_sameDir.bed')
write_starts_bed(VTPdiffdir,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_all_ase_starts_diffDir.bed')



plt.figure(figsize=(10,10))
plot_D_sets([VTPsamedir_sig,VTPdiffdir_sig],['Concordant','Discordant']
            ,D_VAL_FILE='/home/james/Dropbox/Miller/data/Variants/Introgression/ensGene_collapsed_TSS_TAP_20kb.bed')
plt.xlabel('D_Score')
plt.ylabel('Density')
plt.title('shared Cis')
plt.xlim([0,1])
plt.savefig(fig_base + 'Shared_ASE_introgression_TAP.pdf')
plt.show()

plt.figure(figsize=(10,10))
plot_D_sets([VTP_trans_samedir,VTP_trans_diffdir],['Concordant','Discordant']
,D_VAL_FILE='/home/james/Dropbox/Miller/data/Variants/Introgression/ensGene_collapsed_TSS_TAP_20kb.bed')
plt.xlabel('D_Score')
plt.ylabel('Density')
plt.title('shared Trans')
plt.xlim([0,1])
plt.savefig(fig_base + 'Shared_trans_introgression_TAP.pdf')
plt.show()

plt.figure(figsize=(10,10))
plot_D_sets([VTP_expdiff_samedir_sig,VTP_expdiff_diffdir_sig],['Concordant','Discordant']
,D_VAL_FILE='/home/james/Dropbox/Miller/data/Variants/Introgression/ensGene_collapsed_TSS_TAP_20kb.bed')
plt.xlabel('D_Score')
plt.ylabel('Density')
plt.xlim([0,1])
plt.savefig(fig_base + 'Shared_expDiff_introgression_TAP.pdf')
plt.show()


'''
inv_genes = load_inv_genes('/home/james/Dropbox/Miller/data/BED/inversions_ensGene_collapsed_TSS.bed')


plt.figure(figsize=(10,10))
plot_D_sets([VTPsamedir_sig.difference(inv_genes),shared_all_ids.difference(inv_genes)],['same_direction','all'])
plt.xlabel('D_Score')
plt.ylabel('Density')
##plt.savefig(fig_base + 'Shared_ASE_introgression_D.pdf')
##plt.savefig(fig_base + 'Shared_ASE_introgression_D.svg')
plt.show()
'''
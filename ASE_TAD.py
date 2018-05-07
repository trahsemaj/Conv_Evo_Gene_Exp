__author__ = 'james'

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import module
import math
import random

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





##transform ASE dict into log-fold-change ase_dict[ensid] = logFC( F vs M)
def logfc_ase_dict(in_dict):

    logfc_dict = {}

    for ensid in in_dict.keys():
        cur_ase = in_dict[ensid]

        n_reps = len(cur_ase)
        total_log_fc = 0
        for ase_rep in cur_ase:
            cur_fc = ase_rep[0]/(ase_rep[1] + 1e-30)
            cur_log_fc = math.log(cur_fc + 1e-30,2)
            total_log_fc += cur_log_fc
        mean_fc = (total_log_fc / n_reps)


        if mean_fc > 3 or mean_fc < -3:
            continue
        logfc_dict[ensid] = mean_fc

    return logfc_dict

##returns a dict in the form count_dict[ensid] = [F_counts],[M_count]
def M_F_counts_reps(in_dict):
    M_F_dict = {}
    count = 0
    low_cov = 0
    for ensid in in_dict.keys():
        cur_ase = in_dict[ensid]

        n_reps = len(cur_ase)
        total_log_fc = 0
        cur_M = []
        cur_F = []
        for ase_rep in cur_ase:
            cur_F.append(ase_rep[0])
            cur_M.append(ase_rep[1])
            #cur_fc = ase_rep[0]/(ase_rep[1] + 1e-30)
            #cur_log_fc = math.log(cur_fc + 1e-30,2)
            #total_log_fc += cur_log_fc
        #mean_fc = (total_log_fc / n_reps)
        mean_fc = math.log(float(sum(cur_F)+ 1e-20)/(sum(cur_M)+ 1e-20),2)

        if sum(cur_M) + sum(cur_F) < 20:
            low_cov += 1
            continue
        if mean_fc > 5 or mean_fc < -5:
            count += 1
            print cur_M,cur_F,ensid
            continue

        M_F_dict[ensid] = (np.array(cur_F),np.array(cur_M))
    print 'count',count,low_cov
    return M_F_dict

##returns a set of significant ensids
def getSigGenes(gene_exp_diff,condition1,condition2):

    sigset = set()
    total = 0

    f = open(gene_exp_diff,'r')

    ##chomp dat header
    f.readline()

    for line in f:

        line = line.strip()
        split = line.split('\t')

        ##grab the first ensid in the XLOC
        ensids = split[2].split(',')

        ##only look at the given samples
        if (condition1,condition2) != (split[4],split[5]):
            continue

        exp1 = float(split[7]) + 1e-20
        exp2 = float(split[8]) + 1e-20
        total +=1

        ##if the expressions are significantly different
        if split[-1] == 'yes':
            for ensid in ensids:
                sigset.add(ensid)

    return sigset
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

def neighbor_correlation(in_dict,start_dict):

    ##gen neighbors
    start_byPos = {}

    neighbor_cor = []
    total_cor = []

    for ensid in start_dict.keys():
        if ensid not in in_dict:
            continue
        cur_chrome,cur_start = start_dict[ensid]
        if cur_chrome not in start_byPos:
            start_byPos[cur_chrome] = {}

        start_byPos[cur_chrome][cur_start] = ensid

    ##for each chromosomes
    for chrome in start_byPos.keys():
        ##for each position in order
        sorted_pos = sorted(start_byPos[chrome].keys())
        for i,positions in enumerate(sorted_pos):
            ##print positions
            ##if the first gene on the chromosome
            if i == 0:
                continue
            ##if the last on the chromosome
            elif i == (len(sorted_pos)-1):
                continue
            else:
                back_id = start_byPos[chrome][sorted_pos[i-1]]
                cur_id = start_byPos[chrome][sorted_pos[i]]
                fw_id = start_byPos[chrome][sorted_pos[i+1]]

                back_F,back_M = in_dict[back_id]
                cur_F,cur_M = in_dict[cur_id]
                fw_F,fw_M = in_dict[fw_id]


                back_cor,p = stats.pearsonr(back_F/back_M,cur_F/cur_M)
                fw_cor,p = stats.pearsonr(fw_F/fw_M,cur_F/cur_M)

                ##neighbor_cor.append(back_cor)
                ##neighbor_cor.append(fw_cor)
                neighbor_cor.append(max(fw_cor,back_cor))


    ##get total corr
    for c_ensid in in_dict.keys():

        c_F,c_M = in_dict[c_ensid]
        allowed_ensids = in_dict.keys()
        allowed_ensids.remove(c_ensid)
        [b_ensid,f_ensid] = random.sample(allowed_ensids,2)
        b_F,b_M = in_dict[b_ensid]
        f_F,f_M = in_dict[f_ensid]


        b_cor,p = stats.pearsonr(b_F/b_M,c_F/c_M)
        f_cor,p = stats.pearsonr(f_F/f_M,c_F/c_M)

        total_cor.append(max(b_cor,f_cor))

    neighbor_cor = np.array(neighbor_cor)
    total_cor = np.array(total_cor)
    neighbor_cor = neighbor_cor[~np.isnan(neighbor_cor)]
    total_cor = total_cor[~np.isnan(total_cor)]
    return neighbor_cor,total_cor





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
        ##t,binom_p = stats.ttest_rel(cur_F,cur_M)

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

    for ensid in key_list:
        F1,M1 = dict1[ensid]
        F2,M2 = dict2[ensid]

        ##cur_aib1 = float(F1.sum())/(F1.sum() + M1.sum())
        ##cur_aib2 = float(F2.sum())/(F2.sum() + M2.sum())
        cur_aib1 = math.log(float(F1.sum())/M1.sum(),2)
        cur_aib2 = math.log(float(F2.sum())/M2.sum(),2)
        ##print cur_aib1,cur_aib2
        aib1.append(cur_aib1)
        aib2.append(cur_aib2)

        ##if (cur_aib1 > .5 and cur_aib2 > .5) or (cur_aib1 < .5 and cur_aib2 < .5):
        if (cur_aib1 > 0 and cur_aib2 > 0) or (cur_aib1 < 0 and cur_aib2 < 0):
            ##print stickleNames[ensid]
            samedir += 1
        else:
            ##
            diffdir += 1

    aib1 = np.array(aib1)
    aib2 = np.array(aib2)

    print stats.pearsonr(aib1,aib2)
    print samedir,diffdir

##writes a dictionary to a file using the form ensid    log(fc)
def write_ase_file(M_F_dict,outfile,NAMES=False,sig_dict=0):

    w=open(outfile,'w')
    stickleNames = module.ENSid_Human_Dict()


    for ensid in M_F_dict.keys():
        cur_F,cur_M = M_F_dict[ensid]

        cur_aib = math.log(float(cur_F.sum())/cur_M.sum(),2)

        if sig_dict:
            if ensid not in sig_dict:
                continue

        if NAMES:
            ensid = stickleNames[ensid]

        w.write(ensid + '\t' + str(cur_aib) + '\n')

    w.flush()
    w.close()

##writes the mean logfc of each ensid to a tsv file
##ensid mean(logfc)
def write_expDiff_file(in_fpkm,test_index,ref_index,outfile,sig_dict=0):

    f = open(in_fpkm)
    f.readline()

    w = open(outfile,'w')

    for line in f:
        line = line.strip()
        split = line.split('\t')
        cur_id = split[0]
        cur_exps = np.array(map(lambda x: float(x), split[1:]))
        test_exps = cur_exps[test_index]
        ref_exps = cur_exps[ref_index]

        if cur_exps.mean() < .5:
            continue

        if test_exps.mean() < .1 or ref_exps.mean() < .1:
            continue

        cur_logfc = math.log(test_exps.mean()/ref_exps.mean(),2)

        if sig_dict:
            if cur_id not in sig_dict:
                continue

        w.write(cur_id + '\t' + str(cur_logfc) + '\n')

    w.flush()
    w.close()

##writes the mean logfc of each ensid to a tsv file
##ensid mean(logfc)
def write_expDiff_mean(in_fpkm,outfile,sig_dict=0):

    f = open(in_fpkm)
    f.readline()

    w = open(outfile,'w')

    for line in f:
        line = line.strip()
        split = line.split('\t')
        cur_id = split[0]
        cur_exps = np.array(map(lambda x: float(x), split[1:]))


        if cur_exps.mean() < .5:
            continue


        if sig_dict:
            if cur_id not in sig_dict:
                continue

        w.write(cur_id + '\t' + str(cur_exps.mean()) + '\n')

    w.flush()
    w.close()


def plot_ase(ase_dict,start_dict,chrome):

    ensid_toplot = []
    ase_toplot = []
    pos_toplot = []
    print ase_dict.keys()

    for ensid in start_dict.keys():

        cur_chrome,cur_pos = start_dict[ensid]

        if cur_chrome != chrome:
            continue
        if ensid not in ase_dict.keys():
            continue

        ensid_toplot.append(ensid)
        ase_toplot.append(ase_dict[ensid])
        pos_toplot.append(cur_pos)

    ase_toplot = np.array(ase_toplot)
    pos_toplot = np.array(pos_toplot)
    ensid_toplot = np.array(ensid_toplot)


    ##sort everything by postion
    ase_toplot = ase_toplot[pos_toplot.argsort()]
    pos_toplot = pos_toplot[pos_toplot.argsort()]
    ensid_toplot = ensid_toplot[pos_toplot.argsort()]
    '''
    plt.figure()
    plt.scatter(pos_toplot,ase_toplot)
    plt.plot(pos_toplot,np.zeros(len(pos_toplot)))
    plt.ylim([-3,3])

    '''

    ase_resize = np.resize(ase_toplot,(1,len(ase_toplot)))
    ##ase_diff = np.maximum(1-(ase_resize - ase_resize.T)**2,np.zeros( (len(ase_toplot),len(ase_toplot))))
    ase_diff = np.sign(ase_resize) * np.sign(ase_resize).T * np.maximum(1-(ase_resize - ase_resize.T)**2,np.zeros( (len(ase_toplot),len(ase_toplot))))
    np.random.shuffle(ase_toplot)
    ase_resize = np.resize(ase_toplot,(1,len(ase_toplot)))
    ##rand_ase_diff = np.maximum(1-(ase_resize - ase_resize.T)**2,np.zeros( (len(ase_toplot),len(ase_toplot))))
    rand_ase_diff = np.sign(ase_resize) * np.sign(ase_resize).T * np.maximum(1-(ase_resize - ase_resize.T)**2,np.zeros( (len(ase_toplot),len(ase_toplot))))


    plt.figure()
    plt.imshow(ase_diff)
    plt.colorbar(orientation='vertical')
    plt.title('real')

    plt.figure()
    plt.imshow(rand_ase_diff)
    plt.colorbar(orientation='vertical')
    plt.title('random')
    plt.show()

write_expDiff_mean('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table','/home/james/Dropbox/Miller/data/ASE/CRPL_meanExp.l')

ensNames = module.ENSidDict()
start_dict = load_starts()
PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_VTP_M_F_dict = M_F_counts_reps(PAXB_VTP_ase_dict)

CERC_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_AI_all.tsv')
CERC_VTP_M_F_dict = M_F_counts_reps(CERC_VTP_ase_dict)


'''
write_expDiff_mean('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table','/home/james/Dropbox/Miller/data/ASE/CRPL_meanExp.l')

ensNames = module.ENSidDict()
start_dict = load_starts()
PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_VTP_M_F_dict = M_F_counts_reps(PAXB_VTP_ase_dict)
PAXB_VTP_sig_ase_dict = getSig_ASE(PAXB_VTP_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_ASEpval_norm.tsv')
print len(PAXB_VTP_sig_ase_dict.keys())

FTC_BS_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/FTC_BS_AI_all.tsv')
FTC_BS_M_F_dict = M_F_counts_reps(FTC_BS_ase_dict)
write_ase_file(FTC_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/FTC_BS_logAI.l')
write_ase_file(FTC_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/FTC_BS_logAI_names.l',NAMES=True)
FTC_BS_sig_ase_dict = getSig_ASE(FTC_BS_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/FTC_BS_ASEpval.tsv')
write_ase_file(FTC_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/FTC_BS_logAI_sig.l',sig_dict=FTC_BS_sig_ase_dict)
write_ase_file(FTC_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/FTC_BS_logAI_sig_names.l',sig_dict=FTC_BS_sig_ase_dict,NAMES=True)

print 'PAXB_VTP'

##Get PAXB_VTP
PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_VTP_M_F_dict = M_F_counts_reps(PAXB_VTP_ase_dict)
print len(PAXB_VTP_ase_dict)
print len(PAXB_VTP_M_F_dict)
write_ase_file(PAXB_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_logAI.l')
write_ase_file(PAXB_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_logAI_names.l',NAMES=True)
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[6,7,8],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/PAXB_RABS_logfc.l')
sig_PAXB_RABS = getSigGenes('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_gene_exp.diff','RABS','PAXB')
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[6,7,8],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/PAXB_RABS_sig_logfc.l',sig_dict=sig_PAXB_RABS)



PAXB_VTP_sig_ase_dict = getSig_ASE(PAXB_VTP_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_ASEpval.tsv')
write_ase_file(PAXB_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_logAI_sig.l',sig_dict=PAXB_VTP_sig_ase_dict)
write_ase_file(PAXB_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_logAI_sig_names.l',sig_dict=PAXB_VTP_sig_ase_dict,NAMES=True)

##Get PAXB BS
PAXB_BS_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_AI_all.tsv')
PAXB_BS_M_F_dict = M_F_counts_reps(PAXB_BS_ase_dict)
write_ase_file(PAXB_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_logAI.l')
write_ase_file(PAXB_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_logAI_names.l',NAMES=True)
PAXB_BS_sig_ase_dict = getSig_ASE(PAXB_BS_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_BS_ASEpval.tsv')
write_ase_file(PAXB_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_logAI_sig.l',sig_dict=PAXB_BS_sig_ase_dict)
write_ase_file(PAXB_BS_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_logAI_sig_names.l',sig_dict=PAXB_BS_sig_ase_dict,NAMES=True)

##Get PAXB_BS, with D_ref only (all aligned to D personal transcriptome)
PAXB_BS_Dref_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_AI_Dref_all.tsv')
PAXB_BS_Dref_M_F_dict = M_F_counts_reps(PAXB_BS_Dref_ase_dict)
write_ase_file(PAXB_BS_Dref_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_Dref_logAI.l')
PAXB_BS_Dref_sig_ase_dict = getSig_ASE(PAXB_BS_Dref_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/PAXB_BS_Dref_ASEpval.tsv')
write_ase_file(PAXB_BS_Dref_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/PAXB_BS_Dref_logAI_sig.l',sig_dict=PAXB_BS_Dref_sig_ase_dict)

##Get CERC VTP
CERC_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_AI_all.tsv')
CERC_VTP_M_F_dict = M_F_counts_reps(CERC_VTP_ase_dict)
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[0,1,2],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/CERC_RABS_logfc.l')
write_ase_file(CERC_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/CERC_VTP_logAI.l')
write_ase_file(CERC_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/CERC_VTP_logAI_names.l',NAMES=True)
sig_CERC_RABS = getSigGenes('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_gene_exp.diff','CERC','RABS')
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[0,1,2],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/CERC_RABS_sig_logfc.l',sig_dict=sig_CERC_RABS)

CERC_VTP_sig_ase_dict = getSig_ASE(CERC_VTP_M_F_dict,OUTFILE='/home/james/Dropbox/Miller/data/ASE/CERC_VTP_ASEpval.tsv')
write_ase_file(CERC_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/CERC_VTP_logAI_sig.l',sig_dict=CERC_VTP_sig_ase_dict)
write_ase_file(CERC_VTP_M_F_dict,'/home/james/Dropbox/Miller/data/ASE/CERC_VTP_logAI_sig_names.l',sig_dict=CERC_VTP_sig_ase_dict,NAMES=True)


write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[6,7,8],[0,1,2],'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_logfc.l')
sig_PAXB_CERC = getSigGenes('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_gene_exp.diff','CERC','PAXB')
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids.fpkm_table',[6,7,8],[0,1,2],'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_sig_logfc.l',sig_dict=sig_PAXB_CERC)


write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids_bc.fpkm_table',[6,7,8],[0,1,2],'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_bc_logfc.l')
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids_bc.fpkm_table',[6,7,8],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/PAXB_RABS_bc_logfc.l')
write_expDiff_file('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_ids_bc.fpkm_table',[0,1,2],[3,4,5],'/home/james/Dropbox/Miller/data/ASE/CERC_RABS_bc_logfc.l')
'''
'''
PAXB_VTP_sig_genes = set(PAXB_VTP_sig_ase_dict.keys())
PAXB_BS_sig_genes = set(PAXB_BS_sig_ase_dict.keys())
CERC_VTP_sig_genes = set(CERC_VTP_sig_ase_dict.keys())

PAXB_VTP_genes = set(PAXB_VTP_M_F_dict)
PAXB_BS_genes = set(PAXB_BS_M_F_dict)
CERC_VTP_genes = set(CERC_VTP_M_F_dict)

print 'PAXB_VTP',len(PAXB_VTP_M_F_dict),len(PAXB_VTP_sig_genes)
print 'PAXB_BS',len(PAXB_BS_M_F_dict),len(PAXB_BS_sig_genes)
print 'CERC_VTP',len(CERC_VTP_M_F_dict),len(CERC_VTP_sig_genes)

print 'tPvvPb', len(PAXB_VTP_genes.intersection(PAXB_BS_genes))
print 'tPvvCv', len(PAXB_VTP_genes.intersection(CERC_VTP_genes))
print 'tCvvPb', len(CERC_VTP_genes.intersection(PAXB_BS_genes))
print 'tall', len(CERC_VTP_genes.intersection(PAXB_BS_genes).intersection(PAXB_VTP_genes))

print 'sPvvPb', len(PAXB_VTP_sig_genes.intersection(PAXB_BS_sig_genes))
print 'sPvvCv', len(PAXB_VTP_sig_genes.intersection(CERC_VTP_sig_genes))
print 'sall', len(CERC_VTP_sig_genes.intersection(PAXB_BS_sig_genes).intersection(PAXB_VTP_sig_genes))


shared_VTP_ids = PAXB_VTP_genes.intersection(CERC_VTP_genes)
shared_sig_VTP_ids = PAXB_VTP_sig_genes.intersection(CERC_VTP_sig_genes)

shared_all_ids = CERC_VTP_genes.intersection(PAXB_BS_genes).intersection(PAXB_VTP_genes)
shared_all_sig_ids = CERC_VTP_sig_genes.intersection(PAXB_BS_sig_genes).intersection(PAXB_VTP_sig_genes)


cor_2dicts(shared_VTP_ids,PAXB_VTP_M_F_dict,CERC_VTP_M_F_dict)
cor_2dicts(shared_sig_VTP_ids,PAXB_VTP_M_F_dict,CERC_VTP_M_F_dict)

cor_2dicts(shared_all_ids,PAXB_VTP_M_F_dict,CERC_VTP_M_F_dict)
cor_2dicts(shared_all_sig_ids,PAXB_VTP_M_F_dict,CERC_VTP_M_F_dict)

cor_2dicts(shared_all_ids,PAXB_BS_M_F_dict,CERC_VTP_M_F_dict)
cor_2dicts(shared_all_sig_ids,PAXB_BS_M_F_dict,CERC_VTP_M_F_dict)

w = open('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in PAXB_BS_genes:
    w.write(names[ensid] + '\n')

w.close()
w = open('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_sig_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in PAXB_BS_sig_genes:
    w.write(names[ensid] + '\n')

w.close()


w = open('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in PAXB_VTP_genes:
    w.write(names[ensid] + '\n')

w.close()
w = open('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_sig_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in PAXB_VTP_sig_genes:
    w.write(names[ensid] + '\n')

w.close()

w = open('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in CERC_VTP_genes:
    w.write(names[ensid] + '\n')

w.close()
w = open('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_sig_ASE.l','w')
names = module.ENSid_Human_Dict()
for ensid in CERC_VTP_sig_genes:
    w.write(names[ensid] + '\n')

w.close()




nc_cor,total_cor = neighbor_correlation(PAXB_BS_M_F_dict,start_dict)
print nc_cor.mean(),total_cor.mean()
plt.hist(nc_cor,color='g',bins=30)
plt.hist(total_cor,color='b',alpha=.2,bins=30)
plt.show()
plt.violinplot([nc_cor,total_cor])
plt.show()
'''
'''
log_ase_dict = logfc_ase_dict(ase_dict)
print log_ase_dict.keys()
start_dict = load_starts()


w1 = open('/home/james/Dropbox/Miller/data/BED/ensGene_collapsed_TSS.bed','w')
w2 = open('/home/james/Dropbox/Miller/data/BED/ensGene_collapsed_TSS_names.bed','w')
for ensid in start_dict.keys():
    chrome,start =  start_dict[ensid]
    w1.write(chrome + '\t' + str(start) + '\t' + str(int(start)+1) + '\t' + ensid + '\n')
    w2.write(chrome + '\t' + str(start) + '\t' + str(int(start)+1) + '\t' + ensNames[ensid] + '\n')

w1.flush()
w2.flush()
w1.close()
w2.close()
'''

##there doesn't seem to be good evidence for co-regulated TADs. Some thought that the data might be clumpier than random, implying co-regulation between neighbors?
##plot_ase(log_ase_dict,start_dict,'chrXX')
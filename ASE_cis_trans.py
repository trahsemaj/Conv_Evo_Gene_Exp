__author__ = 'james'
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
import module
import math
import random
from scipy.stats import gaussian_kde
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
import subprocess as sp
import matplotlib
##import seaborn as sns
##sns.set(color_codes=True)
##sns.set_style("ticks")
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

##takes a key and a list of ase dicts
##plots the
def plot_key_ase_dict(key,dict_list,LABELS='',COLORS=''):

    plot_vals = []
    plt.figure()
    for ind,ase_dict in enumerate(dict_list):
        cur_counts = np.array(ase_dict[key])
        cur_ase_vals = cur_counts[:,0]/cur_counts[:,1]
        cur_ase_vals = np.log2(cur_ase_vals)
        cur_ase_vals = cur_ase_vals[~np.isnan(cur_ase_vals)]
        plot_vals.append(cur_ase_vals)
        print cur_counts
        xs = np.zeros(len(cur_ase_vals)) + ind
        plt.scatter(xs,cur_ase_vals,marker='.',color='black')
        plt.errorbar(ind,cur_ase_vals.mean(),yerr=stats.sem(cur_ase_vals),fmt='o',ecolor=COLORS[ind])

    plot_vals = np.array(plot_vals)


    xs = np.arange(len(dict_list))
    ##full_xs = np.resize(xs,(len(xs),5))

    ##plt.scatter(full_xs,plot_vals,marker='.')
    ##plt.errorbar(xs,plot_vals.mean(axis=1),yerr=stats.sem(plot_vals,axis=1),fmt='_')


    zero_xs = np.arange(len(dict_list)+1)-.5
    zeros = np.zeros(len(zero_xs))
    plt.plot(zero_xs,zeros,color='grey',ls='dashed')
    plt.xticks(xs,LABELS,rotation=45)
    plt.xlim([-.5,len(dict_list)-.5])
    plt.ylabel('log2(F/M)')


##takes a key, an exp_dict, and a series of indexes
##plots the fc from ref_index mean
def plot_key_expDiff_dict(key,exp_dict,test_indexes,ref_index,LABELS='',COLORS=''):


    plt.figure()
    for ind,cur_pop in enumerate(test_indexes):
        cur_exp = np.array(exp_dict[key][cur_pop])
        cur_ref_mean = np.array(exp_dict[key][ref_index]).mean()
        cur_exp_diff = cur_exp/cur_ref_mean
        cur_exp_diff = np.log2(cur_exp_diff)
        xs = np.zeros(len(cur_exp_diff)) + ind
        plt.scatter(xs,cur_exp_diff,marker='.',color='black')
        plt.errorbar(ind,cur_exp_diff.mean(),yerr=stats.sem(cur_exp_diff),fmt='o',ecolor=COLORS[ind])



    xs = np.arange(len(test_indexes))
    ##full_xs = np.resize(xs,(len(xs),5))

    ##plt.scatter(full_xs,plot_vals,marker='.')
    ##plt.errorbar(xs,plot_vals.mean(axis=1),yerr=stats.sem(plot_vals,axis=1),fmt='_')


    zero_xs = np.arange(len(test_indexes)+1)-.5
    zeros = np.zeros(len(zero_xs))
    plt.plot(zero_xs,zeros,color='grey',ls='dashed')
    plt.xticks(xs,LABELS,rotation=45)
    plt.xlim([-.5,len(test_indexes)-.5])
    plt.ylabel('log2(F/M)')




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



##returns a dict of the mean exp of each gene
def load_mean_exp(exp_file):

    f = open(exp_file)
    f.readline()

    toret = {}

    for line in f:
        line = line.strip()
        split = line.split('\t')

        cur_id = split[0]

        cur_exp = np.array(map(lambda x: float(x), split[1:]))
        toret[cur_id] = cur_exp.mean()

    return toret

##takes x and y vals
##plots a red dashed PCA1 line fo the given data
def plotPCA(xvals,yvals,COLOR='r'):
    xvals = np.array(xvals)
    yvals = np.array(yvals)
    xData = np.reshape(xvals, (len(xvals), 1))
    yData = np.reshape(yvals, (len(yvals), 1))

    data = np.hstack((xData,yData))

    mean = data.mean(axis=0)
    data = data-mean

    eigenVectors, eigenValues, V = np.linalg.svd(
    data.T, full_matrices=False)
    ##print eigenVectors, eigenValues, V

    projected = np.dot(data,eigenVectors)
    sigma = projected.std(axis=0).mean()
    ##print(eigenVectors)
    ##print eigenValues
    ##print data[1]


    xs = np.arange(-10,10,.1)
    ##print eigenVectors
    plt.plot(xs,eigenVectors[1,0]/eigenVectors[0,0]*xs+yvals.mean(),ls='dashed',color=COLOR)
    m,b = np.polyfit(xvals,yvals,1)

    ##plt.plot(xs,m*xs+b,color='green',linestyle='dashed')


##plots correlation between 2 dicts values
def plot_cor(dict1,dict2,DENSITY=False,COLOR='b',XYLIM=[-5,5],UUDD=True,SYM='.'):
    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    vals1 = []
    vals2 = []
    up_up = 0
    down_down=0
    up_down=0
    down_up = 0
    for key in shared_keys:
        vals1.append(dict1[key])
        vals2.append(dict2[key])

        cur_1 = dict1[key]
        cur_2 = dict2[key]

        if cur_1 > 0:
            if cur_2 > 0:
                up_up += 1
            else:
                up_down += 1
        else:
            if cur_2 > 0:
                down_up += 1
            else:
                down_down += 1


        if math.fabs(dict1[key]) > 3:
            pass
            ##print key, dict1[key]
    if UUDD:
        print 'up_up,up_down,down_up,down_down'
        print up_up,up_down,down_up,down_down
    ##print len(shared_keys)
    r,pval = stats.spearmanr(vals1,vals2)
    xy = np.vstack([vals1,vals2])
    psize = 5
    if DENSITY:
        z = gaussian_kde(xy)(xy)
        plt.scatter(vals1,vals2, c=z, s=psize, edgecolor='')
    else:
        plt.scatter(vals1,vals2, color=COLOR,marker=SYM)
    plotPCA(vals1,vals2,COLOR=COLOR)
    r,pval =  stats.spearmanr(vals1,vals2)
    plt.ylim(XYLIM)
    plt.xlim(XYLIM)
    ##plt.tight_layout()

    ##print r,pval
    plt.xlabel(r)
    return r,pval
    ##plt.scatter(vals1,vals2)
    ##plt.show()

##takes a file with ASE data, and a file with
##diff than in ASE_plots
def plot_cis_trans(ASE_list_file,trans_exp_diff_file,sig_ase_file,sig_trans_file):


    plt.figure(figsize=(12,12))
    ##in the form log2(PAXB/(RABS))
    cis_dict = load_logfc(ASE_list_file)
    ##using the PAXB_RABS file, in the form log2(PAXB/RABS)
    trans_dict = load_logfc(trans_exp_diff_file)

    sig_cis_dict = load_logfc(sig_ase_file)
    sig_trans_dict = load_logfc(sig_trans_file)


    ##find the common genes twixt the two

    genes = set(cis_dict.keys()).intersection(trans_dict.keys())

    sig_cis_genes = sig_cis_dict.keys()
    sig_trans_genes = sig_trans_dict.keys()

    sig_cis_trans = set(sig_cis_genes).intersection(sig_trans_genes)


    cis_diff = []
    trans_diff = []

    I = 0
    II = 0
    III=0
    IV = 0
    for gene in genes:

        ##transform ratio to log2
        ase_val = cis_dict[gene]
        ##ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
        trans_val = trans_dict[gene]

        if ase_val > 0:
            if trans_val>0:
                I+=1
            if trans_val<0:
                II+=1
        if ase_val < 0:
            if trans_val>0:
                IV+=1
            if trans_val<0:
                III+=1


        cis_diff.append(ase_val)
        trans_diff.append(trans_val)

    plt.scatter(trans_diff,cis_diff,marker='.',alpha=.3,color='black',label='N.S.')


    cis = []
    trans = []
    for gene in sig_cis_genes:

        if gene not in genes:
            continue

        ase_val = cis_dict[gene]
        ##ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
        trans_val = trans_dict[gene]

        cis.append(ase_val)
        trans.append(trans_val)

    plt.scatter(trans,cis,marker='.',alpha=.3,color='red',label='FDR < .05 displaying ASE')

    cis = []
    trans = []
    for gene in sig_trans_genes:

        if gene not in genes:
            continue

        ase_val = cis_dict[gene]
        ##ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
        trans_val = trans_dict[gene]

        cis.append(ase_val)
        trans.append(trans_val)

    plt.scatter(trans,cis,marker='.',alpha=.3,color='blue',label='Significant between Parentals')

    cis = []
    trans = []
    for gene in sig_cis_trans:

        if gene not in genes:
            continue

        ase_val = cis_dict[gene]
        ##ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
        trans_val = trans_dict[gene]

        cis.append(ase_val)
        trans.append(trans_val)

    plt.scatter(trans,cis,marker='.',alpha=.3,color='green',label='Significant in both')

    xs = np.arange(-5,5,.01)
    zeros = np.zeros(len(xs))

    ##k0 = smooth.NonParamRegression(trans_diff,cis_diff, method=npr_methods.SpatialAverage())
    ##k0.fit()
    m,b = np.polyfit(trans_diff,cis_diff,1)

    plt.plot(xs,m*xs+b,color='red',linestyle='dashed')
    ##plotPCA(cis_diff,trans_diff)
    plt.plot(xs,xs,color='black')
    plt.plot(xs,zeros,color='black')
    plt.plot(zeros,xs,color='black')
    plt.xlim([-5,5])
    plt.ylim([-5,5])
    plt.xlabel('log2( PAXB(FW)_p / RABS(M)_p )\n"trans"',fontsize=26)
    plt.ylabel('log2( PAXB(FW)_F1 / RABS(M)_F1 )\n"cis"',fontsize=26)
    plt.title('PAXB(FW) Changes in Gene Expression',fontsize=40)
    print I,II,III,IV
    plt.legend()

    ##plt.show()



##takes a file with ASE data, and a file with
##diff than in ASE_plots
def plot_cis_trans_v2(f1_dict,expDiff_dict,sig_cis_dict,sig_trans_dict):

    plt.figure(figsize=(12,12))

    shared_ids = set(f1_dict.keys()).intersection(set(expDiff_dict.keys()))



    sig_cis_ids = set(sig_cis_dict.keys()).intersection(shared_ids)
    sig_trans_ids = set(sig_trans_dict.keys()).intersection(shared_ids)


    sig_cis_trans = sig_cis_ids.intersection(sig_trans_ids)
    no_sig_change = shared_ids.difference(sig_cis_ids.union(sig_trans_ids))

    sig_cis_only = sig_cis_ids.difference(sig_cis_trans)
    sig_trans_only = sig_trans_ids.difference(sig_cis_trans)

    concordant_cis_trans = set()
    discordant_cis_trans = set()


    for ensid in sig_cis_trans:
        cur_expDiff = expDiff_dict[ensid]
        cur_ase = f1_dict[ensid]

        if cur_expDiff > 0:
            if cur_ase > 0:
                concordant_cis_trans.add(ensid)
            else:
                discordant_cis_trans.add(ensid)
        else:
            if cur_ase < 0:
                concordant_cis_trans.add(ensid)
            else:
                discordant_cis_trans.add(ensid)

    print 'total',len(shared_ids)
    print 'no_sig',len(no_sig_change)
    print 'cis only', len(sig_cis_only)
    print 'trans only', len(sig_trans_only)
    print 'cis_trans',len(sig_cis_trans)
    print 'discordant', len(discordant_cis_trans)
    print 'concordant', len(concordant_cis_trans)

    total_parental,total_hybrid = select_from_dicts(shared_ids,expDiff_dict,f1_dict)

    no_sig_change_parental,no_sig_change_hybrid = select_from_dicts(no_sig_change,expDiff_dict,f1_dict)
    sig_cis_only_parental,sig_cis_only_hybrid = select_from_dicts(sig_cis_only,expDiff_dict,f1_dict)
    sig_trans_only_parental,sig_trans_only_hybrid = select_from_dicts(sig_trans_only,expDiff_dict,f1_dict)
    concordant_cis_trans_parental,concordant_cis_trans_hybrid = select_from_dicts(concordant_cis_trans,expDiff_dict,f1_dict)
    discordant_cis_trans_parental,discordant_cis_trans_hybrid = select_from_dicts(discordant_cis_trans,expDiff_dict,f1_dict)

    plt.scatter(no_sig_change_hybrid,no_sig_change_parental,marker='.',color='grey',label='No siginificant change')
    plt.scatter(sig_cis_only_hybrid,sig_cis_only_parental,marker='.',color='blue',label='Cis change')
    plt.scatter(sig_trans_only_hybrid,sig_trans_only_parental,marker='.',color='red',label='Trans change')
    plt.scatter(concordant_cis_trans_hybrid,concordant_cis_trans_parental,marker='.',color='purple',label='Concordant cis + trans')
    plt.scatter(discordant_cis_trans_hybrid,discordant_cis_trans_parental,marker='.',color='green',label='Discordant cis + trans')


    plotPCA(total_hybrid,total_parental,COLOR='black')
    xs = np.arange(-5,5,.01)
    zeros = np.zeros(len(xs))
    plt.plot(xs,xs,color='black')
    plt.plot(xs,zeros,color='black')
    plt.plot(zeros,xs,color='black')
    plt.legend()

##takes a set, and two dicts
##selects keys in the set from each dict
##returns two lists, the values of each element in select set in each dict
def select_from_dicts(select_set,dict1,dict2):

    list1 = []
    list2 = []

    for element in select_set:
        list1.append(dict1[element])
        list2.append(dict2[element])

    return list1,list2




def get_sig_trans(exp_dict,fresh_index,marine_index,cis_dict,FDR=.01):

    shared_ids = set(exp_dict.keys()).intersection(set(cis_dict.keys()))

    trans_pvals = []
    trans_ensids = []
    trans_dict = {}
    sig_trans_dict = {}

    for ensid in shared_ids:

        f1_ratio = cis_dict[ensid]

        cur_reads = exp_dict[ensid]
        fresh_reads = cur_reads[fresh_index]
        marine_reads = cur_reads[marine_index]

        fresh_total = float(sum(fresh_reads)) + 1e-20
        marine_total = float(sum(marine_reads)) + 1e-20
        total_reads = fresh_total + marine_total

        fresh_expect = (2. ** f1_ratio)/((2. ** f1_ratio) + 1) *  total_reads
        marine_expect = total_reads - fresh_expect



        cis_trans = math.log(float(fresh_total)/total_reads,2)

        trans = cis_trans - f1_ratio
        trans_dict[ensid] = trans
        G_val = 2 * ((fresh_total * math.log(fresh_total/fresh_expect)) + (marine_total * math.log(marine_total/marine_expect)))
        pval = stats.chisqprob(G_val,1)
        trans_pvals.append(pval)
        trans_ensids.append(ensid)

    trans_pvals = np.array(trans_pvals)
    trans_ensids = np.array(trans_ensids)

    trans_ensids = trans_ensids[np.argsort(trans_pvals)]
    trans_pvals = trans_pvals[np.argsort(trans_pvals)]
    cutoff = 0.
    for k in range(len(trans_pvals)):
        cutoff += 1.
        cur_fdr = (cutoff/len(trans_pvals)) * FDR
        ##print cutoff,len(pval_array)
        ##print cur_fdr
        ##print pval_array[k]
        if cur_fdr < trans_pvals[k]:
            print cutoff
            break

    trans_ensids = list(trans_ensids)
    for ensid in trans_ensids[:int(cutoff)]:
        sig_trans_dict[ensid] = trans_dict[ensid]

    return sig_trans_dict


##takes two dicts - parental fc and hybrid AIB
##returns the calculated trans values
def calc_trans(parental,hybrid,outfile='',total=True):


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

        if not total:
            percent_cis = cis/cis_trans


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






def plot_percent_cis(dict1,dict2,name1,name2):
    dict1_percent_cis = []
    dict2_percent_cis = []


    shared_keys = set(dict1.keys()).intersection(dict2.keys())
    print len(shared_keys)

    for key in shared_keys:
        dict1_percent_cis.append(dict1[key])
        dict2_percent_cis.append(dict2[key])

    dict1_percent_cis = np.array(dict1_percent_cis)
    dict2_percent_cis = np.array(dict2_percent_cis)
    print dict1_percent_cis
    plotSmoothedHist(dict1_percent_cis,name1,x_range=([-2,2]))
    plotSmoothedHist(dict2_percent_cis,name2,x_range=([-2,2]))
    plt.legend()
    ##plt.show()

def hist_dict_list(inlist,namelist):

    for cur_dict,cur_name in zip(inlist,namelist):

        cur_vals = []

        for key in cur_dict.keys():
            cur_vals.append(cur_dict[key])

        cur_vals = np.array(cur_vals)
        ##print cur_vals[np.isnan(cur_vals)]
        ##print cur_vals[np.isinf(cur_vals)]
        ##cur_vals = cur_vals[~np.isnan(cur_vals)]
        ##cur_vals = cur_vals[~np.isinf(cur_vals)]
        ##cur_vals = np.reshape(cur_vals,(len(cur_vals),1))

        print stats.normaltest(cur_vals)
        ##print stats.shapiro(cur_vals[:,None])
        plotSmoothedHist(cur_vals,cur_name,x_range=[-5,5])


    ##plt.show()

##works with any tsv exp file (assumes 3 reps per condition)
def load_exp_table_reps(exp_file,CUTOFF=.1,REPS=3):
    f = open(exp_file, 'r')

    f.readline()

    exp_dict = {}

    for line in f:
        line = line.strip()
        split = line.split('\t')
        name = split[0]
        cur_exp = map(lambda x: float(x), split[1:])
        ##filter out non-expressed genes
        if cur_exp.count(0) > 4:
            continue

        if sum(cur_exp)/len(cur_exp) < CUTOFF:
            continue

        exp_reps = []

        ##conditions = len(cur_exp)/3
        conditions = range(len(cur_exp)/REPS)
        for i in conditions:
            cur_condition = cur_exp[i*REPS:(i+1)*REPS]
            cur_condition = cur_exp[i*REPS:(i+1)*REPS]
            exp_reps.append(cur_condition)

        exp_dict[name] = exp_reps

    return exp_dict

##returns a set of all genes in bitecode
def loadBiteCode(fpath='/home/james/Dropbox/Miller/data/BED/BITECODE.l'):
    toret = set()

    f = open(fpath, 'r')

    for line in f:
        line = line.strip()
        gene = line.upper()
        toret.add(gene)

    return toret

##takes a gene expression dictionary (in name form), a list of gene sets (sets of gene names)
##the empty list will test every gene in the genome
def plotCoV(exp_dict,gset_list,totest,ref_index):

    ##will be a list of arrays of CoV values for each gene set
    total_CoV = []

    for gset in gset_list:

        gset_CoV = []
        mean_exp1 = []
        mean_exp2 = []
        mean_ref = []

        for name in exp_dict.keys():

            if (name.upper().split(' ')[0] in gset) or not gset:
                expression = exp_dict[name]
                exp1 = list(expression[totest[0]])
                exp2 = list(expression[totest[1]])
                ref_exp = np.array(expression[ref_index])
                if (sum(exp1) < .1) or (sum(exp2) < .1) or (sum(list(ref_exp)) < .1):
                    continue
                combined = exp1 + exp2
                combined = [sum(exp1)/len(exp1),sum(exp2)/len(exp2)]
                combined = np.array(combined)
                combined = combined/ref_exp.mean()

                cv = combined.std()/combined.mean()
                cor_cv = cv * (13./12)
                gset_CoV.append(cor_cv)

                mean_exp1.append(sum(exp1)/len(exp1))
                mean_exp2.append(sum(exp2)/len(exp2))
                mean_ref.append(ref_exp.mean())


        mean_exp1 = np.array(mean_exp1)
        mean_exp2 = np.array(mean_exp2)
        mean_ref = np.array(mean_ref)

        mean_exp1 = np.log2(mean_exp1)
        mean_exp2 = np.log2(mean_exp2)
        mean_ref = np.log2(mean_ref)

        print stats.pearsonr(mean_exp1,mean_exp2)
        print stats.pearsonr(mean_exp1,mean_ref)
        print stats.pearsonr(mean_exp2,mean_ref)
        print
        print stats.spearmanr(mean_exp1,mean_exp2)
        print stats.spearmanr(mean_exp1,mean_ref)
        print stats.spearmanr(mean_exp2,mean_ref)
        print
        print stats.spearmanr(mean_exp1-mean_ref,mean_exp2-mean_ref)
        print stats.spearmanr(mean_ref-mean_exp1,mean_exp2-mean_exp1)
        print stats.spearmanr(mean_exp1-mean_exp2,mean_ref-mean_exp2)
        print
        print stats.pearsonr(mean_exp1-mean_ref,mean_exp2-mean_ref)
        print stats.pearsonr(mean_ref-mean_exp1,mean_exp2-mean_exp1)
        print stats.pearsonr(mean_exp1-mean_exp2,mean_ref-mean_exp2)
        print

        gset_CoV = np.array(gset_CoV)
        print np.median(gset_CoV)
        print gset_CoV.std()
        total_CoV.append(gset_CoV)
    w,pval= stats.mannwhitneyu(total_CoV[0],total_CoV[1])
    print pval
    print np.median(total_CoV[0]),np.median(total_CoV[1])
    plt.boxplot(total_CoV)
    plt.xlabel('p = %5f' % pval)

##takes a dictionary, and returns keys that are in filter (or not in filter in infilter=False)
def filter_dict(indict,filter,infilter=True):

    toret = {}

    for key in indict.keys():
        if (key in filter) * infilter:
            toret[key] = indict[key]

    return toret

##plots a histogram of ASE values
def plot_ASE_vals(ase_dict_list,label_list):

    plt.figure()

    for ase_dict,cur_label in zip(ase_dict_list,label_list):

        cur_ase_vals = np.array(ase_dict.values())
        print cur_ase_vals.mean()
        print np.median(cur_ase_vals)
        print stats.ttest_1samp(cur_ase_vals,0)
        print stats.wilcoxon(cur_ase_vals)
        plotSmoothedHist(cur_ase_vals,cur_label,x_range=[-3,3])
    plt.legend()

##calculates the mean exp in a given condition
##retuns of a dict, fpkm cutoff=.1
def mean_exp_cond_cor(gene_dict,cond,REPS=5,LOG=10):
    toret = {}
    for gene in gene_dict.keys():
        cur_exp = gene_dict[gene]
        cur_cond_exp = cur_exp[cond]
        cur_cond_exp_mean = sum(cur_cond_exp)/len(cur_cond_exp)
        if cur_cond_exp_mean < .1:
            continue
        if LOG:
            cur_cond_exp_mean = math.log(cur_cond_exp_mean,LOG)
        toret[gene] = cur_cond_exp_mean
    return toret

##takes 4 dicts, of mean expression and ASE from tissues A and B
##plots the ABS(ase) for each
def tissue_ase_exp(meanA,meanB,aseA,aseB,DENSITY=False,FIG_BASE=''):

    shared_keys = set(meanA.keys()).intersection(meanB.keys()).intersection(aseA.keys()).intersection(aseB.keys())
    print len(shared_keys)
    percent_change_vals = []
    ase_A_vals = []
    ase_B_vals = []
    ase_diff_vals = []
    for gene in shared_keys:
        cur_a_exp = meanA[gene]
        cur_b_exp = meanB[gene]
        cur_a_ase = aseA[gene]
        cur_b_ase = aseB[gene]
        cur_percent_change = math.fabs(cur_a_exp-cur_b_exp)/(cur_a_exp+cur_b_exp)
        cur_percent_change = math.fabs( math.log(cur_a_exp,2) - math.log(cur_b_exp,2))
        cur_ase_diff = math.fabs(cur_a_ase-cur_b_ase)
        percent_change_vals.append(cur_percent_change)
        ase_A_vals.append(cur_a_ase)
        ase_B_vals.append(cur_b_ase)
        ase_diff_vals.append(cur_ase_diff)

    percent_change_vals = np.array(percent_change_vals)
    ase_A_vals = np.array(ase_A_vals)
    ase_B_vals = np.array(ase_B_vals)
    ase_diff_vals = np.array(ase_diff_vals)
    if DENSITY:
        f, ax = plt.subplots(figsize=(10,10))
        sns.kdeplot(ase_A_vals, ase_B_vals, ax=ax,alpha=.7,cmap=sns.dark_palette("red",as_cmap=True),n_levels=7, shade=False)
        #sns.rugplot(ase_A_vals, color="g", ax=ax)
        #sns.rugplot(ase_B_vals, vertical=True, ax=ax)
        ax.scatter(ase_A_vals,ase_B_vals,c='k',marker='.')
        plotPCA(ase_A_vals,ase_B_vals,COLOR='blue')
        ax.set_xlim([-4,4])
        ax.set_xlabel('ASE_1')
        ax.set_ylim([-4,4])
        ax.set_ylabel('ASE_2')
        if FIG_BASE:
            plt.savefig(FIG_BASE + 'shared_ase.pdf')

        sns.jointplot(percent_change_vals,ase_diff_vals,kind='reg',marker='.',color='m')
        plt.xlim([0,np.max(percent_change_vals)+.1])
        plt.ylim([0,np.max(ase_diff_vals)+.1])
        plt.xlabel('Mean Expression Difference (log2)')
        plt.ylabel('Mean Cis Difference (log2)')
        plt.tight_layout()
        if FIG_BASE:
            plt.savefig(FIG_BASE + 'tissue_ase_exp.pdf')

        ##ax.scatter(percent_change_vals,ase_diff_vals,color='r')
    else:
        plt.scatter(percent_change_vals,ase_diff_vals,color='r')
    print stats.linregress(percent_change_vals,ase_diff_vals)
    print stats.linregress(ase_A_vals,ase_B_vals)
    m,b,r,pval,stderr =  stats.linregress(percent_change_vals,ase_diff_vals)
    xs = np.arange(0,int(np.max(percent_change_vals)),.05)
    ys = xs*m + b
    ##ax.plot(xs,ys)
    ##plt.colorbar()
    plt.show()

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

##SupFig CoV smaller for BC genes


'''
parental_exp = load_exp_table_reps('/home/james/Dropbox/Miller/data/RNA_seq_exp/CRPL_norm_2kmask_names.fpkm_table')
biteCode_set = loadBiteCode()
rand_gset = []
for key in random.sample(parental_exp.keys(),1000):
    if key.startswith('ENS'):
        continue
    new_key = key.upper().split(' ')[0]
    rand_gset.append(new_key)
plotCoV(parental_exp,[rand_gset,biteCode_set],[0,2],1)
plt.savefig(fig_base + 'PAXB_CERC_BC_CoV.pdf')

fpath_base = '/home/james/Dropbox/Miller/data/ASE/'

PAXB_VTP_ase = load_logfc(fpath_base + 'PAXB_VTP_logAI.l')
CERC_VTP_ase = load_logfc(fpath_base + 'CERC_VTP_logAI.l')
plot_ASE_vals([PAXB_VTP_ase,CERC_VTP_ase],['PAXB','CERC'])
plt.savefig(fig_base + 'PAXB_CERC_ASE_hist.pdf')
plt.show()
'''
'''
if biteCode:
        in_biteCode = []
        bc_genes = loadBiteCode()
        for ensID in ENSIDS:
            name = ensNames[ensID].split(' ')[0]
            if ensID in hnames:
                name = hnames[ensID]
            in_biteCode.append((name.upper() in bc_genes) or (ensNames[ensID].split(' ')[0].upper() in bc_genes))
        in_biteCode = np.array(in_biteCode)
'''

fpath_base = '/home/james/Dropbox/Miller/data/ASE/'
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE/'
parental_exp = load_exp_table_reps('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table')

PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_BS_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_AI_all.tsv')
CERC_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_AI_all.tsv')


fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Old_6/'
'''
plot_key_ase_dict('ENSGACT00000013425',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'snai1a_CERC_PAXB_VTP_ase.pdf')
plot_key_ase_dict('ENSGACT00000024635',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'pdgfrl_CERC_PAXB_VTP_ase.pdf')
plot_key_ase_dict('ENSGACT00000012993',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'znf106_CERC_PAXB_VTP_ase.pdf')
plot_key_ase_dict('ENSGACT00000004902',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'plod2_CERC_PAXB_VTP_ase.pdf')
plot_key_ase_dict('ENSGACT00000014468',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'bmp1a_CERC_PAXB_VTP_ase.pdf')
plot_key_ase_dict('ENSGACT00000008011',[PAXB_VTP_ase_dict,CERC_VTP_ase_dict],LABELS=['PAXB','CERC'],COLORS=['b','g'])
plt.savefig(fig_base + 'ctnnb1_CERC_PAXB_VTP_ase.pdf')
##snai1a - ENSGACT00000013425
##pdgfrl - ENSGACT00000024635
##ZNF106 - ENSGACT00000012993
##Plod2 - ENSGACT00000004902
##bmp1a - ENSGACT00000014468
##ctnnb1 - ENSGACT00000008011

plot_key_expDiff_dict('ENSGACT00000013425',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'snai1a_CERC_PAXB_VTP_expDiff.pdf')
plot_key_expDiff_dict('ENSGACT00000024635',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'pdgfrl_CERC_PAXB_VTP_expDiff.pdf')
plot_key_expDiff_dict('ENSGACT00000012993',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'znf106_CERC_PAXB_VTP_expDiff.pdf')
plot_key_expDiff_dict('ENSGACT00000004902',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'Plod2_CERC_PAXB_VTP_expDiff.pdf')
plot_key_expDiff_dict('ENSGACT00000014461',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'bmp1a_CERC_PAXB_VTP_expDiff.pdf')
plot_key_expDiff_dict('ENSGACT00000008011',parental_exp,[2,0],1,COLORS=['b','g'],LABELS=['PAXB','CERC'])
plt.savefig(fig_base + 'ctnnb1_CERC_PAXB_VTP_expDiff.pdf')
'''
##snai1a - ENSGACT00000013425
##pdgfrl - ENSGACT00000024635
##ZNF106 - ENSGACT00000012993
##Plod2 - ENSGACT00000004902
##bmp1a - ENSGACT00000014461
##ctnnb1 - ENSGACT00000008011
plt.show()

fpath_base = '/home/james/Dropbox/Miller/data/ASE/'
parental_exp = load_exp_table_reps('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table')

PAXB_VTP_expdiff = load_logfc(fpath_base + 'PAXB_RABS_logfc.l')
PAXB_VTP_expdiff_sig = load_logfc(fpath_base + 'PAXB_RABS_sig_logfc.l')
PAXB_VTP_ase = load_logfc(fpath_base + 'PAXB_VTP_logAI.l')
PAXB_VTP_sig_ase = load_logfc(fpath_base + 'PAXB_VTP_logAI_sig.l')
PAXB_VTP_ase_sigExpdiff = filter_dict(PAXB_VTP_ase,PAXB_VTP_expdiff_sig)
mean_exp = load_mean_exp('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table')

CERC_VTP_expdiff = load_logfc(fpath_base + 'CERC_RABS_logfc.l')
CERC_VTP_expdiff_sig = load_logfc(fpath_base + 'CERC_RABS_sig_logfc.l')
CERC_VTP_ase = load_logfc(fpath_base + 'CERC_VTP_logAI.l')
CERC_VTP_sig_ase = load_logfc(fpath_base + 'CERC_VTP_logAI_sig.l')
CERC_VTP_ase_sigExpdiff = filter_dict(CERC_VTP_ase,CERC_VTP_expdiff_sig)

PAXB_BS_ase = load_logfc(fpath_base + 'PAXB_BS_logAI.l')
PAXB_BS_sig_ase = load_logfc(fpath_base + 'PAXB_BS_logAI_sig.l')

PAXB_BS_Dref_ase = load_logfc(fpath_base + 'PAXB_BS_Dref_logAI.l')
PAXB_BS_Dref_sig_ase = load_logfc(fpath_base + 'PAXB_BS_Dref_logAI_sig.l')

FTC_BS_ase = load_logfc(fpath_base + 'FTC_BS_logAI.l')
FTC_BS_sig_ase = load_logfc(fpath_base + 'FTC_BS_logAI_sig.l')

PAXB_CERC_expdiff = load_logfc(fpath_base + 'PAXB_CERC_logfc.l')
PAXB_CERC_expdiff_sig = load_logfc(fpath_base + 'PAXB_CERC_sig_logfc.l')


PAXB_CERC_bc_expdiff = load_logfc(fpath_base + 'PAXB_CERC_bc_logfc.l')
PAXB_RABS_bc_expdiff = load_logfc(fpath_base + 'PAXB_RABS_bc_logfc.l')
CERC_RABS_bc_expdiff = load_logfc(fpath_base + 'CERC_RABS_bc_logfc.l')

RABS_VTP_expdiff = load_logfc(fpath_base + 'PAXB_CERC_bc_logfc.l')




print 'CERC trans all'
CERC_trans,CERC_percent_cis = calc_trans(CERC_VTP_expdiff,CERC_VTP_ase)
print 'CERC trans - ExpDiff sig'
sig_CERC_trans,sig_CERC_percent_cis = calc_trans(CERC_VTP_expdiff_sig,CERC_VTP_ase)
sig_CERC_trans_net,sig_CERC_percent_cis_net = calc_trans(CERC_VTP_expdiff_sig,CERC_VTP_ase,total=False)

print 'PAXB trans all'
PAXB_trans,PAXB_percent_cis = calc_trans(PAXB_VTP_expdiff,PAXB_VTP_ase)
print 'PAXB trans - ExpDiff sig'
sig_PAXB_trans,sig_PAXB_percent_cis = calc_trans(PAXB_VTP_expdiff_sig,PAXB_VTP_ase)
sig_PAXB_trans_net,sig_PAXB_percent_cis_net = calc_trans(PAXB_VTP_expdiff_sig,PAXB_VTP_ase,total=False)


PAXB_VTP_expdiff_filtered = filter_dict(PAXB_VTP_expdiff,PAXB_VTP_ase)
CERC_VTP_expdiff_filtered = filter_dict(CERC_VTP_expdiff,CERC_VTP_ase)

PAXB_VTP_expdiff_sig_filtered = filter_dict(PAXB_VTP_expdiff_sig,PAXB_VTP_ase)
CERC_VTP_expdiff_sig_filtered = filter_dict(CERC_VTP_expdiff_sig,CERC_VTP_ase)



CERC_trans_sigOnly = get_sig_trans(parental_exp,0,1,CERC_VTP_ase)
PAXB_trans_sigOnly = get_sig_trans(parental_exp,2,1,PAXB_VTP_ase)

PAXB_VTP_ase_vals = np.fabs(np.array(PAXB_VTP_ase.values()))
CERC_VTP_ase_vals = np.fabs(np.array(CERC_VTP_ase.values()))
PAXB_VTP_sig_ase_vals = np.fabs(np.array(PAXB_VTP_sig_ase.values()))
CERC_VTP_sig_ase_vals = np.fabs(np.array(CERC_VTP_sig_ase.values()))


print 'ase'
print PAXB_VTP_ase_vals.mean()
print np.median(PAXB_VTP_ase_vals)
print CERC_VTP_ase_vals.mean()
print np.median(CERC_VTP_ase_vals)
print PAXB_VTP_sig_ase_vals.mean()
print np.median(PAXB_VTP_sig_ase_vals)
print CERC_VTP_sig_ase_vals.mean()
print np.median(CERC_VTP_sig_ase_vals)

print 'CxP diff',stats.mannwhitneyu(PAXB_VTP_ase_vals,CERC_VTP_ase_vals)
print stats.mannwhitneyu(PAXB_VTP_sig_ase_vals,CERC_VTP_sig_ase_vals)

print plot_cor(PAXB_VTP_expdiff,CERC_VTP_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='blue')
print plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='red')

bc_keys = PAXB_RABS_bc_expdiff.keys()
bc_length = len(bc_keys)


fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/SupFigs/'

##suppFigs - BC cor
'''
plt.figure()
print 'CERC'
print plot_cor(CERC_VTP_expdiff,PAXB_CERC_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='blue')
print plot_cor(CERC_RABS_bc_expdiff,PAXB_CERC_bc_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='red')
plt.savefig(fig_base + 'CERC_RABS_PAXB_CERC_bc_cor.pdf')
plt.figure()
print 'PAXB'
print plot_cor(PAXB_VTP_expdiff,PAXB_CERC_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='blue')
print plot_cor(PAXB_RABS_bc_expdiff,PAXB_CERC_bc_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='red')
plt.savefig(fig_base + 'PAXB_RABS_PAXB_CERC_bc_cor.pdf')



ct_labels = ['Opposing cis','Opposing trans','Concordant cis','Concordant trans']
CERC_ct = np.array([2587., 3934., 806., 986.])
CERC_ct = CERC_ct/np.sum(CERC_ct)
PAXB_ct = np.array([2281. ,5108., 1283., 2337.])
PAXB_ct = PAXB_ct/np.sum(PAXB_ct)
ys = np.arange(len(ct_labels))
plt.figure()
plt.barh(ys, CERC_ct, align='center',color='g')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'CERC_cis_trans_bar_all.pdf')
plt.savefig(fig_base + 'CERC_cis_trans_bar_all.svg')
plt.show()

plt.barh(ys, PAXB_ct, align='center',color='b')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'PAXB_cis_trans_bar_all.pdf')
plt.savefig(fig_base + 'PAXB_cis_trans_bar_all.svg')
plt.show()
'''



'''
rval_list = []
for i in range(10000):

    rand_keys = random.sample(PAXB_VTP_expdiff.keys(),bc_length)
    rand_keys_dict = {}
    for rkey in rand_keys:
        rand_keys_dict[rkey] = PAXB_VTP_expdiff[rkey]
    r,pval =  plot_cor(rand_keys_dict,CERC_VTP_expdiff)
    rval_list.append(r)

rval_list.sort()
print rval_list
print plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff,XYLIM=[-7,7],DENSITY=False,COLOR='red')
'''
'''
print plot_cor(PAXB_BS_ase,PAXB_VTP_ase)
plt.show()
print plot_cor(PAXB_VTP_ase,CERC_VTP_ase)
plt.show()
print plot_cor(PAXB_BS_ase,CERC_VTP_ase)
plt.show()
print plot_cor(PAXB_BS_sig_ase,PAXB_VTP_sig_ase)
plt.show()
print plot_cor(PAXB_VTP_sig_ase,CERC_VTP_sig_ase)
plt.show()
print plot_cor(PAXB_BS_sig_ase,CERC_VTP_sig_ase)
plt.show()
'''
##supFig - ase vs expDiff
'''
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/SupFigs/'
r,pval = plot_cor(PAXB_VTP_expdiff_sig,PAXB_VTP_ase)
plt.xlabel('Expression Divergence - PAXB vs RABS, log2')
plt.ylabel('cis - PAXB vs RABS, log2')
plt.savefig(fig_base + 'PAXB_expDiv_cis_cor.pdf')
plt.show()
r,pval = plot_cor(PAXB_VTP_expdiff_sig,PAXB_trans)
plt.xlabel('Expression Divergence - PAXB vs RABS, log2')
plt.ylabel('trans - PAXB vs RABS, log2')
plt.savefig(fig_base + 'PAXB_expDiv_trans_cor.pdf')
plt.show()
r,pval = plot_cor(PAXB_VTP_expdiff_sig,CERC_trans)
plt.xlabel('Expression Divergence - PAXB vs RABS, log2')
plt.ylabel('trans - CERC vs RABS, log2')
plt.savefig(fig_base + 'PAXB_expDiv_CERCtrans_cor.pdf')
plt.show()
r,pval = plot_cor(CERC_VTP_expdiff_sig,CERC_VTP_ase)
plt.xlabel('Expression Divergence - CERC vs RABS, log2\n r=%s' % r)
plt.ylabel('cis - CERC vs RABS, log2')
plt.savefig(fig_base + 'CERC_expDiv_cis_cor.pdf')
plt.show()
r,pval = plot_cor(CERC_VTP_expdiff_sig,CERC_trans)
plt.xlabel('Expression Divergence - CERC vs RABS, log2\n r=%s' % r)
plt.ylabel('trans - CERC vs RABS, log2')
plt.savefig(fig_base + 'CERC_expDiv_trans_cor.pdf')
plt.show()
r,pval = plot_cor(CERC_VTP_expdiff_sig,PAXB_trans)
plt.xlabel('Expression Divergence - PAXB vs RABS, log2\n r=%s' % r)
plt.ylabel('trans - PAXB vs RABS, log2')
plt.savefig(fig_base + 'CERC_expDiv_PAXBtrans_cor.pdf')
plt.show()
'''
'''
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/SupFigs/'
print plot_cor(sig_PAXB_percent_cis,mean_exp)
plt.xlim([0,1])
plt.ylim([0,3000])
plt.savefig(fig_base + 'PAXB_ase_sigexp_cor.pdf')
plt.show()
print plot_cor(sig_CERC_percent_cis,mean_exp)
plt.xlim([0,1])
plt.ylim([0,3000])
plt.savefig(fig_base + 'CERC_ase_sigexp_cor.pdf')
plt.show()
'''
'''
print 'BC', len(PAXB_RABS_bc_expdiff)
import random
r_list = []
for i in range(100):
    random_dict =  filter_dict(PAXB_VTP_expdiff,random.sample(PAXB_VTP_expdiff,439))
    r,pval = plot_cor(random_dict,CERC_VTP_expdiff)
    r_list.append(r)

r_list.sort(reverse=True)
print r_list

print plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff)
plt.show()
##plot ASE by tissue
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Fig6/'
F1_exp = load_exp_table_reps('/home/james/Dropbox/Miller/data/RNA_seq_exp/F1all_norm_ids.fpkm_table',REPS=5)
PAXB_BS_mean_exp = mean_exp_cond_cor(F1_exp,2,LOG=0)
PAXB_VTP_mean_exp = mean_exp_cond_cor(F1_exp,0,LOG=0)
CERC_VTP_mean_exp = mean_exp_cond_cor(F1_exp,1,LOG=0)
print sum(PAXB_VTP_ase.values())/len(PAXB_VTP_ase.values())
print stats.wilcoxon(PAXB_VTP_ase.values())
print stats.ttest_1samp(PAXB_VTP_ase.values(),0)
print sum(CERC_VTP_ase.values())/len(CERC_VTP_ase.values())
print stats.wilcoxon(CERC_VTP_ase.values())
print stats.ttest_1samp(CERC_VTP_ase.values(),0)


plot_cor(CERC_VTP_ase,PAXB_VTP_ase)
plt.show()
plot_cor(CERC_VTP_ase,PAXB_VTP_ase)
plt.show()
plot_cor(PAXB_VTP_mean_exp,CERC_VTP_mean_exp)
plt.show()
'''
##tissue_ase_exp(PAXB_VTP_mean_exp,PAXB_BS_mean_exp,PAXB_VTP_sig_ase,PAXB_BS_sig_ase,DENSITY=True,FIG_BASE=fig_base + 'PAXB_VTP_PAXB_BS_')
##tissue_ase_exp(CERC_VTP_mean_exp,PAXB_BS_mean_exp,CERC_VTP_sig_ase,PAXB_BS_sig_ase,DENSITY=True,FIG_BASE=fig_base + 'CERC_VTP_PAXB_BS_')
##tissue_ase_exp(PAXB_VTP_mean_exp,CERC_VTP_mean_exp,PAXB_VTP_sig_ase,CERC_VTP_sig_ase,DENSITY=True,FIG_BASE=fig_base + 'PAXB_VTP_CERC_VTP_')

def blizzard_plot(inarray,LABEL,COLOR='b'):
    steps = len(inarray)

    ys = np.arange(0,steps)/float(steps)
    plt.plot(inarray,ys,marker=',',label=LABEL,color=COLOR)



##plotSmoothedHist(sig_CERC_percent_cis.values(),'CERC',x_range=[0,1])
'''
plt.figure(figsize=(10,1))
blizzard_plot(sorted(sig_PAXB_percent_cis.values()),'PAXB',COLOR='b')
blizzard_plot(sorted(sig_CERC_percent_cis.values()),'CERC',COLOR='g')

plt.legend()
print stats.mannwhitneyu(sig_PAXB_percent_cis.values(),sig_CERC_percent_cis.values())
plt.show()
'''
def filter_for_bitecode(indict):
    snames = module.ENSidDict()
    toret = {}
    bc_genes = loadBiteCode()
    bc_set = set()
    for gene in bc_genes:
        bc_set.add(gene.upper())
    print len(bc_set)
    for ensid in indict.keys():
        cur_sname = snames[ensid].split()[0].upper()
        if (cur_sname + 'A' in bc_set):
            cur_sname = cur_sname + 'A'
        if (cur_sname + 'B' in bc_set):
            cur_sname = cur_sname + 'B'
        if (cur_sname in bc_set):
            toret[ensid] = indict[ensid]
    return toret

##writes a table of ASE results in the form
#name   ensid   logAI   sig(Y/N)
def write_ASE_table(ai_dict,sig_ai_dict,outfile):

    w = open(outfile,'w')
    snames = module.ENSidDict()

    table_base = '{0}\t{1}\t{2}\t{3}'
    header = table_base.format('Name','ENSID','log2(F/M)','isSig?')
    w.write(header + '\n')
    for cur_ensid in ai_dict:
        cur_name = snames[cur_ensid]
        cur_AI = ai_dict[cur_ensid]
        isSig = 'N'
        if cur_ensid in sig_ai_dict:
            isSig = 'Y'
        newrow = table_base.format(cur_name,cur_ensid,cur_AI,isSig)
        w.write(newrow + '\n')
    w.flush()
    w.close()

write_ASE_table(PAXB_VTP_expdiff,PAXB_VTP_expdiff_sig,'/home/james/Dropbox/Miller/data/ASE/PAXB_RABS_VTP_expDiff.tsv')
write_ASE_table(CERC_VTP_expdiff,CERC_VTP_expdiff_sig,'/home/james/Dropbox/Miller/data/ASE/CERC_RABS_VTP_expDiff.tsv')
write_ASE_table(PAXB_CERC_expdiff,PAXB_CERC_expdiff_sig,'/home/james/Dropbox/Miller/data/ASE/PAXB_CERC_VTP_expDiff.tsv')

write_ASE_table(PAXB_VTP_ase,PAXB_VTP_sig_ase,'/home/james/Dropbox/Miller/data/ASE/PAXB_RABS_VTP_AI.tsv')
write_ASE_table(CERC_VTP_ase,CERC_VTP_sig_ase,'/home/james/Dropbox/Miller/data/ASE/CERC_RABS_VTP_AI.tsv')

'''
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/SupFigs/'
plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_expdiff_filtered,CERC_VTP_expdiff_filtered,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'expDiff', r,pval
#bc_r,bc_pval = plot_cor(filter_for_bitecode(PAXB_VTP_expdiff_filtered),filter_for_bitecode(CERC_VTP_expdiff_filtered),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')
#print 'expDiffBC',bc_r,bc_pval

plt.title('PAXB_CERC_expDiff_correlation')
#plt.xlabel('PAXB expression change: log2(PAXB/RABS)\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC expression change: log2(CERC/RABS)')
plt.axhline(c='black')
plt.axvline(c='black')
#plt.savefig(fig_base + 'CERC_PAXB_shared_expDiff_all_withBC.pdf')

plt.figure(figsize=(10,10))

r,pval =  plot_cor(PAXB_trans,CERC_trans,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'trans',r,pval
#bc_r,bc_pval =  plot_cor(filter_for_bitecode(PAXB_trans),filter_for_bitecode(CERC_trans),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')

#print 'transBC',bc_r,bc_pval
plt.title('PAXB_CERC_trans_correlation')
#plt.xlabel('PAXB trans change: log2\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC trans change: log2')
plt.axhline(c='black')
plt.axvline(c='black')
#plt.savefig(fig_base + 'CERC_PAXB_shared_trans_all_withBC.pdf')

plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_VTP_ase,CERC_VTP_ase,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'cis',r,pval
#bc_r,bc_pval =  plot_cor(filter_for_bitecode(PAXB_VTP_ase),filter_for_bitecode(CERC_VTP_ase),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')
#print 'cisBC',bc_r,bc_pval
plt.title('PAXB_CERC_cis_correlation')
#plt.xlabel('PAXB cis change: log2(PAXB_F1/RABS_F1)\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC cis change: log2(CERC_F1/RABS_F1)')
plt.axhline(c='black')
plt.axvline(c='black')
#plt.savefig(fig_base + 'CERC_PAXB_shared_cis_all_withBC.pdf')
plt.show()

plt.figure(figsize=(10,10))
plot_percent_cis(PAXB_percent_cis,CERC_percent_cis,'PAXB','CERC')
plt.xlabel('Percent Cis, Total Displacement')
plt.ylabel('Density')
plt.xlim([0,1])
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_all_total.pdf')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_all_total.svg')
plt.show()
'''


plt.figure(figsize=(10,10))
sig_PAXB_percent_cis_bc = filter_for_bitecode(sig_PAXB_percent_cis)
sig_CERC_percent_cis_bc = filter_for_bitecode(sig_CERC_percent_cis)
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Fig4/'
print len(PAXB_VTP_expdiff)
print len(CERC_VTP_expdiff)
print len(sig_PAXB_percent_cis)
print len(sig_PAXB_percent_cis_bc)
print len(sig_CERC_percent_cis)
print len(sig_CERC_percent_cis_bc)

PAXB_percent_cis_vals =  np.array(sig_PAXB_percent_cis.values())
print np.median(PAXB_percent_cis_vals)
print np.median(np.array(sig_PAXB_percent_cis_bc.values()))

CERC_percent_cis_vals =  np.array(sig_CERC_percent_cis.values())
print np.median(CERC_percent_cis_vals)
print np.median(np.array(sig_CERC_percent_cis_bc.values()))

print 'Opposing cis','Opposing trans','Concordant cis','Concordant trans'
calc_trans(PAXB_VTP_expdiff_sig,PAXB_VTP_ase)
a,b = calc_trans(filter_for_bitecode(PAXB_VTP_expdiff_sig),PAXB_VTP_ase)
c,d = calc_trans(CERC_VTP_expdiff_sig,CERC_VTP_ase)
c,d = calc_trans(filter_for_bitecode(CERC_VTP_expdiff_sig),CERC_VTP_ase)
print len(a),len(b),len(c),len(d)


##new plots, with BC lists fig 5+6
blizzard_plot(sorted(sig_PAXB_percent_cis.values()),'PAXB',COLOR='blue')
blizzard_plot(sorted(sig_PAXB_percent_cis_bc.values()),'PAXB_BC',COLOR='purple')
blizzard_plot(sorted(sig_CERC_percent_cis.values()),'CERC',COLOR='green')
blizzard_plot(sorted(sig_CERC_percent_cis_bc.values()),'CERC_BC',COLOR='orange')

print 'PAXB,CERC'
print stats.mannwhitneyu(sig_PAXB_percent_cis.values(),sig_CERC_percent_cis.values())


print stats.mannwhitneyu(sig_PAXB_percent_cis.values(),sig_PAXB_percent_cis_bc.values())
print stats.mannwhitneyu(sig_PAXB_percent_cis_bc.values(),sig_PAXB_percent_cis.values())
print stats.mannwhitneyu(sig_CERC_percent_cis.values(),sig_CERC_percent_cis_bc.values())
plt.xlabel('Percent cis')
plt.ylabel('Cumulative Proportion of genes')
plt.legend()
plt.savefig(fig_base + 'Cum_percent_cis_withBC.pdf')
plt.show()

'''
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Fig5/'
plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_expdiff_sig_filtered,CERC_VTP_expdiff_sig_filtered,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'expDiff', r,pval
bc_r,bc_pval = plot_cor(filter_for_bitecode(PAXB_VTP_expdiff_sig_filtered),filter_for_bitecode(CERC_VTP_expdiff_sig_filtered),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')
print 'expDiffBC',bc_r,bc_pval

plt.title('PAXB_CERC_expDiff_correlation')
plt.xlabel('PAXB expression change: log2(PAXB/RABS)\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC expression change: log2(CERC/RABS)')
plt.axhline(c='black')
plt.axvline(c='black')
plt.savefig(fig_base + 'CERC_PAXB_shared_expDiff_withBC.pdf')

plt.figure(figsize=(10,10))

r,pval =  plot_cor(sig_PAXB_trans,sig_CERC_trans,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'trans',r,pval
bc_r,bc_pval =  plot_cor(filter_for_bitecode(sig_PAXB_trans),filter_for_bitecode(sig_CERC_trans),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')

print 'transBC',bc_r,bc_pval
plt.title('PAXB_CERC_trans_correlation_sig')
plt.xlabel('PAXB trans change: log2\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC trans change: log2')
plt.axhline(c='black')
plt.axvline(c='black')
plt.savefig(fig_base + 'CERC_PAXB_shared_trans_withBC.pdf')

plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_VTP_ase_sigExpdiff,CERC_VTP_ase_sigExpdiff,XYLIM=[-7,7],DENSITY=True,COLOR='k')
print 'cis',r,pval
bc_r,bc_pval =  plot_cor(filter_for_bitecode(PAXB_VTP_ase_sigExpdiff),filter_for_bitecode(CERC_VTP_ase_sigExpdiff),XYLIM=[-7,7],DENSITY=False,COLOR='k',SYM='*')
print 'cisBC',bc_r,bc_pval
plt.title('PAXB_CERC_cis_correlation_sig')
plt.xlabel('PAXB cis change: log2(PAXB_F1/RABS_F1)\nr = %6.3f bc_r = %6.3f'% (r,bc_r))
plt.ylabel('CERC cis change: log2(CERC_F1/RABS_F1)')
plt.axhline(c='black')
plt.axvline(c='black')
plt.savefig(fig_base + 'CERC_PAXB_shared_cis_withBC.pdf')
plt.show()

'''



'''
plot_cor(PAXB_VTP_mean_exp,PAXB_BS_mean_exp)
plt.show()
plot_cor(CERC_VTP_mean_exp,PAXB_BS_mean_exp)
plt.show()
plot_cor(PAXB_VTP_mean_exp,CERC_VTP_mean_exp)
plt.show()
'''

fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Fig4/'


##make plots for Fig4
'''
plot_cis_trans_v2(PAXB_VTP_ase,PAXB_VTP_expdiff,PAXB_VTP_sig_ase,PAXB_trans_sigOnly)
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.xlabel('Hybrid Expression Ratio (log2)')
plt.ylabel('Parental Expression Ratio (log2)')
plt.savefig(fig_base + 'PAXB_wittkopp.pdf')
plt.show()


plot_cis_trans_v2(CERC_VTP_ase,CERC_VTP_expdiff,CERC_VTP_sig_ase,CERC_trans_sigOnly)
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.xlabel('Hybrid Expression Ratio (log2)')
plt.ylabel('Parental Expression Ratio (log2)')
plt.savefig(fig_base + 'CERC_wittkopp.pdf')
plt.show()


plotSmoothedHist(PAXB_VTP_ase.values(),'PAXB',x_range=[-3,3])
plotSmoothedHist(CERC_VTP_ase.values(),'CERC',x_range=[-3,3])
print sum(PAXB_VTP_ase.values())/len(PAXB_VTP_ase.values()),stats.ttest_1samp(PAXB_VTP_ase.values(),0)
print sum(CERC_VTP_ase.values())/len(CERC_VTP_ase.values()), stats.ttest_1samp(CERC_VTP_ase.values(),0)
plt.legend()
plt.savefig(fig_base + 'CERC_PAXB_ASE_density.pdf')
plt.show()

##make cis/trans comparison for genes with sigExpDiff
ct_labels = ['Opposing cis','Opposing trans','Concordant cis','Concordant trans']
CERC_ct = np.array([289., 807., 326., 536.])
##stats.fisher_exact([ [333,303],[116,165]]) = (1.5632468419255718, 0.0020716688232271271)
##stats.binom_test( (449,468)) = .552

CERC_ct = CERC_ct/np.sum(CERC_ct)
PAXB_ct = np.array([391., 2035., 715., 1711.  ])
##stats.fisher_exact([ [1162,1183],[198,439]]) = (2.177813639352101, 5.0068318038781472e-17)
##stats.binom_test( (1360,1662)) = 4.2572384455236682e-08
PAXB_ct = PAXB_ct/np.sum(PAXB_ct)
ys = np.arange(len(ct_labels))
plt.barh(ys, CERC_ct, align='center',color='g')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'CERC_cis_trans_bar.pdf')
plt.savefig(fig_base + 'CERC_cis_trans_bar.svg')
plt.show()

plt.barh(ys, PAXB_ct, align='center',color='b')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'PAXB_cis_trans_bar.pdf')
plt.savefig(fig_base + 'PAXB_cis_trans_bar.svg')
plt.show()

##plot net percetn cis (change for total)
plt.figure(figsize=(10,10))
plot_percent_cis(sig_PAXB_percent_cis,sig_CERC_percent_cis,'PAXB','CERC')
plt.xlabel('Percent Cis, Total Displacement')
plt.ylabel('Density')
plt.xlim([0,1])
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_total.pdf')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_total.svg')
plt.show()

plt.figure(figsize=(10,10))
plot_percent_cis(sig_PAXB_percent_cis_net,sig_CERC_percent_cis_net,'PAXB','CERC')
plt.xlabel('Percent Cis, Net Displacement')
plt.ylabel('Density')
plt.xlim([-2,2])
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_net.pdf')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_net.svg')
plt.show()






##plots for Fig5

fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/Main_figs/Fig5/'
shared_labels = ['Concordant Up','Discordant Up/Down','Discordant Down/Up','Concordant Down']
Shared_ExpDiff = [523,58,144,500]
Shared_Trans = [490,117,190,428]
Shared_Cis = [302,257,289,377]
xs = np.arange(len(shared_labels))
plt.bar(xs,Shared_ExpDiff,color='grey',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'ExpDiff_shared_sig_bar.pdf')
plt.show()

plt.bar(xs,Shared_Trans,color='grey',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'Trans_shared_sig_bar.pdf')
plt.show()

plt.bar(xs,Shared_Cis,color='grey',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'Cis_shared_sig_bar.pdf')
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_expdiff_filtered,CERC_VTP_expdiff_filtered,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_expDiff_correlation')
plt.xlabel('PAXB expression change: log2(PAXB/RABS)\nr = %6.3f'%r)
plt.ylabel('CERC expression change: log2(CERC/RABS)')
plt.savefig(fig_base + 'PAXB_CERC_expDiff_cor.pdf')
plt.show()
plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_trans,CERC_trans,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_trans_correlation')
plt.xlabel('PAXB trans change: log2\nr = %6.3f'%r)
plt.ylabel('CERC trans change: log2')
plt.savefig(fig_base + 'PAXB_CERC_trans_cor.pdf')
plt.show()
plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_VTP_ase,CERC_VTP_ase,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_cis_correlation')
plt.xlabel('PAXB cis change: log2(PAXB_F1/RABS_F1)\nr = %6.3f'%r)
plt.ylabel('CERC cis change: log2(CERC_F1/RABS_F1)')
plt.savefig(fig_base + 'PAXB_CERC_cis_cor.pdf')
plt.show()


##this looks like a fantastic way to show trans effects dominate
##make shared cis into suppFig??
plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_VTP_expdiff_sig_filtered,CERC_VTP_expdiff_sig_filtered,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_expDiff_correlation_sig')
plt.xlabel('PAXB expression change: log2(PAXB/RABS)\nr = %6.3f'%r)
plt.ylabel('CERC expression change: log2(CERC/RABS)')
plt.savefig(fig_base + 'PAXB_CERC_expDiff_cor_sig.pdf')
plt.show()
plt.figure(figsize=(10,10))
r,pval =  plot_cor(sig_PAXB_trans,sig_CERC_trans,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_trans_correlation_sig')
plt.xlabel('PAXB trans change: log2\nr = %6.3f'%r)
plt.ylabel('CERC trans change: log2')
plt.savefig(fig_base + 'PAXB_CERC_trans_cor_sig.pdf')
plt.show()
plt.figure(figsize=(10,10))
r,pval =  plot_cor(PAXB_VTP_ase_sigExpdiff,CERC_VTP_ase_sigExpdiff,XYLIM=[-7,7],DENSITY=True,COLOR='r')
plt.title('PAXB_CERC_cis_correlation_sig')
plt.xlabel('PAXB cis change: log2(PAXB_F1/RABS_F1)\nr = %6.3f'%r)
plt.ylabel('CERC cis change: log2(CERC_F1/RABS_F1)')
plt.savefig(fig_base + 'PAXB_CERC_cis_cor_sig.pdf')
plt.show()
'''
'''
##supp fig ASE correleation plots
print plot_cor(CERC_VTP_sig_ase,PAXB_VTP_sig_ase,XYLIM=[-7,7])
plt.show()
print plot_cor(PAXB_BS_sig_ase,PAXB_VTP_sig_ase,XYLIM=[-7,7])
plt.show()
print plot_cor(PAXB_BS_sig_ase,FTC_BS_sig_ase,XYLIM=[-7,7])
plt.show()
print plot_cor(PAXB_VTP_sig_ase,FTC_BS_sig_ase,XYLIM=[-7,7])
plt.show()
print plot_cor(CERC_VTP_sig_ase,FTC_BS_sig_ase,XYLIM=[-7,7])
plt.show()
'''

'''
##plot_cis_trans(fpath_base + 'PAXB_VTP_logAI.l',fpath_base + 'PAXB_RABS_logfc.l',fpath_base + 'PAXB_RABS_sig_logfc.l',fpath_base + 'PAXB_VTP_logAI_sig.l')
##wittkopp plot
plot_cis_trans(fpath_base + 'CERC_VTP_logAI.l',fpath_base + 'CERC_RABS_logfc.l',fpath_base + 'CERC_RABS_sig_logfc.l',fpath_base + 'CERC_VTP_logAI_sig.l')
plt.xlabel('log2( CERC(FW)_p / RABS(M)_p )\n',fontsize=26)
plt.ylabel('log2( CERC(FW)_F1 / RABS(M)_F1 )\n',fontsize=26)
plt.title('',fontsize=40)
plt.savefig(fig_base + 'CERC_wittkopp.pdf')
plt.savefig(fig_base + 'CERC_wittkopp.svg')
plt.show()



##find trans and shared values for a bar graph
print 'CERC trans all'
CERC_trans,CERC_percent_cis = calc_trans(CERC_VTP_expdiff,CERC_VTP_ase)
print 'CERC trans - ase Sig'
CERC_trans_sig,CERC_sig_percent_cis = calc_trans(CERC_VTP_expdiff,CERC_VTP_sig_ase)
print 'CERC trans - ExpDiff sig'
sig_CERC_trans,sig_CERC_percent_cis = calc_trans(CERC_VTP_expdiff_sig,CERC_VTP_ase)

print 'PAXB trans all'
PAXB_trans,PAXB_percent_cis = calc_trans(PAXB_VTP_expdiff,PAXB_VTP_ase)
print 'PAXB trans - ase Sig'
PAXB_trans_sig,PAXB_sig_percent_cis = calc_trans(PAXB_VTP_expdiff,PAXB_VTP_sig_ase)
print 'PAXB trans - ExpDiff sig'
sig_PAXB_trans,sig_PAXB_percent_cis = calc_trans(PAXB_VTP_expdiff_sig,PAXB_VTP_ase)

print 'VTP: PAXB and CERC'
r,pval = plot_cor(PAXB_VTP_sig_ase,CERC_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
print 'BS: PAXB and CERC'
r,pval = plot_cor(PAXB_BS_sig_ase,CERC_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
print 'BS: PAXB and PAXB'
r,pval = plot_cor(PAXB_BS_sig_ase,PAXB_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
'''
##ASE_list_file,trans_exp_diff_file,sig_ase_file,sig_trans_file
##plot net percetn cis (change for total)
##make bars for cis/trans opposing
'''
shared_labels = ['Concordant Up','Discordant Up/Down','Discordant Down/Up','Concordant Down']
VTPp_VTPc = [223, 176, 178 ,267]
BSp_VTPc = [250, 201, 169, 273]
BSp_VTPp = [753, 119, 90, 781]
xs = np.arange(len(shared_labels))
plt.bar(xs,VTPp_VTPc,color='b',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'VTPp_VTPc_shared_sig_ase_bar.pdf')
plt.show()

plt.bar(xs,BSp_VTPc,color='r',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'BSp_VTPc_shared_sig_ase_bar.pdf')
plt.show()

plt.bar(xs,BSp_VTPp,color='g',align='center')
plt.xticks(xs,shared_labels)
plt.ylabel('Number of Genes')
plt.tight_layout()
plt.savefig(fig_base + 'BSp_VTPp_shared_sig_ase_bar.pdf')
plt.show()
'''
'''
ct_labels = ['Opposing cis','Opposing trans','Concordant cis','Concordant trans']
CERC_ct = np.array([1126., 1679., 426., 579.])
CERC_ct = CERC_ct/np.sum(CERC_ct)
PAXB_ct = np.array([1263. ,2970., 816., 1638.])
PAXB_ct = PAXB_ct/np.sum(PAXB_ct)
ys = np.arange(len(ct_labels))
plt.barh(ys, CERC_ct, align='center',color='g')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'CERC_cis_trans_bar.pdf')
plt.savefig(fig_base + 'CERC_cis_trans_bar.svg')
plt.show()

plt.barh(ys, PAXB_ct, align='center',color='b')
plt.yticks(ys,ct_labels)
plt.xlabel('Proportion')
plt.tight_layout()
plt.xlim([0,.5])
plt.savefig(fig_base + 'PAXB_cis_trans_bar.pdf')
plt.savefig(fig_base + 'PAXB_cis_trans_bar.svg')
plt.show()
'''
'''
plt.figure(figsize=(10,10))
plot_percent_cis(sig_PAXB_percent_cis,sig_CERC_percent_cis,'PAXB','CERC')
plt.xlabel('Percent Cis, Total Displacement')
plt.ylabel('Density')
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_total.pdf')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_total.svg')
plt.show()
'''
'''
plot_cis_trans(fpath_base + 'PAXB_VTP_logAI.l',fpath_base + 'PAXB_RABS_logfc.l',fpath_base + 'PAXB_RABS_sig_logfc.l',fpath_base + 'PAXB_VTP_logAI_sig.l')
##wittkopp plot
##plot_cis_trans(fpath_base + 'CERC_VTP_logAI.l',fpath_base + 'CERC_RABS_logfc.l',fpath_base + 'CERC_RABS_sig_logfc.l',fpath_base + 'CERC_VTP_logAI_sig.l')
##plt.xlabel('log2( CERC(FW)_p / RABS(M)_p )\n"trans"',fontsize=26)
##plt.ylabel('log2( CERC(FW)_F1 / RABS(M)_F1 )\n"cis"',fontsize=26)
##plt.title('CERC(FW) Changes in Gene Expression',fontsize=40)
plt.savefig(fig_base + 'PAXB_wittkopp.pdf')
plt.savefig(fig_base + 'PAXB_wittkopp.svg')
plt.show()



##plot correlated expDiff
plt.figure(figsize=(10,10))
allr,pval = plot_cor(PAXB_CERC_expdiff,PAXB_VTP_expdiff,XYLIM=[-7,7]) ##.72
bcr,pval = plot_cor(PAXB_CERC_bc_expdiff,PAXB_RABS_bc_expdiff,XYLIM=[-7,7],COLOR='r') ##.70
plt.xlabel('Freshwater expression change: log(PAXB(fw)/CERC(fw)' +  '\n'+ str(allr) + ',' + str(bcr),fontsize='x-large')
plt.ylabel('PAXB(fw) expression change: log(PAXB(fw)/RABS(m))',fontsize='x-large')
plt.title('PAXBxCERC vs PAXBxRABS, all genes and BC')
plt.savefig(fig_base + 'PxC_PxR_geneExp_BC.pdf')
plt.savefig(fig_base + 'PxC_PxR_geneExp_BC.svg')
plt.show()

plt.figure(figsize=(10,10))
allr,pval = plot_cor(PAXB_CERC_expdiff,CERC_VTP_expdiff,XYLIM=[-7,7]) ##.72
bcr,pval = plot_cor(PAXB_CERC_bc_expdiff,CERC_RABS_bc_expdiff,XYLIM=[-7,7],COLOR='r') ##.70
plt.xlabel('Freshwater expression change: log(PAXB(fw)/CERC(fw)' +  '\n'+ str(allr) + ',' + str(bcr),fontsize='x-large')
plt.ylabel('CERC(fw) expression change: log(CERC(fw)/RABS(m))',fontsize='x-large')
plt.title('PAXBxCERC vs CERCxRABS, all genes and BC')
plt.savefig(fig_base + 'PxC_CxR_geneExp_BC.pdf')
plt.savefig(fig_base + 'PxC_CxR_geneExp_BC.svg')
plt.show()

plt.figure(figsize=(10,10))
allr,pval = plot_cor(PAXB_VTP_expdiff,CERC_VTP_expdiff,XYLIM=[-7,7]) ##.72
bcr,pval = plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff,XYLIM=[-7,7],COLOR='r') ##.70
plt.xlabel('PAXB(fw) expression change: log(PAXB(fw)/RABS(fw)' +  '\n'+ str(allr) + ',' + str(bcr),fontsize='x-large')
plt.ylabel('CERC(fw) expression change: log(CERC(fw)/RABS(m))',fontsize='x-large')
plt.title('PAXBxRABS vs CERCxRABS, all genes and BC')
plt.savefig(fig_base + 'PxR_CxR_geneExp_BC.pdf')
plt.savefig(fig_base + 'PxR_CxR_geneExp_BC.svg')
plt.show()

##plot_cor(PAXB_CERC_expdiff,CERC_VTP_expdiff) ##-.18
##plot_cor(PAXB_VTP_expdiff,CERC_VTP_expdiff) ##.41

##plot_cor(PAXB_CERC_bc_expdiff,PAXB_RABS_bc_expdiff) ##.70
##plot_cor(PAXB_CERC_bc_expdiff,CERC_RABS_bc_expdiff) ##.042
##plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff) ##.68

plt.figure(figsize=(10,10))
plot_cor(PAXB_VTP_expdiff,CERC_VTP_expdiff,COLOR='b') ##.41
plot_cor(PAXB_RABS_bc_expdiff,CERC_RABS_bc_expdiff,COLOR='r') ##.68
plt.xlabel('PAXB(fw) expression change: log(PAXB(fw)/RABS(m))',fontsize='x-large')
plt.ylabel('PAXB(fw) expression change: log(PAXB(fw)/RABS(m))',fontsize='x-large')

plt.show()

plt.figure(figsize=(10,10))
hist_dict_list([PAXB_VTP_ase,CERC_VTP_ase,PAXB_BS_ase],['PxR','CxR','BS'])
plt.xlim([-3,3])
plt.xlabel('AIB')
plt.ylabel('Density')
plt.legend()
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'ASE_hist.pdf')
plt.savefig(fig_base + 'ASE_hist.svg')
plt.show()
##plot_cor(PAXB_VTP_ase,CERC_VTP_ase)
##plot_cor(CERC_VTP_expdiff_sig,PAXB_VTP_expdiff_sig)

##plot_cor(CERC_trans,PAXB_trans)
##plot_cor(CERC_trans_sig,PAXB_trans_sig)

##plot_cor(CERC_VTP_sig_ase,PAXB_VTP_sig_ase)
##plot_cor(CERC_VTP_sig_ase,PAXB_VTP_sig_ase)
##plot_cor(CERC_VTP_sig_ase,PAXB_BS_sig_ase)

##plot shared ASE
plt.figure(figsize=(10,10))
##plot_cor(PAXB_BS_ase,PAXB_VTP_ase,COLOR='b')
r,pval = plot_cor(PAXB_BS_sig_ase,PAXB_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Branchial AIB - PAXB vs RABS' + '\n' + str(r))
plt.ylabel('Toothplate AIB - PAXB vs RABS')
plt.savefig(fig_base + 'BS_VTP_shared_ase.pdf')
plt.savefig(fig_base + 'BS_VTP_shared_ase.svg')
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_BS_sig_ase,CERC_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Branchial AIB - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate AIB - CERC vs RABS')
plt.savefig(fig_base + 'BS_CERC_VTP_shared_ase.pdf')
plt.savefig(fig_base + 'BS_CERC_VTP_shared_ase.svg')
plt.show()

plt.figure(figsize=(10,10))
##plot_cor(PAXB_BS_ase,PAXB_VTP_ase,COLOR='b')
r,pval = plot_cor(PAXB_BS_ase,PAXB_VTP_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Branchial AIB - PAXB vs RABS' + '\n' + str(r))
plt.ylabel('Toothplate AIB - PAXB vs RABS')
plt.savefig(fig_base + 'BS_VTP_shared_ase_all.pdf')
plt.savefig(fig_base + 'BS_VTP_shared_ase_all.svg')
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_BS_ase,CERC_VTP_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Branchial AIB - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate AIB - CERC vs RABS')
plt.savefig(fig_base + 'BS_CERC_VTP_shared_ase_all.pdf')
plt.savefig(fig_base + 'BS_CERC_VTP_shared_ase_all.svg')
plt.show()


plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_sig_ase,CERC_VTP_sig_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Toothplate AIB - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate AIB - CERC vs RABS')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_ase.pdf')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_ase.svg')
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_trans_sig,CERC_trans_sig,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Toothplate trans effects - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate trans effects - CERC vs RABS')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_trans.pdf')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_trans.svg')
plt.show()


plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_ase,CERC_VTP_ase,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Toothplate AIB - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate AIB - CERC vs RABS')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_ase_all.pdf')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_ase_all.svg')
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_trans,CERC_trans,DENSITY=True,XYLIM=[-4,4])
plt.xlabel('Toothplate trans effects - PAXB vs RABS'+ '\n' + str(r))
plt.ylabel('Toothplate trans effects - CERC vs RABS')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_trans_all.pdf')
plt.savefig(fig_base + 'PAXB_VTP_CERC_VTP_shared_trans_all.svg')
plt.show()

##plot_cor(PAXB_VTP_expdiff_sig,PAXB_trans_sig)
##plot_cor(PAXB_VTP_expdiff_sig,PAXB_VTP_sig_ase)
##plot_cor(PAXB_VTP_ase,PAXB_trans)

##plot_cor(PAXB_BS_ase,PAXB_BS_Dref_ase)
##plot_cor(PAXB_BS_sig_ase,PAXB_BS_Dref_sig_ase)

##plot net percetn cis (change for total)
plt.figure(figsize=(10,10))
plot_percent_cis(sig_PAXB_percent_cis,sig_CERC_percent_cis,'PAXB','CERC')
plt.xlabel('Percent Cis, Net Displacement')
plt.ylabel('Density')
plt.axvline(x=0,ls='dashed',color='black')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_net.pdf')
plt.savefig(fig_base + 'CERC_PAXB_percentCis_net.svg')
plt.show()


plot_cor(PAXB_VTP_sig_ase,PAXB_VTP_expdiff)
plt.show()
plot_cor(PAXB_VTP_sig_ase,PAXB_trans_sig)
plt.show()
plot_cor(PAXB_trans_sig,PAXB_VTP_expdiff)
plt.show()

plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_ase,PAXB_VTP_expdiff,DENSITY=True)
plt.xlabel('PAXB ToothPlate AIB, PAXB/RABS' + '\n' + str(r),fontsize='x-large')
plt.ylabel('PAXB ToothPlate expression difference, log_2(PAXB/RABS)',fontsize='x-large')
plt.tight_layout()
plt.savefig(fig_base + 'PAXB_AIB_exp_cor.pdf')
plt.savefig(fig_base + 'PAXB_AIB_exp_cor.svg')
plt.show()
plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_VTP_ase,PAXB_trans,DENSITY=True)
plt.xlabel('PAXB ToothPlate AIB, PAXB/RABS' + '\n' + str(r),fontsize='x-large')
plt.ylabel('PAXB ToothPlate trans difference, log_2(PAXB/RABS)',fontsize='x-large')
plt.tight_layout()
plt.savefig(fig_base + 'PAXB_AIB_trans_cor.pdf')
plt.savefig(fig_base + 'PAXB_AIB_trans_cor.svg')
plt.show()
plt.figure(figsize=(10,10))
r,pval = plot_cor(PAXB_trans,PAXB_VTP_expdiff,DENSITY=True)
plt.xlabel('PAXB ToothPlate trans difference, log_2(PAXB/RABS)' + '\n' + str(r),fontsize='x-large')
plt.ylabel('PAXB ToothPlate expression difference, log_2(PAXB/RABS)',fontsize='x-large')
plt.tight_layout()
plt.savefig(fig_base + 'PAXB_trans_exp_cor.pdf')
plt.savefig(fig_base + 'PAXB_trans_exp_cor.svg')
plt.show()
'''
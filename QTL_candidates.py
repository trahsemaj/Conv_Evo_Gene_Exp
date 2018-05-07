__author__ = 'james'

import pybedtools
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib,math,module
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rc('font', family='sans-serif')
matplotlib.rc('font', serif='Arial')
matplotlib.rc('text', usetex='false')
matplotlib.rc('pdf', fonttype=42)



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

##loads the chr,start, and stop of the TSS
def load_TSS(fpath='/home/james/Dropbox/Miller/data/BED/ensGene_collapsed_TSS.bed'):


    f = open(fpath,'r')

    chrpos_dict = {}
    snames = module.ENSidDict()
    for line in f:
        line = line.strip()
        split = line.split('\t')
        chrome = split[0]
        start = split[1]
        stop = split[2]

        chrpos = chrome + ':' + start + '-' + stop

        ensid = split[3]

        chrpos_dict[ensid] = chrpos
    return chrpos_dict


##takes a gene_exp.diff file
##returns a dict of all ensids, with values as a list of every ensid with the same xloc (compare expDiff ensids to collapsed ensids)
def load_xloc_synonyms(in_diff):


    syn_dict = {}

    f = open(in_diff,'r')

    for line in f:
        line = line.strip()
        split = line.split('\t')
        cur_ensids = split[2].split(',')

        for ensid in cur_ensids:
            syn_dict[ensid] = cur_ensids

    return syn_dict


##writes a BED output for the given score dict, giving each ensid a TSS and a name
##will match using the name of the ensid instead of the ensid itself if the MATCH_NAMES is set to True
def write_diff_bed(score_dict,outbed,TSS_PATH='/home/james/Dropbox/Miller/data/BED/ensGene_collapsed_TSS.bed',
                   IN_DIFF='/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_LITC_gene_exp.diff',NAMES=True,USE_SYN=True,TRACK='',BGA=False):

    tss_chrpos = load_TSS(fpath=TSS_PATH)
    snames = module.ENSidDict()
    syn_dict = load_xloc_synonyms(IN_DIFF)

    w = open(outbed,'w')
    if TRACK:
        w.write(TRACK + '\n')

    for ensid in score_dict.keys():

        chrpos = ''
        name = ''
        ##use all ensids at each recorded transcribed locus
        if USE_SYN:
            all_syn = syn_dict[ensid]
            for cur_syn in all_syn:
                if cur_syn in tss_chrpos:
                    chrpos = tss_chrpos[cur_syn]
                    name = snames[cur_syn]


        else:
            chrpos = tss_chrpos[ensid]
            name = snames[ensid]

        ##skip ensids not found in the collapsed TSS
        if not chrpos:
            continue



        chrome = chrpos.split(':')[0]
        start = chrpos.split(':')[1].split('-')[0]
        stop = chrpos.split(':')[1].split('-')[1]

        if not NAMES:
            name = ensid

        score = str(score_dict[ensid])


        if not BGA:
            towrite = '\t'.join([chrome,start,stop,name,score])
        else:
            towrite = '\t'.join([chrome,start,stop,score])
        w.write(towrite + '\n')

    w.flush()
    w.close()
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
##takes a key and a list of ase dicts
##plots the
def plot_key_ase_dict(keys,dict_list,LABELS='',COLORS=''):

    keys_ase =[]
    key_names = []
    plt.figure()
    snames = module.ENSidDict()
    for key in keys:
        cur_name = snames[key]
        key_names.append(cur_name)
        gene_ase = []
        for ind,ase_dict in enumerate(dict_list):
            cur_counts = np.array(ase_dict[key])
            print cur_name
            print cur_counts
            if np.count_nonzero(cur_counts[:,0]) < 3 or np.count_nonzero(cur_counts[:,0]) < 3:
                gene_ase.append(np.array([0]))
                continue
            cur_ase_vals = cur_counts[:,0]/cur_counts[:,1]
            cur_ase_vals = np.log2(cur_ase_vals)
            cur_ase_vals = cur_ase_vals[~np.isnan(cur_ase_vals)]
            gene_ase.append(cur_ase_vals)
        keys_ase.append(gene_ase)

    for ind_gene,cur_gene in enumerate(keys_ase):
        for ind_pop,cur_pop_ase in enumerate(cur_gene):
            if len(cur_pop_ase) < 3:
                continue
            cur_index = ind_gene * 3 + ind_pop
            xs = np.zeros(len(cur_pop_ase)) + cur_index
            plt.scatter(xs,cur_pop_ase,marker='.',color='black')
            plt.errorbar(cur_index,cur_pop_ase.mean(),yerr=stats.sem(cur_pop_ase),fmt='o',ecolor=COLORS[ind_pop],color=COLORS[ind_pop])


    zero_xs = np.arange(len(keys)+1) * 3
    zeros = np.zeros(len(zero_xs))
    plt.plot(zero_xs,zeros,color='grey',ls='dashed')
    plt.xlim([-.5,(len(keys)*3)-.5])
    x_ticks = np.arange(len(keys)) * 3
    plt.xticks(zero_xs ,key_names,rotation=45)
    plt.ylabel('log2(F/M)')
    plt.tight_layout()
    #xs = np.zeros(len(cur_ase_vals)) + ind
    #plt.scatter(xs,cur_ase_vals,marker='.',color='black')
    #plt.errorbar(ind,cur_ase_vals.mean(),yerr=stats.sem(cur_ase_vals),fmt='o',ecolor=COLORS[ind])
    '''
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
    '''
##takes a key, an exp_dict, and a series of indexes
##plots the fc from ref_index mean
def plot_key_expDiff_dict(keys,exp_dict,test_indexes,ref_index,LABELS='',COLORS=''):


    plt.figure()
    keys_exp_diff =[]
    key_names = []
    snames = module.ENSidDict()
    for key in keys:
        cur_name = snames[key]
        key_names.append(cur_name)
        gene_exp_diff = []
        for cur_pop in test_indexes:
            cur_exp = np.array(exp_dict[key][cur_pop])
            cur_ref_mean = np.array(exp_dict[key][ref_index]).mean()
            cur_exp_diff = cur_exp/cur_ref_mean
            cur_exp_diff = np.log2(cur_exp_diff)
            gene_exp_diff.append(cur_exp_diff)
        keys_exp_diff.append(gene_exp_diff)


    for ind_gene,cur_gene in enumerate(keys_exp_diff):
        for ind_pop,cur_pop_diff in enumerate(cur_gene):
            cur_index = ind_gene * 3 + ind_pop
            xs = np.zeros(len(cur_pop_diff)) + cur_index
            plt.scatter(xs,cur_pop_diff,marker='.',color='black')
            plt.errorbar(cur_index,cur_pop_diff.mean(),yerr=stats.sem(cur_pop_diff),fmt='o',ecolor=COLORS[ind_pop],color=COLORS[ind_pop])


    zero_xs = np.arange(len(keys)+1) * 3
    zeros = np.zeros(len(zero_xs))
    plt.plot(zero_xs,zeros,color='grey',ls='dashed')
    plt.xlim([-.5,(len(keys)*3)-.5])
    x_ticks = np.arange(len(keys)) * 3
    plt.xticks(zero_xs ,key_names,rotation=45)
    plt.ylabel('log2(F/M)')
    plt.tight_layout()
    ##full_xs = np.resize(xs,(len(xs),5))

    ##plt.scatter(full_xs,plot_vals,marker='.')
    ##plt.errorbar(xs,plot_vals.mean(axis=1),yerr=stats.sem(plot_vals,axis=1),fmt='_')


    #zero_xs = np.arange(len(test_indexes)+1)-.5
    #zeros = np.zeros(len(zero_xs))
    #plt.plot(zero_xs,zeros,color='grey',ls='dashed')
    #plt.xticks(xs,LABELS,rotation=45)
    #plt.xlim([-.5,len(test_indexes)-.5])
    #plt.ylabel('log2(F/M)')



fpath_base = '/home/james/Dropbox/Miller/data/ASE/'
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE/Candidates'



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


parental_exp = load_exp_table_reps('/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table')
PAXB_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_VTP_AI_all.tsv')
PAXB_BS_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/PAXB_BS_AI_all.tsv')
CERC_VTP_ase_dict = load_ase_dict('/home/james/Dropbox/Miller/data/ASE/CERC_VTP_AI_all.tsv')

fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/SupFigs/'
plot_key_expDiff_dict(['ENSGACT00000004032','ENSGACT00000021002','ENSGACT00000008405','ENSGACT00000010071','ENSGACT00000004441','ENSGACT00000008512'],
                      parental_exp,[0,2],1,COLORS=['g','b'],LABELS=['CERC','PAXB'])
plt.savefig(fig_base + 'Taste_Bud_expDiff.pdf')
plt.show()
##calb2a is up in cis in CERC,PAXB vs RABS
##calb2a is up in CERC vs PAXB
##Tas1r1 is down in CERC and PAXB vs RABS
##PLCB2 is up in CERC and PAXB
##PKD2L1 is down in CERC
'''
##Plot all exp data for each gene
fig_base = '/home/james/Dropbox/Miller/figures/Candidates/'

##PITX2, ZNF106, Plod2, snai1a
plot_key_ase_dict(['ENSGACT00000021804','ENSGACT00000012993','ENSGACT00000004902','ENSGACT00000013425'],
                  [CERC_VTP_ase_dict,PAXB_VTP_ase_dict],LABELS=['CERC','PAXB'],COLORS=['g','b'])
plt.savefig(fig_base + 'CERC_candidate_ASE.pdf')
plt.show()
##snai1a - ENSGACT00000013425
##pdgfrl - ENSGACT00000024635
##ZNF106 - ENSGACT00000012993
##Plod2 - ENSGACT00000004902
##bmp1a - ENSGACT00000014468
##ctnnb1 - ENSGACT00000008011
## PITX2 - ENSGACT00000021804

##PITX2, ZNF106, Plod2, snai1a
plot_key_expDiff_dict(['ENSGACT00000021804','ENSGACT00000012993','ENSGACT00000004902','ENSGACT00000013425'],
                      parental_exp,[0,2],1,COLORS=['g','b'],LABELS=['CERC','PAXB'])
plt.savefig(fig_base + 'CERC_candidate_expDiff.pdf')
plt.show()

write_diff_bed(PAXB_VTP_ase,fpath_base + 'Candidates/PAXB_VTP_ase_all.bed',USE_SYN=False)
write_diff_bed(PAXB_VTP_expdiff,fpath_base + 'Candidates/PAXB_VTP_expDiff_all.bed',USE_SYN=True)

write_diff_bed(PAXB_VTP_ase,fpath_base + 'Candidates/PAXB_VTP_ase_all_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="PAXB_VTP_ase_all" description="PAXB_VTP_ase_all"')

write_diff_bed(PAXB_VTP_expdiff,fpath_base + 'Candidates/PAXB_VTP_expDiff_all_ids.bed',USE_SYN=True,NAMES=False,BGA=True,
                TRACK='track type=bedGraph name="PAXB_VTP_expDiff_all" description="PAXB_VTP_expDiff_all"')

write_diff_bed(PAXB_BS_ase,fpath_base + 'Candidates/PAXB_BS_ase_all.bed',USE_SYN=False)
write_diff_bed(PAXB_BS_ase,fpath_base + 'Candidates/PAXB_BS_ase_all_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="PAXB_BS_ase_all" description="PAXB_BS_ase_all"')

write_diff_bed(CERC_VTP_ase,fpath_base + 'Candidates/CERC_VTP_ase_all.bed',USE_SYN=False)
write_diff_bed(CERC_VTP_expdiff,fpath_base + 'Candidates/CERC_VTP_expDiff_all.bed',USE_SYN=True)

write_diff_bed(CERC_VTP_ase,fpath_base + 'Candidates/CERC_VTP_ase_all_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="CERC_VTP_ase_all" description="CERC_VTP_ase_all"')
write_diff_bed(CERC_VTP_expdiff,fpath_base + 'Candidates/CERC_VTP_expDiff_all_ids.bed',USE_SYN=True,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="CERC_VTP_expDiff_all" description="CERC_VTP_expDiff_all"')


write_diff_bed(PAXB_VTP_sig_ase,fpath_base + 'Candidates/PAXB_VTP_ase_sig.bed',USE_SYN=False)
write_diff_bed(PAXB_VTP_expdiff_sig,fpath_base + 'Candidates/PAXB_VTP_expDiff_sig.bed',USE_SYN=True)

write_diff_bed(PAXB_VTP_sig_ase,fpath_base + 'Candidates/PAXB_VTP_ase_sig_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="PAXB_VTP_ase_sig" description="PAXB_VTP_ase_sig"')
write_diff_bed(PAXB_VTP_expdiff_sig,fpath_base + 'Candidates/PAXB_VTP_expDiff_sig_ids.bed',USE_SYN=True,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="PAXB_VTP_expDiff_sig" description="PAXB_VTP_expDiff_sig"')

write_diff_bed(PAXB_BS_sig_ase,fpath_base + 'Candidates/PAXB_BS_ase_sig.bed',USE_SYN=False)
write_diff_bed(PAXB_BS_sig_ase,fpath_base + 'Candidates/PAXB_BS_ase_sig_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="PAXB_BS_ase_sig" description="PAXB_BS_ase_sig"')

write_diff_bed(CERC_VTP_sig_ase,fpath_base + 'Candidates/CERC_VTP_ase_sig.bed',USE_SYN=False)
write_diff_bed(CERC_VTP_expdiff_sig,fpath_base + 'Candidates/CERC_VTP_expDiff_sig.bed',USE_SYN=True)

write_diff_bed(CERC_VTP_sig_ase,fpath_base + 'Candidates/CERC_VTP_ase_sig_ids.bed',USE_SYN=False,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="CERC_VTP_ase_sig" description="CERC_VTP_ase_sig"')
write_diff_bed(CERC_VTP_expdiff_sig,fpath_base + 'Candidates/CERC_VTP_expDiff_sig_ids.bed',USE_SYN=True,NAMES=False,BGA=True,
               TRACK='track type=bedGraph name="CERC_VTP_expDiff_sig" description="CERC_VTP_expDiff_sig"')
'''
PAXB_QTLs_fpath = '/home/james/Dropbox/Miller/data/QTL_candidates/PAXB_nVTP_QTLs.bed'
CERC_QTLs_fpath = '/home/james/Dropbox/Miller/data/QTL_candidates/CERC_nVTP_QTLs.bed'


##snai1a - ENSGACT00000013425
##pdgfrl - ENSGACT00000024635
##ZNF106 - ENSGACT00000012993
##Plod2 - ENSGACT00000004902
##bmp1a - ENSGACT00000014468
##ctnnb1 - ENSGACT00000008011

print 'name','PAXB_ASE','CERC_ASE','sig in PAXB','sig in CERC'
print 'snai1a', PAXB_VTP_ase['ENSGACT00000013425'],CERC_VTP_ase['ENSGACT00000013425'], 'ENSGACT00000013425' in PAXB_VTP_sig_ase,'ENSGACT00000013425' in CERC_VTP_sig_ase
print 'pdgfrl', PAXB_VTP_ase['ENSGACT00000024635'],'NA', 'ENSGACT00000024635' in PAXB_VTP_sig_ase,'ENSGACT00000024635' in CERC_VTP_sig_ase
print 'ZNF106', PAXB_VTP_ase['ENSGACT00000012993'],CERC_VTP_ase['ENSGACT00000012993'], 'ENSGACT00000012993' in PAXB_VTP_sig_ase,'ENSGACT00000012993' in CERC_VTP_sig_ase
print 'Plod2', PAXB_VTP_ase['ENSGACT00000004902'],CERC_VTP_ase['ENSGACT00000004902'], 'ENSGACT00000004902' in PAXB_VTP_sig_ase,'ENSGACT00000004902' in CERC_VTP_sig_ase
print 'bmp1a', PAXB_VTP_ase['ENSGACT00000014468'],'NA', 'ENSGACT00000014468' in PAXB_VTP_sig_ase,'ENSGACT00000014468' in CERC_VTP_sig_ase
print 'ctnnb1', PAXB_VTP_ase['ENSGACT00000008011'],CERC_VTP_ase['ENSGACT00000008011'], 'ENSGACT00000008011' in PAXB_VTP_sig_ase,'ENSGACT00000008011' in CERC_VTP_sig_ase
print 
print 'name','PAXB_expDiff','CERC_expDiff','sig in PAXB','sig in CERC'
print 'snai1a', PAXB_VTP_expdiff['ENSGACT00000013425'],CERC_VTP_expdiff['ENSGACT00000013425'], 'ENSGACT00000013425' in PAXB_VTP_expdiff_sig,'ENSGACT00000013425' in CERC_VTP_expdiff_sig
print 'pdgfrl', PAXB_VTP_expdiff['ENSGACT00000024635'],CERC_VTP_expdiff['ENSGACT00000024635'], 'ENSGACT00000024635' in PAXB_VTP_expdiff_sig,'ENSGACT00000024635' in CERC_VTP_expdiff_sig
print 'ZNF106', PAXB_VTP_expdiff['ENSGACT00000012993'],CERC_VTP_expdiff['ENSGACT00000012993'], 'ENSGACT00000012993' in PAXB_VTP_expdiff_sig,'ENSGACT00000012993' in CERC_VTP_expdiff_sig
print 'Plod2', PAXB_VTP_expdiff['ENSGACT00000004902'],CERC_VTP_expdiff['ENSGACT00000004902'], 'ENSGACT00000004902' in PAXB_VTP_expdiff_sig,'ENSGACT00000004902' in CERC_VTP_expdiff_sig
print 'bmp1a', PAXB_VTP_expdiff['ENSGACT00000014461'],CERC_VTP_expdiff['ENSGACT00000014461'], 'ENSGACT00000014461' in PAXB_VTP_expdiff_sig,'ENSGACT00000014461' in CERC_VTP_expdiff_sig
print 'ctnnb1', PAXB_VTP_expdiff['ENSGACT00000008011'],CERC_VTP_expdiff['ENSGACT00000008011'], 'ENSGACT00000008011' in PAXB_VTP_expdiff_sig,'ENSGACT00000008011' in CERC_VTP_expdiff_sig
print 'pitx2', PAXB_VTP_expdiff['ENSGACT00000021804'],CERC_VTP_expdiff['ENSGACT00000021804'], 'ENSGACT00000021804' in PAXB_VTP_expdiff_sig,'ENSGACT00000021804' in CERC_VTP_expdiff_sig



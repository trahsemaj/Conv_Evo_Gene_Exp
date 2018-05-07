__author__ = 'james'

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



##takes a toothcode genelist
##returns a dictionary of dict[gene_name] = +/-1, depending on direction of effect
def load_ToothCode(infile):
    tc_dict = {}

    f = open(infile, 'r')

    for line in f:
        line = line.strip()
        split = line.split('\t')

        name = split[0].upper()
        sign = split[1]

        direction = 1
        if sign == '-':
            direction = -1

        tc_dict[name] = direction

    return tc_dict

##loads the suptable_1 from Kandyba 2012 into a dict
##keys are Gene Symbols, values are a list of [Bulge vs HG  , Bulge vs Epidermis/ORS  ,Fold-Change hfSCs BMPR1A cKO-RU vs. CON-RU]
def load_Kandyba_array(fpath='/home/james/Dropbox/Miller/data/BED/Kandyba_2012_st01_marray.tsv'):


    array_dict = {}

    f = open(fpath)
    f.readline()
    for line in f:
        line = line.strip()
        split = line.split('\t')

        gene_symbol = split[2].upper()
        bulge_HG = float(split[3])
        bulge_Epi = float(split[4])
        KO_WT = float(split[5])
        sign = 1
        if KO_WT < 0:
            sign = -1
        array_dict[gene_symbol] = sign##[bulge_HG,bulge_Epi,KO_WT]

    return array_dict

##counts the number of concordant and discordant cis and trans changes in the given tc_dict
def calc_HF_pos(cis_dict,trans_dict,tc_dict):

    hnames = module.ENSid_Human_Dict()
    snames = module.ENSidDict()
    shared_ids = set(cis_dict.keys()).intersection(set(trans_dict.keys()))


    ##counts the cis_trans changes (up or down, binary)
    up_up = 0
    up_down = 0
    down_up = 0
    down_down = 0
    for ensid in shared_ids:
        hortho = ''
        if ensid in hnames:
            hortho = hnames[ensid]
        if hortho in tc_dict:

            cis_lfc = cis_dict[ensid]
            trans_lfc = trans_dict[ensid]
            ##print hortho,cis_lfc,trans_lfc
            if cis_lfc > 0:

                if trans_lfc > 0:
                    up_up += 1
                else:
                    up_down += 1
            else:
                if trans_lfc > 0:
                    down_up += 1
                else:
                    down_down += 1

    return up_up,up_down,down_up,down_down



##calculates HF signs of selection for a list of dictionaries
def calc_all_HF(cis_dict,trans_dict,tc_list,name_list=[]):

    shared_ids = set(cis_dict.keys()).intersection(set(trans_dict.keys()))

    total_up_up = 0
    total_up_down = 0
    total_down_up = 0
    total_down_down = 0

    for ensid in shared_ids:

        cis = cis_dict[ensid]

        trans = trans_dict[ensid]


        if cis > 0 :
            if trans > 0:
                total_up_up += 1
            else:
                total_up_down += 1
        else:
            if trans > 0:
                total_down_up += 1
            else:
                total_down_down += 1


    total_up = total_up_up + total_up_down
    total_down = total_down_down + total_down_up
    print 'Total:'
    print total_up,total_down

    for ind,tc_dict in enumerate(tc_list):
        uu,ud,du,dd =  calc_HF_pos(cis_dict,trans_dict,tc_dict)
        test_up = uu + ud
        test_down = dd + du


        odds,pval = stats.fisher_exact([[total_up,total_down],[test_up,test_down]])
        if pval < .05 and name_list:
            print test_up,test_down
            print pval, name_list[ind]
        ##print stats.hypergeom.sf(test_up,total_up + total_down, total_up, test_up + test_down)


##looks for signs of concordant selection on a given pathway
def calc_all_concord(cis_dict,trans_dict,tc_list,name_list=[]):

    shared_ids = set(cis_dict.keys()).intersection(set(trans_dict.keys()))

    total_up_up = 0
    total_up_down = 0
    total_down_up = 0
    total_down_down = 0

    for ensid in shared_ids:

        cis = cis_dict[ensid]

        trans = trans_dict[ensid]


        if cis > 0 :
            if trans > 0:
                total_up_up += 1
            else:
                total_up_down += 1
        else:
            if trans > 0:
                total_down_up += 1
            else:
                total_down_down += 1


    total_concord = total_up_up + total_down_down
    total_discord = total_up_down + total_down_up
    print 'Total:'
    print total_concord,total_discord

    for ind,tc_dict in enumerate(tc_list):
        uu,ud,du,dd =  calc_HF_pos(cis_dict,trans_dict,tc_dict)
        test_concord = uu + dd
        test_discord = ud + du


        odds,pval = stats.fisher_exact([[total_concord,total_discord],[test_concord,test_discord]])
        ##print stats.fisher_exact([[total_concord,total_discord],[test_concord,test_discord]])
        if pval < .05 and name_list:
            print test_concord,test_discord
            print pval, name_list[ind]
        ##print stats.hypergeom.sf(test_up,total_up + total_down, total_up, test_up + test_down)




fpath_base = '/home/james/Dropbox/Miller/data/ASE/'
fig_base = '/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE/'



bmp_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Bmp.tsv')
fgf_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Fgf.tsv')
shh_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Shh.tsv')
wnt_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Wnt.tsv')
act_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Activin.tsv')
tgf_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Tgfb.tsv')
##del tgf_tc['DCN']
notch_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Notch.tsv')
eda_tc = load_ToothCode('/home/james/Dropbox/Miller/data/BED/ToothCode_Eda.tsv')

upSC = load_Kandyba_array()


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


calc_HF_pos(PAXB_VTP_ase,PAXB_trans,fgf_tc)
##Note: all FGF genes are up in cis
##calc_all_HF(PAXB_VTP_ase,PAXB_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])
##calc_all_HF(CERC_VTP_ase,CERC_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])

##nothing survives, but all fgf components are up in cis
##calc_all_HF(PAXB_VTP_sig_ase,PAXB_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])
##calc_all_HF(CERC_VTP_sig_ase,CERC_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])

pathways = module.loadGeneSets(fpath='/home/james/Dropbox/Miller/data/BED/human_gene_sets.gmt')
kegg_pathways = []
kegg_names = []
for key in pathways.keys():
    if key.startswith('KEGG'):
        kegg_pathways.append(set(pathways[key]))
        kegg_names.append(key)

##nothing survives
##calc_all_HF(PAXB_VTP_sig_ase,PAXB_trans,kegg_pathways,name_list=kegg_names)
##calc_all_HF(CERC_VTP_sig_ase,CERC_trans,kegg_pathways,name_list=kegg_names)

##nothing survives
##calc_all_HF(PAXB_VTP_ase,PAXB_trans,kegg_pathways,name_list=kegg_names)
##calc_all_HF(CERC_VTP_ase,CERC_trans,kegg_pathways,name_list=kegg_names)

##nothing survives
##calc_all_concord(PAXB_VTP_ase,PAXB_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])
##calc_all_concord(CERC_VTP_ase,CERC_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])

##nothing survives
##calc_all_concord(PAXB_VTP_sig_ase,PAXB_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])
##calc_all_concord(CERC_VTP_sig_ase,CERC_trans,[bmp_tc,fgf_tc,shh_tc,wnt_tc,act_tc,tgf_tc,notch_tc,eda_tc,upSC])

print len(kegg_names)
##bonferroni is a bitch, though closer this time
calc_all_concord(PAXB_VTP_ase,PAXB_trans,kegg_pathways,name_list=kegg_names)
calc_all_concord(CERC_VTP_ase,CERC_trans,kegg_pathways,name_list=kegg_names)


##nothing survives
##calc_all_concord(PAXB_VTP_sig_ase,PAXB_trans,kegg_pathways,name_list=kegg_names)
##calc_all_concord(CERC_VTP_sig_ase,CERC_trans,kegg_pathways,name_list=kegg_names)
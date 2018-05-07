'''
general plotting functions for ASE data
'''
import module
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np
import math
from scipy import stats
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods
import subprocess as sp
    
##takes an array and a sample_name
##plots the smoothed histogram (density) of the data
def plotSmoothedHist(array,sample_name,x_range=[-2,8],cur_color=''):

    density = gaussian_kde(array)

    xmin = x_range[0]
    xmax = x_range[1]
    print xmin,xmax
    xs = np.linspace(xmin,xmax,200)
    density.covariance_factor = lambda : .25
    density._compute_covariance()
    if not cur_color:
	    plt.plot(xs,density(xs),label=sample_name)
    else:
	    plt.plot(xs,density(xs),label=sample_name,color=cur_color)

##input: a file containing the list of the ASE values for each gene in the genome
##output: a density histogram of the ASE values across the whole genome
def ASE_density(ase_l_fpath):
    
    ase_dict = {}

    f = open(ase_l_fpath,'r')


    for line in f:
        line = line.strip()
        split = line.split('\t')
        ensid = split[0]
        ase_val = float(split[1])
        ase_val = math.log((ase_val+1e-20)/((1-ase_val)+1e-20),2)
        ase_dict[ensid] = ase_val


    all_ase_values = ase_dict.values()
    all_ase_values = np.array(all_ase_values)

    print all_ase_values.mean()
    plotSmoothedHist(all_ase_values,'VTP_ASE',x_range=[-4,4])
    plt.xlabel('log2(F/M)',fontsize='large')
    plt.ylabel('Density',fontsize='large')
    plt.ylim([0,1])
    xs = [0] * len(range(0,8))
    ys = range(0,8)
    plt.plot(xs,ys,ls='dashed',color='red')
    plt.show()


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

##input: a ASE.l file
##output: the dictionary form of the file
def loadASE_names(ase_file):

    toret = {}
    names = module.ENSidDict()
    f = open(ase_file,'r')

    for line in f:
        line = line.strip()
        split = line.split('\t')
        ensid = split[0]
        gname = names[ensid].upper()
        ase = float(split[1])
        toret[gname]=ase

    return toret


##input: a called ASE file, and a true ASE file
##output: the density plot of the ASE error (abs(True_ase - Called_ase))
def ASEerror_density(called_ase,true_ase):

    called_ase_dict = loadASE_list(called_ase)
    true_ase_dict = loadASE_list(true_ase)

    ##to be filled in with error values later (not paired with names)
    ase_errors = []

    ##called_ase genes are a subset of true_ase genes
    for ensid in called_ase_dict.keys():

        if ensid not in true_ase_dict:
            continue
        true_ase = true_ase_dict[ensid]
        called_ase = called_ase_dict[ensid]

        true_ase = math.log( (true_ase/(1-true_ase)),2)
        called_ase = math.log( (called_ase/(1-called_ase)),2)
        ase_err = math.fabs(true_ase - called_ase)
        ase_errors.append(ase_err)

    ase_errors = np.array(ase_errors)
    print ase_errors.mean()

    plotSmoothedHist(ase_errors,'log2(ASE_error)',x_range=[0,5])
    plt.xlabel('log2(ASE Error)')
    plt.ylabel('Density')
    plt.show()

##input: a called ASE file, and a true ASE file
##output: plots ase_error vs trueASE
def ASEerror_ase(called_ase,true_ase):

    called_ase_dict = loadASE_list(called_ase)
    true_ase_dict = loadASE_list(true_ase)

    ##to be filled in with error values later (not paired with names)
    ase_errors = []
    ase_values = []

    ##called_ase genes are a subset of true_ase genes
    for ensid in called_ase_dict.keys():
        true_ase = true_ase_dict[ensid]
        called_ase = called_ase_dict[ensid]


        true_ase = math.log( (true_ase/(1-true_ase)),2)
        called_ase = math.log( (called_ase/(1-called_ase)),2)

        ase_values.append(true_ase)

        ase_err = math.fabs(true_ase - called_ase)
        ase_errors.append(ase_err)

    plt.scatter(ase_values,ase_errors,alpha=.3,marker='o')
    plt.xlabel('log2(ASE)')
    plt.ylabel('log2 ASE_errors')
    plt.ylim([0,3])
    plt.xlim([-3,3])
    plt.show()


##input: a genes.fpkm_tracking file from cufflinks (not normalized with cuffnorm)
##output: a dictionary of gene expression in the form 
def loadFPKM(filepath):
    f = open(filepath,'r')
    
    f.readline()
    
    toret = {}
    for line in f:
        line = line.strip()
        split = line.split('\t')
        chr = split[6].split(':')[0]
        ##Get rid of ones place, sloppy start
        start = split[6].split(':')[1].split('-')[0][:-1]
        end = split[6].split(':')[1].split('-')[1][:-1]
        key = split[0]
        fpkm = float(split[9])
        if fpkm == 0.0:
            fpkm += 1e-25
        
        if split[0].startswith('CUFF'):
            continue
        
        
        toret[key] = fpkm
        
        
    f.close()
    return toret

##input: a cufflinks genes.fpkm_tracking file, a ref ASE and a called ASE file
##output: plots the ase error vs gene expression
def aseErr_expression(genes_fpkm,called_ase_file,true_ase_file):
    
    called_ase_dict = loadASE_list(called_ase_file)
    true_ase_dict = loadASE_list(true_ase_file)

    gene_expression = loadFPKM(genes_fpkm)

    ##to be filled in with error values later (not paired with names)
    ase_errors = []
    ##filled with fpkms of expressed genes, in step with ase_errors
    expression_values = []

    ##called_ase genes are a subset of true_ase genes
    for ensid in called_ase_dict.keys():
        true_ase = true_ase_dict[ensid]
        called_ase = called_ase_dict[ensid]

        fpkm = gene_expression[ensid]

        true_ase = math.log( (true_ase/(1-true_ase)),2)
        called_ase = math.log( (called_ase/(1-called_ase)),2)

        expression_values.append(math.log(fpkm,2))

        ase_err = math.fabs(true_ase - called_ase)
        ase_errors.append(ase_err)

    plt.scatter(expression_values,ase_errors,alpha=.3,marker='o')
    r,pval =  stats.spearmanr(expression_values,ase_errors)
    plt.xlabel('log2(FPKM)\nr=%f'%r)
    plt.ylabel('log2 ASE_errors')
    plt.ylim([0,2.5])
    plt.xlim([-1,12])


    k0 = smooth.NonParamRegression(expression_values, ase_errors, method=npr_methods.SpatialAverage())
    k0.fit()

    xs = np.arange(-1,12,.01)

    plt.plot(xs, k0(xs), linewidth=2)

    plt.show()




##input: a called ASE file and a cufflinks genes.fpkm_tracking file
##plots ase vs expression
def ase_expression(genes_fpkm,called_ase_file):
    called_ase_dict = loadASE_list(called_ase_file)

    gene_expression = loadFPKM(genes_fpkm)

    ##to be filled in with error values later (not paired with names)
    ase_values = []
    ##filled with fpkms of expressed genes, in step with ase_values
    expression_values = []

    ##called_ase genes are a subset of true_ase genes
    for ensid in called_ase_dict.keys():
        called_ase = called_ase_dict[ensid]

        fpkm = gene_expression[ensid]

        called_ase = math.fabs(math.log( (called_ase/(1-called_ase)),2))

        expression_values.append(math.log(fpkm,2))

        ase_values.append(called_ase)

    plt.scatter(expression_values,ase_values,alpha=.3,marker='o')
    r,pval =  stats.spearmanr(expression_values,ase_values)
    plt.xlabel('log2(FPKM)\nr=%f' % r)
    plt.ylabel('log2 ASE')
    plt.title('VTP ASE vs Expression')
    plt.ylim([0,2.5])
    plt.xlim([-1,12])




    k0 = smooth.NonParamRegression(expression_values, ase_values, method=npr_methods.SpatialAverage())
    k0.fit()

    xs = np.arange(-1,12,.01)

    plt.plot(xs, k0(xs), linewidth=2)

    plt.show()


    
def allPathways(called_ase_file,outfile):

    pathways = module.loadGeneSets()
    ase_dict = loadASE_names(called_ase_file)
    w = open(outfile,'w')
    path_data = []
    for poi in pathways.keys():
        ##only look at KEGG pathways
        if not poi.startswith('KEGG'):
            continue
        up_genes = []
        down_genes = []
        poi_genes = pathways[poi]
        ##add expression in sample1, subtract in sample2
        exp_diff = 0
        for gene in poi_genes[:]:
            if gene not in ase_dict:
                poi_genes.remove(gene)

            else:
                bias = ase_dict[gene]
                exp_diff += bias
                ##up in sample1
            if bias > 0:
                up_genes.append((gene,bias))
            else:
                down_genes.append((gene,bias))

        num_dn = len(down_genes)
        num_up = len(up_genes)

        pval = stats.binom_test(num_up,num_dn+num_up,.5)

        ##towrite =  '\t'.join(map(lambda x: str(x),[poi,len(up_genes),len(down_genes), exp_diff,pval]))
        path_data.append([poi,len(up_genes),len(down_genes), exp_diff,pval])
        ##w.write(towrite + '\n')
        path_data.sort(key = lambda x: x[4])
        for ind,row in enumerate(path_data):
            pval = row[4]
            qval = len(path_data)/(ind+ 1) * pval
            towrite =  '\t'.join(map(lambda x: str(x), row + [qval]))
            w.write(towrite + '\n')
            towrite = ''

##Takes cis and trans gene expression changes, and a fold change cutoff for each (log-scaled)
##filters for genes with expression and cis changes in same direction, and writes them to out_gtf
def findQTLcandidates(ASE_list_file,trans_exp_diff_file,out_gtf,ASE_CUTOFF=.5,TRANS_CUTOFF=.5,REF_GTF='/data/James/BED/ensGene.gtf'):
    
    ##candidate gene list
    candidate_ensids = []

    cis_dict = loadASE_list(ASE_list_file)
    ##using the PAXB_RABS file, in the form log2(PAXB/RABS)
    trans_dict = loadASE_list(trans_exp_diff_file)

    

    ##find the common genes twixt the two

    genes = set(cis_dict.keys()).intersection(trans_dict.keys())


    for gene in genes:

        ##transform ratio to log2
        ase_val = cis_dict[gene]
        ase_val = math.log(ase_val/(1-ase_val+1e-25),2)
        trans_val = trans_dict[gene]

        ##if above cutoff in both
        if math.fabs(ase_val) > ASE_CUTOFF and math.fabs(trans_val) > TRANS_CUTOFF:

            ##if in same direction
            if (ase_val > 0 and trans_val > 0) or (ase_val < 0 and trans_val < 0):
                candidate_ensids.append(gene)
    ##clear the file
    sp.call('> {0}'.format(out_gtf),shell=True)

    for gene in candidate_ensids:
	    sp.call('grep {0} {1} >> {2}'.format(gene,REF_GTF,out_gtf),shell=True)

    print len(candidate_ensids)

##takes a gtf file, and prints out a list of gene names from the file
def gtfToGnames(ingtf):

    gene_names = set()

    ensidNames = module.ENSidDict()
    
    f = open(ingtf,'r')

    for line in f:
        line = line.strip()
        split = line.split('\t')

        cur_ensid = split[8].split('"')[1]
        cur_name = ensidNames[cur_ensid]
        gene_names.add(cur_name)

    for name in gene_names:
	    print name



##input: two ase.l files containing ase ratios for a series of genes
##output: the correlation between the two
##log_vals is true if the input values are already in log2 form
def ase_concordance(ase_file_1,ase_file_2,log_vals=False):

    ase_dict_1 = loadASE_list(ase_file_1)
    ase_dict_2 = loadASE_list(ase_file_2)
    ##get a list of genes shared between the two

    genes = set(ase_dict_1.keys()).intersection(ase_dict_2.keys())

    ##a = open('/data/James/AlleleSpecificExpression/VTP/UUDD_CxP_ASE.bed','w')
    ##b = open('/data/James/AlleleSpecificExpression/VTP/DUUD_CxP_ASE.bed','w')

    positions = module.ensidStartStop()
    names = module.ENSidDict()

    ase_list_1 = []
    ase_list_2 = []
    up_up = 0
    down_down = 0
    up_down = 0
    down_up = 0
    
    for gene in genes:

	    ##log2 the ratio value
	    ##print ase_dict_1[gene],ase_dict_2[gene]
        if not log_vals:
            ase_val_1 = math.log((ase_dict_1[gene] + 1e-30)/(1-ase_dict_1[gene]+1e-30),2)
            ase_val_2 = math.log((ase_dict_2[gene] + 1e-30)/(1-ase_dict_2[gene]+1e-30),2)
        else:
            ase_val_1 = ase_dict_1[gene]
            ase_val_2 = ase_dict_2[gene]
        ##print ase_dict_1[gene]
        ##print ase_val_1
        if math.fabs(ase_val_1) > 5 or math.fabs(ase_val_2) > 5:
            ##print ase_dict_1[gene],ase_dict_2[gene], ase_val_1,ase_val_2
            continue

	    ase_list_1.append(ase_val_1)
	    ase_list_2.append(ase_val_2)
        ##if gene not in positions:
            ##continue
        ##chrome,start,stop = positions[gene]
        template = '%s\t%s\t%s\t%s\n'
        if ase_val_1 > 0:
            if ase_val_2 >0:
                ##print gene
                ##a.write(template % (chrome, str(start),str(stop),names[gene]))
                up_up += 1
            else:
                ##b.write(template % (chrome, str(start),str(stop),names[gene]))
                up_down +=1
        else:
            if ase_val_2 > 0:
                ##b.write(template % (chrome, str(start),str(stop),names[gene]))
                down_up +=1
            else:
                ##print gene
                ##a.write(template % (chrome, str(start),str(stop),names[gene]))
                down_down +=1
    print up_up,up_down,down_up,down_down

    ##a.flush()
    ##a.close()
    ##b.flush()
    ##b.close()
    ##plt.show()
    plt.figure(figsize=(15,15))
    plt.scatter(ase_list_1,ase_list_2,alpha=.5,marker='o')
    ##plt.plot(xs,eigenVectors[1,0]/eigenVectors[0,0]*xs+ase_list_2.mean(),ls='dashed',color='red')
    plotPCA(ase_list_1,ase_list_2)
    r,pval =  stats.pearsonr(ase_list_1,ase_list_2)
    ##plt.xlabel('PAXB (FW) Allelic Bias\nr=%f' % r,fontsize=22)
    ##plt.ylabel('CERC (FW) Allelic Bias',fontsize=22)
    ##plt.title('Allelic Bias in PAXB (FW) vs CERC (FW)',fontsize=30)
    plt.ylim([-5,5])
    plt.xlim([-5,5])
    ##plt.tight_layout()
    return r
    ##plt.show()
##ase_concordance('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l')

##takes a file with ASE data, and a file with 
def plot_cis_trans(ASE_list_file,trans_exp_diff_file,sig_ase_file,sig_trans_file):
    

    plt.figure(figsize=(12,12))
    ##in the form PAXB/(PAXB+RABS)
    cis_dict = loadASE_list(ASE_list_file)
    ##using the PAXB_RABS file, in the form log2(PAXB/RABS)
    trans_dict = loadASE_list(trans_exp_diff_file)

    sig_cis_dict = loadASE_list(sig_ase_file)
    sig_trans_dict = loadASE_list(sig_trans_file)
    

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
        ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
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
        ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
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
        ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
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
        ase_val = math.log((ase_val+1e-30)/(1-ase_val+1e-30),2)
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


##takes x and y vals
##plots a red dashed PCA1 line fo the given data
def plotPCA(xvals,yvals):
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

    
    xs = np.arange(-5,5,.1)
    plt.plot(xs,eigenVectors[1,0]/eigenVectors[0,0]*xs+xvals.mean(),ls='dashed',color='red')
    m,b = np.polyfit(xvals,yvals,1)

    ##plt.plot(xs,m*xs+b,color='green',linestyle='dashed')
 

##input: a ASE.l file
##output: the dictionary form of the file
def writeASE_names(ase_file,names_file):

    names = module.ENSidDict()
    f = open(ase_file,'r')
    w = open(names_file,'w')

    for line in f:
        line = line.strip()
        split = line.split('\t')
        ensid = split[0]
        gname = names[ensid]
        ase = split[1]
        w.write(gname + '\t' + ase + '\n')

##takes two ASE files (with ASE = p(FW allele)
##computes the predicted imbalence between the two freshwater alleles 
def recalc_ASE(ase_file_1,ase_file_2,outfile):
    ase_1 = loadASE_list(ase_file_1)
    ase_2 = loadASE_list(ase_file_2)

    w = open(outfile,'w')

    ##keys shared among both 
    shared = set(ase_1.keys()).intersection(ase_2.keys())


    for gene in shared:
        ase_val_1 = ase_1[gene]
        ase_val_2 = ase_2[gene]

        r_1 = ase_val_1/(1-ase_val_1)
        r_2 = ase_val_2/(1-ase_val_2)

        ##recalc_ase_val = ase_val_1/( ase_val_1 + ( (1-ase_val_1)/(1-ase_val_2)) * ase_val_2)
        recalc_ase_val = r_1/(r_1 + r_2)
        ##print ase_val_1,ase_val_2,recalc_ase_val
        w.write(gene + '\t' + str(recalc_ase_val) + '\n')

        
##builds a dictionary with the mean values of each ensid at each index sotred with ensid as a key
##indices start at 0, but 0 is the 1st column
def load_cuffnorm_indices(infile,indices):

    f = open(infile,'r')

    toret = {}
    
    f.readline()

    for line in f:
        line = line.strip()
        split = line.split('\t')
        ensid = split[0]
        fpkm_vals = map(lambda x: float(x), split[1:])
        total_fpkm = 0.
        for i in indices:
            total_fpkm += fpkm_vals[i]
        if total_fpkm < .0001:
            continue
        toret[ensid] = total_fpkm/len(indices)

    return toret

##takes a parental expDiff file (in log2 difference) and a hybrid ASE (in ASB (0,1)) file, and a gene exp file
##plots the relative cis and trans contributions for each gene
def cis_vs_trans(parental_file,hybrid_file,exp_file,indices):


    hybrid = loadASE_list(hybrid_file)
    parental = loadASE_list(parental_file)

    parental_exp = load_cuffnorm_indices(exp_file,indices)

	

    ##convert hybrid ase to log2,remove outliers

    ensids = hybrid.keys()
    for ensid in ensids:
        ase_val = hybrid[ensid]
        if ase_val > .99 or ase_val < .01:
            hybrid.pop(ensid)
        else:
            hybrid[ensid] = math.log((ase_val/(1-ase_val)),2)


    shared_ids = set(hybrid.keys()).intersection(set(parental.keys())).intersection(set( parental_exp.keys()))

    ##stores abs value of cis and trans exp changes
    cis_vals = []
    trans_vals = []
    percent_cis_vals = []
    exp_vals = []

    for ensid in shared_ids:
        cis_trans = parental[ensid]
        cis = hybrid[ensid]

        trans = cis_trans - cis

        percent_cis = math.fabs(cis)/(math.fabs(cis) + math.fabs(trans))

        exp = parental_exp[ensid]


        cis_vals.append(math.fabs(cis))
        trans_vals.append(math.fabs(trans))
        percent_cis_vals.append(percent_cis)
        exp_vals.append(math.log(exp,10))

    cis_vals = np.array(cis_vals)
    trans_vals = np.array(trans_vals)
    percent_cis_vals = np.array(percent_cis_vals)
    exp_vals = np.array(exp_vals)

    ##100 must be divisable by BINS
    BINS = 10

    ##add a leftmost bin
    bin_range = [[] for i in range(BINS)]
    
    cur_percent = 0
    ##min leftmost bin
    for i in range(BINS):


        cur_percent += 100/BINS
        ##add a small number so that the max_val is less than the rightmost bin max 
        bin_range[i] = np.percentile(exp_vals,cur_percent) + 1e-10
    print bin_range

    binned_indices = np.digitize(exp_vals,bin_range)
    ##print binned_indices

    binned_cis_vals = [[] for i in range(BINS)]

    for bin_index,cis_percent in zip(binned_indices,percent_cis_vals):
        binned_cis_vals[bin_index].append(cis_percent)

    print len(binned_cis_vals[0]),len(binned_cis_vals[1])

    plt.figure()
    plt.boxplot(binned_cis_vals,notch=True)
        
    print cis_vals.mean(),trans_vals.mean(),percent_cis_vals.mean()
    plt.figure()
    plotSmoothedHist(cis_vals,'cis',x_range=[0,6])
    plotSmoothedHist(trans_vals,'trans',x_range=[0,6])
    plt.legend()
    
    plt.figure()
    plotSmoothedHist(percent_cis_vals,'percent_cis',x_range=[0,1])
    plt.figure()
    plt.scatter(exp_vals,percent_cis_vals,marker='.')
    
    xs = np.arange(-5,5,.1)
    ##plt.plot(xs,eigenVectors[1,0]/eigenVectors[0,0]*xs+xvals.mean(),ls='dashed',color='red')
    m,b = np.polyfit(exp_vals,percent_cis_vals,1)

    plt.plot(xs,m*xs+b,color='green',linestyle='dashed')
    plt.xlim([-1.5,4.5])
    print stats.pearsonr(exp_vals,percent_cis_vals)
    plt.show()
    
##cis_vs_trans('PAXB_RABS_expDiff.l','loose_PAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[3,4,5,6,7,8])
##cis_vs_trans('old_expDiff/PAXB_RABS_expDiff.l','loose_PAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[3,4,5,6,7,8])

##takes a Parental_fc file, and a hybrid expression file (hybrid not in log2)
##also takes an expression file in fpkm, and the indices used within
##writes the calculated trans effects to outfile
def calc_trans(parental_file,hybrid_file,exp_file,exp_ind,outfile,CUTOFF=.1):


    w = open(outfile,'w')

    hybrid = loadASE_list(hybrid_file)
    parental = loadASE_list(parental_file)
    parental_exp = load_cuffnorm_indices(exp_file,exp_ind)

    ensids = hybrid.keys()
    for ensid in ensids:
        ase_val = hybrid[ensid]
        if ase_val > .99 or ase_val < .01:
            hybrid.pop(ensid)
        else:
            hybrid[ensid] = math.log((ase_val/(1-ase_val)),2)


    shared_ids = set(hybrid.keys()).intersection(set(parental.keys())).intersection(set( parental_exp.keys()))

    ##stores abs value of cis and trans exp changes
    cis_vals = []
    trans_vals = []
    percent_cis_vals = []

    for ensid in shared_ids:

        exp = parental_exp[ensid]
        if exp < CUTOFF:
            continue

        cis_trans = parental[ensid]
        cis = hybrid[ensid]

        trans = cis_trans - cis

        percent_cis = math.fabs(cis)/(math.fabs(cis) + math.fabs(trans))


        cis_vals.append(math.fabs(cis))
        trans_vals.append(math.fabs(trans))
        percent_cis_vals.append(percent_cis)

        w.write(ensid + '\t' + str(trans) + '\n')
    w.flush()
    w.close()



##takes a list of differential expression files in lsit form (ENSID\t##)
##takes a list of tuples [(hybrid,parental,population)]
def plot_percent_cis(diff_files,exp_file,indices):

    parental_exp = load_cuffnorm_indices(exp_file,indices)
    plt.figure()

    for hybrid_file,parental_file,population in diff_files:

        hybrid = loadASE_list(hybrid_file)
        parental = loadASE_list(parental_file)
        ##convert hybrid ase to log2,remove outliers

        ensids = hybrid.keys()
        for ensid in ensids:
            ase_val = hybrid[ensid]
            if ase_val > .99 or ase_val < .01:
                hybrid.pop(ensid)
            else:
                hybrid[ensid] = math.log((ase_val/(1-ase_val)),2)


        shared_ids = set(hybrid.keys()).intersection(set(parental.keys())).intersection(set( parental_exp.keys()))

        ##stores abs value of cis and trans exp changes
        cis_vals = []
        trans_vals = []
        percent_cis_vals = []
        exp_vals = []

        for ensid in shared_ids:
            cis_trans = parental[ensid]
            cis = hybrid[ensid]

            trans = cis_trans - cis

            percent_cis = math.fabs(cis)/(math.fabs(cis) + math.fabs(trans))

            exp = parental_exp[ensid]

            if exp < .1:
                continue


            cis_vals.append(math.fabs(cis))
            trans_vals.append(math.fabs(trans))
            percent_cis_vals.append(percent_cis)
            exp_vals.append(math.log(exp,10))

        cis_vals = np.array(cis_vals)
        trans_vals = np.array(trans_vals)
        percent_cis_vals = np.array(percent_cis_vals)
        exp_vals = np.array(exp_vals)



        plotSmoothedHist(percent_cis_vals,population,x_range=[0,1])

##takes a list of tuples in the form [(fpath,name,Bool (T/F log vals use in fpath),color]
def plot_density(flist):

    plt.figure()
    
    
    for fpath,name,logvals,color in flist:
    
        ase_dict = {}

        f = open(fpath,'r')


        for line in f:
            line = line.strip()
            split = line.split('\t')
            ensid = split[0]
            ase_val = float(split[1])
            if not logvals:
                ase_val = math.log((ase_val+1e-20)/((1-ase_val)+1e-20),2)
                ase_dict[ensid] = ase_val


        all_ase_values = ase_dict.values()
        all_ase_values = np.array(all_ase_values)

        print all_ase_values.mean()
        plotSmoothedHist(all_ase_values,name,x_range=[-4,4],cur_color=color)


    plt.xlabel('log2(F/M)',fontsize='large')
    plt.ylabel('Density',fontsize='large')
    plt.ylim([0,1])
    xs = [0] * len(range(0,8))
    ys = range(0,8)
    plt.plot(xs,ys,ls='dashed',color='red')
    plt.legend()
    ##plt.show()


##plots the concordance between two .l files
##color-coded density options
def density_concord(ase_file_1,ase_file_2,log_vals=False,psize=5):
    ase_dict_1 = loadASE_list(ase_file_1)
    ase_dict_2 = loadASE_list(ase_file_2)
    ##get a list of genes shared between the two

    genes = set(ase_dict_1.keys()).intersection(ase_dict_2.keys())

    ##a = open('/data/James/AlleleSpecificExpression/VTP/UUDD_CxP_ASE.bed','w')
    ##b = open('/data/James/AlleleSpecificExpression/VTP/DUUD_CxP_ASE.bed','w')

    positions = module.ensidStartStop()
    names = module.ENSidDict()

    ase_list_1 = []
    ase_list_2 = []
    up_up = 0
    down_down = 0
    up_down = 0
    down_up = 0
    
    for gene in genes:

        ##log2 the ratio value
        ##print ase_dict_1[gene],ase_dict_2[gene]
        if not log_vals:
            ase_val_1 = math.log((ase_dict_1[gene] + 1e-30)/(1-ase_dict_1[gene]+1e-30),2)
            ase_val_2 = math.log((ase_dict_2[gene] + 1e-30)/(1-ase_dict_2[gene]+1e-30),2)
        else:
            ase_val_1 = ase_dict_1[gene]
            ase_val_2 = ase_dict_2[gene]

        if math.fabs(ase_val_1) > 5 or math.fabs(ase_val_2) > 5:
            ##print ase_dict_1[gene],ase_dict_2[gene], ase_val_1,ase_val_2
            continue

        ase_list_1.append(ase_val_1)
        ase_list_2.append(ase_val_2)

        template = '%s\t%s\t%s\t%s\n'
        if ase_val_1 > 0:
            if ase_val_2 >0:
                ##print gene
                ##a.write(template % (chrome, str(start),str(stop),names[gene]))
                up_up += 1
            else:
                ##b.write(template % (chrome, str(start),str(stop),names[gene]))
                up_down +=1
        else:
            if ase_val_2 > 0:
                ##b.write(template % (chrome, str(start),str(stop),names[gene]))
                down_up +=1
            else:
                ##print gene
                ##a.write(template % (chrome, str(start),str(stop),names[gene]))
                down_down +=1
    print up_up,up_down,down_up,down_down

    plt.figure(figsize=(15,15))
    ##plt.scatter(ase_list_1,ase_list_2)
    xy = np.vstack([ase_list_1,ase_list_2])
    z = gaussian_kde(xy)(xy)
    plt.scatter(ase_list_1,ase_list_2, c=z, s=psize, edgecolor='')
    plotPCA(ase_list_1,ase_list_2)
    r,pval =  stats.pearsonr(ase_list_1,ase_list_2)
    plt.ylim([-5,5])
    plt.xlim([-5,5])
    ##plt.tight_layout()
    return r
    ##plt.show()

r = density_concord('PAXB_RABS_expDiff.l','CERC_RABS_expDiff.l',log_vals=True)
plt.xlabel('log2 (PAXB / RABS) expression\nr=%f'%r,size='x-large')
plt.ylabel('log2 (CERC / RABS) expression',size='x-large')
plt.show()
    
'''
plot_percent_cis([('loose_PAXB_VTP_ASE.l','PAXB_RABS_expDiff.l','PAXBxRABS'),('CERCxPAXB_VTP_ASE.l','CERC_PAXB_expDiff.l','CERCxPAXB'),('CERC_VTP_ASE.l','CERC_RABS_expDiff.l','CERCxRABS')],'/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,3,4,5,6,7,8])
plt.legend()
plt.xlabel('Percent Cis',size='x-large')
plt.ylabel('Density',size='x-large')
plt.savefig('/home/james/Dropbox/Miller/figures/Convergent_Evolution/Fig_4/Percent_cis_histogram_combined.png')
plt.show()
'''

'''

calc_trans('PAXB_RABS_expDiff.l','loose_PAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[3,4,5,6,7,8],'PAXB_RABS_transDiff.l')
calc_trans('CERC_RABS_expDiff.l','CERC_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,3,4,5],'CERC_RABS_transDiff.l')
calc_trans('CERC_PAXB_expDiff.l','CERCxPAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,6,7,8],'CERC_PAXB_transDiff.l')

print ase_concordance('loose_PAXB_VTP_ASE.l','CERC_VTP_ASE.l')
print ase_concordance('PAXB_RABS_transDiff.l','CERC_RABS_transDiff.l',log_vals=True)
print ase_concordance('PAXB_RABS_expDiff.l','CERC_RABS_expDiff.l',log_vals=True)

print ase_concordance('PAXB_RABS_transDiff.l','CERC_RABS_transDiff.l',log_vals=True)
print ase_concordance('PAXB_RABS_transDiff.l','CERC_PAXB_transDiff.l',log_vals=True)
print ase_concordance('CERC_PAXB_transDiff.l','CERC_RABS_transDiff.l',log_vals=True)

calc_trans('PAXB_RABS_sigDiff.l','loose_PAXB_VTP_sigASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[3,4,5,6,7,8],'PAXB_RABS_sigtransDiff.l')
calc_trans('CERC_RABS_sigDiff.l','CERC_VTP_sigASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,3,4,5],'CERC_RABS_sigtransDiff.l')

print ase_concordance('loose_PAXB_VTP_sigASE.l','CERC_VTP_sigASE.l')
print ase_concordance('PAXB_RABS_sigtransDiff.l','CERC_RABS_sigtransDiff.l',log_vals=True)
print ase_concordance('PAXB_RABS_sigDiff.l','CERC_RABS_sigDiff.l',log_vals=True)

ase_concordance('PAXB_RABS_sigDiff.l','CERC_PAXB_sigDiff.l',log_vals=True)
ase_concordance('CERC_PAXB_sigDiff.l','CERC_RABS_sigDiff.l',log_vals=True)

print ase_concordance('loose_PAXB_VTP_ASE.l','CERC_VTP_ASE.l')
print ase_concordance('CERC_VTP_ASE.l','CERCxPAXB_VTP_ASE.l')
print ase_concordance('loose_PAXB_VTP_ASE.l','CERCxPAXB_VTP_ASE.l')

print ase_concordance('loose_PAXB_VTP_sigASE.l','CERC_VTP_sigASE.l')
print ase_concordance('CERC_VTP_sigASE.l','CERCxPAXB_VTP_sigASE.l')
print ase_concordance('loose_PAXB_VTP_sigASE.l','CERCxPAXB_VTP_sigASE.l')
'''
'''
    
cis_vs_trans('PAXB_RABS_expDiff.l','loose_PAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[3,4,5,6,7,8])
##0.345843744095 0.641302713558 0.397060711531
##(0.034637228373993918, 0.0074641598991580559)
cis_vs_trans('CERC_RABS_expDiff.l','CERC_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,3,4,5])
##0.389294047974 0.610795582681 0.436915938611
##(0.17348594754195201, 1.3818962939391299e-15)
cis_vs_trans('CERC_PAXB_expDiff.l','CERCxPAXB_VTP_ASE.l','/home/james/Dropbox/Miller/data/RNA_seq_exp/CERC_RABS_PAXB_ids.fpkm_table',[0,1,2,6,7,8])
##0.494933188825 0.778767303571 0.431127140648
##(0.081887997053347095, 0.00029478767015112624)
'''
'''
recalc_ASE('CERC_VTP_sigASE.l','loose_PAXB_VTP_sigASE.l','CERCxPAXB_VTP_sigASE.l')
recalc_ASE('CERC_VTP_ASE.l','loose_PAXB_VTP_ASE.l','CERCxPAXB_VTP_ASE.l')
plot_cis_trans('CERCxPAXB_VTP_ASE.l','CERC_PAXB_expDiff.l','CERCxPAXB_VTP_sigASE.l','CERC_PAXB_sigDiff.l')
plt.xlabel('log2( CERC(FW)_p / PAXB(M)_p )\n"trans"',fontsize=26)
plt.ylabel('log2( CERC(FW)_F1 / PAXB(M)_F1 (est.) )\n"cis"',fontsize=26)
plt.title('PAXB(FW) x CERC (FW) Changes in Gene Expression',fontsize=30)
plt.savefig('/home/james/Dropbox/Miller/figures/Convergent_Evolution/CERCxPAXB_cis_trans.png')
##print I,II,III,IV
##575 448 597 474

r = ase_concordance('loose_PAXB_VTP_sigASE.l','CERC_VTP_sigASE.l')
plt.xlabel('PAXB (FW) VTP Allelic Bias\nr=%f' % r,fontsize=22)
plt.ylabel('CERC(FW) VTP Allelic Bias',fontsize=22)
plt.title('Allelic Bias in PAXB (FW) VTP vs CERC (FW) VTP',fontsize=22)
plt.savefig('/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE_concord_P-VTP_C-VTP.png')
plt.clf()
##180 138 127 190
r = ase_concordance('loose_PAXB_BS_sigASE.l','CERC_VTP_sigASE.l')
plt.xlabel('PAXB (FW) BS Allelic Bias\nr=%f' % r,fontsize=22)
plt.ylabel('CERC (FW) VTP Allelic Bias',fontsize=22)
plt.title('Allelic Bias in PAXB (FW) BS vs CERC (FW) VTP',fontsize=30)
plt.savefig('/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE_concord_P-BS_C-VTP.png')
plt.show()
plt.clf()
##182 147 138 194

r = ase_concordance('loose_PAXB_BS_sigASE.l','loose_PAXB_VTP_sigASE.l')
plt.xlabel('PAXB (FW) BS Allelic Bias\nr=%f' % r,fontsize=22)
plt.ylabel('PAXB (FW) VTP Allelic Bias',fontsize=22)
plt.title('Allelic Bias in PAXB (FW) BS vs PAXB (FW) VTP',fontsize=30)
plt.savefig('/home/james/Dropbox/Miller/figures/Convergent_Evolution/ASE_concord_P-BS_P-VTP.png')
plt.show()
plt.clf()
##689 82 87 690
writeASE_names('CERC_VTP_ASE.l','CERC_VTP_ASE_names.l')

ASE_density('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_ASE.l')
ASE_density('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l')

ASE_density('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_ASE.l')
ASE_density('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l')

ase_concordance('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/PAXB_BS/loose_PAXB_BS_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/PAXB_BS/loose_PAXB_BS_sigASE.l','/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l')

'''
'''
plot_cis_trans('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_ASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/PAXB_RABS_expDiff.l','/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/PAXB_RABS_sigDiff.l')
##plot_cis_trans('/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_ASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_expDiff.l','/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_sigDiff.l')
##plot_cis_trans('/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_expDiff.l','/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_sigDiff.l')
plot_cis_trans('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_expDiff.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_sigDiff.l')

gtfToGnames('PAXB_nVTP_candidates.gtf')
gtfToGnames('CERC_nVTP_candidates.gtf')

findQTLcandidates('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_diff/PAXB_RABS_sigDiff.l','/data/James/QTL_candidates/PAXB_candidates_sig.gtf')
findQTLcandidates('/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_PAXB_diff/CERC_RABS_sigDiff.l','/data/James/QTL_candidates/CERC_candidates_sig.gtf')
'''
'''
ase_concordance('/data/James/AlleleSpecificExpression/PAXB_BS/loose_PAXB_BS_sigASE.l','/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/PAXB_BS/loose_PAXB_BS_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_sigASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/VTP/loose_PAXB_VTP_sigASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_sigASE.l')
ase_concordance('/data/James/AlleleSpecificExpression/VTP/CERC/loose_CERC_VTP_ASE.l','/data/James/AlleleSpecificExpression/VTP/CERC/CERC_VTP_ASE.l')
'''


##allPathways('loose_VTP_ASE.l','test.out')
##ASE_density('/data/James/AlleleSpecificExpression/PAXB_BS/BS_ASE_all.l')
##ASEerror_ase('/data/James/AlleleSpecificExpression/simVTP/simVTP_ASE_all.l','/data/James/AlleleSpecificExpression/simVTP/JCH003A_sim_ASE.l')
##aseErr_expression('/data/James/AlleleSpecificExpression/VTP/JCH003A_forced_collapsed_cuff/genes.fpkm_tracking','/data/James/AlleleSpecificExpression/simVTP/simVTP_ASE_all.l','/data/James/AlleleSpecificExpression/simVTP/JCH003A_sim_ASE.l')
##ase_expression('/data/James/AlleleSpecificExpression/VTP/JCH003A_forced_collapsed_cuff/genes.fpkm_tracking','/data/James/AlleleSpecificExpression/VTP/VTP_ASE_all.l')


##ASEerror_density('/data/James/AlleleSpecificExpression/simVTP/VTP_sim_noASE_all.l','/data/James/AlleleSpecificExpression/simVTP/JCH003A_sim_noASE_ASE.l')
##ASEerror_density('/data/James/AlleleSpecificExpression/simVTP/simVTP_ASE_all.l','/data/James/AlleleSpecificExpression/simVTP/JCH003A_sim_ASE.l')

##ase_concordance('/data/James/AlleleSpecificExpression/VTP/VTP_ASE_all.l','/data/James/AlleleSpecificExpression/PAXB_BS/BS_ASE_all.l')
##ASEerror_density('/data/James/AlleleSpecificExpression/VTP/VTP_ASE_all.l','/data/James/AlleleSpecificExpression/PAXB_BS/BS_ASE_all.l')

##plot_cis_trans('/data/James/AlleleSpecificExpression/VTP/VTP_ASE_all.l','/data/James/Parental_RNA_seq/CERC_RABS_diff/PAXB_RABS_expDiff.l','/data/James/AlleleSpecificExpression/VTP/VTP_sigASE.l','/data/James/Parental_RNA_seq/CERC_RABS_diff/PAXB_RABS_sigDiff.l')

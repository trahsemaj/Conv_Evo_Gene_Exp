#! /usr/bin/env python

import sys,getopt
import subprocess as sp
import module


'''
personalTranscriptome2.py takes a phased (or partially phased) VCF, with an optional ref genome / ref transcriptome gtf

It outputs a phased transcriptome, with two copies of each transcript (one per haplotype)
It will only write perfectly phased variants, unphased variants will appear as refseq

INPUT: [options] Phased VCF
Phased vcf is a phased vcf file
refTrascriptome.fa is the .fa file of the reference
POS is the position (starting a 0) of the child in the phased haplotype.  Usually 2
hapPairs is the outfile for haplotype-file, needed for express
OUTPUT: personalTransciptome.fa to STDOUT

will use a gtf by default, but given the -b optin will use a BED instead

'''

##loads a phased vcf into a dict in the form vcfdict[chrpos] = [allele0,allele1,len(ref) the length of the reference allele])
def loadPhasedVCF(VCF,POS):
    
    vcf = open(VCF,'r')
    HAPLOVARS = {}
    for line in vcf:
        
        if line.startswith('#'):
            continue
        
        split = line.split()
        if len(split) < 9:
            continue
        genotypes = split[9+POS][:3].split('|')
        
        
        chr = split[0]
        ##vcfs are 1-based
        pos = str(int(split[1])-1)
        ##alleles = ref,alt (selected by 0,1 genotypes), only take the first multiallele
        alleles = [split[3],split[4].split(',')[0]]
        
        ##if both alleles are ref, don't return anything
        if int(genotypes[0]) == 0 and int(genotypes[1]) == 0:
            continue
        
        phasedAlleles = [alleles[int(genotypes[0])],alleles[int(genotypes[1])],len(split[3])]
        
        HAPLOVARS[chr+':'+pos] = phasedAlleles
        
    return HAPLOVARS


def writeDiploid(genome_file, gtf_file, vcf_index, vcf_file):
    ##load the reference
    REF = module.loadRef(genome_file)

    
    ##load phased variants from the vcf file
    variants = loadPhasedVCF(vcf_file,vcf_index)
    
    ##loadGTF, print alleles gene by gene
    f = open(gtf_file)
    
    haplist = []
    
    fasta_template = '>{0}_1\n{1}\n>{0}_2\n{2}'
    haptemp = '{0}_1,{0}_2'
    chr = ''
    oldId = ''
    id = ''
    seq0 = ''
    seq1 = ''
    start = -1
    end = -1
    strand = ''
    i = 0
    firststart = 0
    for line in f:
        split = line.split('\t')


        ##only look at 'exons'
        if split[2] != 'exon':
            continue
    
        ##initialize
        if start == -1:
            start = int(split[3])
            end =  int(split[4])
            chr = split[0]
            strand = split[6]
            oldId = split[8].split('"')[3]
            firststart = start

        id = split[8].split('"')[3]
    
        ##new transcript
        if id != oldId:
            
            
            i+=1
            ##print >> sys.stderr, start,end
            ##print the sequence if it isn't all masked
            if seq0.count('N') < len(seq0):
                print fasta_template.format(oldId,seq0,seq1)
                haplist.append(haptemp.format(oldId))
            oldId = id
            start = int(split[3])
            end = int(split[4])
            chr = split[0]
            strand = split[6]
            seq0 = ''
            seq1 = ''
            firststart = start
        
        start = int(split[3])    
        end = int(split[4])
        strand = split[6]
    
        ##gtfs are 1-indexed, inclusive on both ends
        i = start-1
        while i < end:
            chrpos = chr + ':' + str(i)
            ##if a variant position
            if chrpos in variants:
                var =  variants[chrpos]
                ##print var
                seq0 += var[0]
                seq1 += var[1]
                i += variants[chrpos][2]
            else:
                seq0 += REF[chr][i]
                seq1 += REF[chr][i]
                i += 1
        ##seq = seq + REF[chr][start-1:end]
        

    ##print the last stored sequence
    if seq0.count('N') < len(seq0):
        print fasta_template.format(oldId,seq0,seq1)
        haplist.append(haptemp.format(oldId))

    return haplist



##will generate the RNA-seq experiment
def main(args):
    
    
    REF_GENOME_FILE = '/home/ctmiller/ReferenceGenomes/StickleUnmasked/stickleUnmasked.fa'
    ##REF_GENOME_FILE ='/home/ctmiller/ReferenceGenomes/StickleMasked/stickleMasked.fa'
    ##REF_GTF_FILE = '/data/James/BED/ensGene.gtf'
    ##REF_GTF_FILE = '/data/James/AlleleSpecificExpression/simulationBS/chrXXI/JCH001D_XXI.gtf'
    ##REF_GTF_FILE = '/data/James/AlleleSpecificExpression/PAXB_BS/JCH001D_collapsed_forced_cuff/transcripts.gtf'
    REF_GTF_FILE = '/data/James/BED/ensGene_collapsed.gtf'
    HAP_FILE = ''
    VCF_INDEX = 0
    VCF_FILE = ''
    
    optlist,args = getopt.getopt(args[1:],'r:g:h:p:')

    for a in optlist:
        if a[0] == '-r':
            REF_GENOME_FILE = a[1]
        if a[0] == '-g':
            REF_GTF_FILE = a[1]
        if a[0] == '-h':
            HAP_FILE = a[1]
        if a[0] == '-p':
            VCF_INDEX = int(a[1])


    VCF_FILE = args[0]
    
    ##writes the diploid genome given the phased variants to stdout, stores the pairs of haplotypes in haplist
    haplist = writeDiploid(REF_GENOME_FILE, REF_GTF_FILE, VCF_INDEX, VCF_FILE)
    
    if HAP_FILE != '':
        w = open(HAP_FILE,'w')
        w.write('\n'.join(haplist))
        w.flush()
        w.close()
    ##write('\n'.join(haplist))
    
    


















if __name__ == '__main__':
    sys.exit(main(sys.argv))
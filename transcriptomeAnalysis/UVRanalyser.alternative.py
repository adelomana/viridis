###
### This script compares the UVR-response genes in HC vs LC.
### It compares if there is consistent (exp AND sta) differential expression in epoch 0 vs 1, and this set of UVR-responding genes is different in HC vs LC.
###

import sys
import library

def consistencyFinder(label):

    '''
    This function reads the output of cuffdiff and finds sets of genes in common with other cuffdiff runs.
    '''

    dir=cuffdiffDir+'uvr.alternative.{}.exp.AM.epoch0.vs.epoch1'.format(label)
    a=cuffdiffReader(dir)

    dir=cuffdiffDir+'uvr.alternative.{}.exp.PM.epoch0.vs.epoch1'.format(label)
    b=cuffdiffReader(dir)

    dir=cuffdiffDir+'uvr.alternative.{}.sta.AM.epoch0.vs.epoch1'.format(label)
    c=cuffdiffReader(dir)

    dir=cuffdiffDir+'uvr.alternative.{}.sta.PM.epoch0.vs.epoch1'.format(label)
    d=cuffdiffReader(dir)

    e1=list(set(a) & set(b))
    e2=list(set(c) & set(d))
    e=list(set(e1) & set(e2))

    print(len(a),len(b),len(c),len(d),len(e))

    return a,b,c,d,e

def cuffdiffReader(dir):

    '''
    This function reads cuffdiff output and selects DETs.
    '''

    DETs=[]
    with open(dir+'/isoform_exp.diff','r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')

            transcriptName=vector[0]
            expressionA=vector[7]
            expressionB=vector[8]
            foldChange=vector[9]
            pValue=vector[11]
            significance=vector[13].replace('\n','')
                
            if significance == 'yes':
                DETs.append(transcriptName)

    return DETs

def hypothesisTestingRunner(label):

    '''
    This function calls cuffdiff to run sample comparisons.
    '''

    # f.1.
    samplesA=[]; samplesB=[]
    for sampleID in metadata.keys():
        if metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'exp':
            samplesA.append(sampleID)
        elif metadata[sampleID]['co2'] == 1000 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'exp':
            samplesB.append(sampleID)

    print(samplesA,samplesB)
    bamFilesA,bamFilesB=library.samples2bamfiles(samplesA,samplesB,bamFilesDir)
    library.cuffdiffCaller(bamFilesA,bamFilesB,'uvr.alternative.epoch{}.exp.AM.LC.vs.HC'.format(label),cuffdiffDir,gtfFile,fastaFile,numberOfThreads)

    # f.2.
    samplesA=[]; samplesB=[]
    for sampleID in metadata.keys():
        if metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'PM' and metadata[sampleID]['growth'] == 'exp':
            samplesA.append(sampleID)
        elif metadata[sampleID]['co2'] == 1000 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'PM' and metadata[sampleID]['growth'] == 'exp':
            samplesB.append(sampleID)

    print(samplesA,samplesB)
    bamFilesA,bamFilesB=library.samples2bamfiles(samplesA,samplesB,bamFilesDir)
    library.cuffdiffCaller(bamFilesA,bamFilesB,'uvr.alternative.epoch{}.exp.PM.LC.vs.HC'.format(label),cuffdiffDir,gtfFile,fastaFile,numberOfThreads)

    # f.3.
    samplesA=[]; samplesB=[]
    for sampleID in metadata.keys():
        if metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'sta':
            samplesA.append(sampleID)
        elif metadata[sampleID]['co2'] == 1000 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'sta':
            samplesB.append(sampleID)

    print(samplesA,samplesB)
    bamFilesA,bamFilesB=library.samples2bamfiles(samplesA,samplesB,bamFilesDir)
    library.cuffdiffCaller(bamFilesA,bamFilesB,'uvr.alternative.epoch{}.sta.AM.LC.vs.HC'.format(label),cuffdiffDir,gtfFile,fastaFile,numberOfThreads)

    # f.4.
    samplesA=[]; samplesB=[]
    for sampleID in metadata.keys():
        if metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'PM' and metadata[sampleID]['growth'] == 'sta':
            samplesA.append(sampleID)
        elif metadata[sampleID]['co2'] == 1000 and metadata[sampleID]['epoch'] == label and metadata[sampleID]['diurnal'] == 'PM' and metadata[sampleID]['growth'] == 'sta':
            samplesB.append(sampleID)

    print(samplesA,samplesB)
    bamFilesA,bamFilesB=library.samples2bamfiles(samplesA,samplesB,bamFilesDir)
    library.cuffdiffCaller(bamFilesA,bamFilesB,'uvr.alternative.epoch{}.sta.PM.LC.vs.HC'.format(label),cuffdiffDir,gtfFile,fastaFile,numberOfThreads)

    
    return None

def uvrRespondingGenesReader():

    '''
    This function retrieves the names of UVR-responding genes.
    '''

    names=[]
    with open(uvrRespondingGenesFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split(',')
            names.append(vector[1])

    names.sort()

    return names

###
### MAIN
###

# 0. user defined variables
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
uvrRespondingGenesFile='/Volumes/omics4tb/alomana/projects/dtp/data/functionalAnnotation/UVR_Responsive_Genes.csv'
bamFilesDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cuffdiffDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
gtfFile='/Volumes/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
fastaFile='/Volumes/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.dna.genome.fa'
numberOfThreads=16 

# 1. read data
print('reading data...')
# 1.1. reading metadata
metadata=library.metadataReader(metaDataFile)

# 1.2. reading expression
expression=library.expressionReader(expressionFile)
sortedGeneNames=list(expression.keys())
sortedGeneNames.sort()

# 1.3. read UVR-responding genes
uvrGenes=uvrRespondingGenesReader()

# 2. run hypotheses testing
print('running hypotheses testing...')
hypothesisTestingRunner(0)
hypothesisTestingRunner(1)

sys.exit()

# 3. make the intersect
print('defining consistency and intersect...')

# 3.1. define consistency in exp and sta for each condition
LCa,LCb,LCc,LCd,LCe=consistencyFinder(300)
HCa,HCb,HCc,HCd,HCd=consistencyFinder(1000)

# 3.2.1. for all genes

# 3.2.2. for UVR genes


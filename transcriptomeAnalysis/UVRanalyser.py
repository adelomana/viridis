###
### This script compares the UVR-response genes in HC vs LC.
### It compares if there is consistent (exp AND sta) differential expression in epoch 0 vs 1, and this set of UVR-responding genes is different in HC vs LC.
###

import sys
import library

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
expressionFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
uvrRespondingGenesFile='/proj/omics4tb/alomana/projects/dtp/data/functionalAnnotation/UVR_Responsive_Genes.csv'
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cuffdiffDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
fastaFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.dna.genome.fa'
numberOfThreads=4 

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

# 2.1. detect DETs between epoch 0 and 1 in LC

# 2.1.1. define DETs between exp
samplesA=[]; samplesB=[]
for sampleID in metadata.keys():
    if metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == 0 and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'exp':
        samplesA.append(sampleID)
    elif metadata[sampleID]['co2'] == 300 and metadata[sampleID]['epoch'] == 1 and metadata[sampleID]['diurnal'] == 'AM' and metadata[sampleID]['growth'] == 'exp':
        samplesB.append(sampleID)

bamFilesA,bamFilesB=library.samples2bamfiles(samplesA,samplesB,bamFilesDir)
library.cuffdiffCaller(bamFilesA,bamFilesB,'uvr.300.exp.epoch0.vs.epoch1',cuffdiffDir,gtfFile,fastaFile,numberOfThreads)

# 2.1.2. define DETs in sta



# 2.2. detect DETs between epoch 0 and 1 in HC

# 3. make the intersect
print('defining consistency and intersect...')

# 3.1. define consistency in exp and sta for each condition

# 3.2. define comparison for LC and HC


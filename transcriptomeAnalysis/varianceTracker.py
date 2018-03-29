###
### This script tracks the variability of gene expression.
### We look for genes that have increase in variability. The pattern should be 01 and 001 for LC and HC respectively.
###

import sys,numpy
import library

###
### MAIN
###

print('\nwelcome to varianceTracker...\n')

# 0. user defined variables
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'

# 1. read data
print('reading data...')
# 1.1. reading metadata
metadata=library.metadataReader(metaDataFile)

# 1.2. reading expression
expression=library.expressionReader(expressionFile)
sortedGeneNames=list(expression.keys())
sortedGeneNames.sort()

# 2. analysis
print('performing analysis...')

# 2.1. define the ordered set of conditions
samples300=library.sampleOrderer(300,metadata)
samples1000=library.sampleOrderer(1000,metadata)

# 2.2. iterate over genes
print('\t iterating {} genes...'.format(len(sortedGeneNames)))
candidates={}
for geneName in sortedGeneNames:

    # define LC expression
    y300=[]
    for condition in sorted(samples300.keys()):
        y=[]
        for sampleID in samples300[condition]:
            expressionValue=expression[geneName][sampleID]
            y.append(expressionValue)
        y300.append(y)

    # define HC expression
    y1000=[]
    for condition in sorted(samples1000.keys()):
        y=[]
        for sampleID in samples1000[condition]:
            expressionValue=expression[geneName][sampleID]
            y.append(expressionValue)
        y1000.append(y)

    # format expressions
    E300=numpy.array(y300)
    E1000=numpy.array(y1000)

    print(E300)
    print(E1000)
    
    mE300=numpy.mean(E300,axis=1)
    mE1000=numpy.mean(E1000,axis=1)

    print(E300,mE300)
    print(E1000,mE1000)
    
    #sE300=numpy.sds(E300,axis=1); sE1000=numpy.sds(E1000,axis=1)

    #maxE300=numpy.max(mE300); maxE1000=numpy.max(mE1000)

    # filter out genes whose maximum expression is below 10. FPKMs
    if max([maxE300,maxE1000]) > 10:
        print(geneName)
        print(mE300)
        print(mE1000)
        print()
    else:
        print('Excluding {} because low expression...')
        print(mE300)
        print(mE1000)
        print()
    
    
    # filter out genes whose abs(log2 FC) < 1
    
    

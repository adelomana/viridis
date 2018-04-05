###
### This script compares the UVR-response genes in HC vs LC.
### It compares if there is consistent expression difference in epoch 0 vs 1 with respect to HC and LC on UVR-response genes.
###

import sys,numpy
import scipy,scipy.stats
import matplotlib,matplotlib.pyplot
import library

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def histogrammer(theData):

    '''
    this function creates a histogram.
    '''

    n,bins=numpy.histogram(theData,bins=int(numpy.sqrt(len(theData))))

    x=[]
    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        x.append(center)
    x.pop()

    y=[]
    y=numpy.array(n)
    y=list(y/float(sum(y)))

    return x,y

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
epochs=[0,1]
growths=['exp','sta']
diurnals=['AM','PM']

diffUVR=[]; diffnoUVR=[]

for gene in uvrGenes:
    name='gene:Thaps{}'.format(gene)
    print('working with gene {}...'.format(name))

    # 2.1. find trajectory for LC
    orderedExpressionLC=[]
    co2level=300

    for epoch in epochs:
        for growth in growths:
            for diurnal in diurnals:
                
                localExpression=[]
                for sampleID in metadata.keys():
                    if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:
                        value=expression[name][sampleID]
                        localExpression.append(value)
                orderedExpressionLC.append(localExpression)
                
    # 2.2. find trajectory for HC
    orderedExpressionHC=[]
    co2level=1000

    for epoch in epochs:
        for growth in growths:
            for diurnal in diurnals:
                
                localExpression=[]
                for sampleID in metadata.keys():
                    if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:
                        value=expression[name][sampleID]
                        localExpression.append(value)
                orderedExpressionHC.append(localExpression)

    # 2.3. find differences
    LC=numpy.array(orderedExpressionLC)
    HC=numpy.array(orderedExpressionHC)

    mLC=numpy.mean(LC,axis=1)
    mHC=numpy.mean(HC,axis=1)

    log2fcLC=numpy.log2(mLC/mLC[0])
    log2fcHC=numpy.log2(mHC/mLC[0])

    # adding differences
    noUVR=log2fcHC[0:4]-log2fcLC[0:4]
    UVR=log2fcHC[4:]-log2fcLC[4:]

    for i in range(len(noUVR)):
        if numpy.isnan(UVR[i]) == False:
            diffUVR.append(UVR[i])
        if numpy.isnan(noUVR[i]) == False:
            diffnoUVR.append(noUVR[i])

# 3. compute the difference between UVR and no UVR and plot a histogram
statistic,pvalue=scipy.stats.ttest_ind(diffUVR,diffnoUVR)
print(statistic,pvalue)

x,y=histogrammer(diffnoUVR)
matplotlib.pyplot.plot(x,y,'-',color='green',lw=2,label='Stage 1')

x,y=histogrammer(diffUVR)
matplotlib.pyplot.plot(x,y,'-',color='magenta',lw=2,label='Stage 2')

matplotlib.pyplot.xlabel('Carbon response [log$_2$ FC (HC) - log$_2$ FC (LC)]')
matplotlib.pyplot.ylabel('Probability density function')

matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=1,fontsize=18)

matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('histogram.pdf')
matplotlib.pyplot.clf()

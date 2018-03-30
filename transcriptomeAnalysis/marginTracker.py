###
### This script tracks the margin of gene expression for each gene along different conditions.
###

import sys,numpy
import matplotlib,matplotlib.pyplot
import library

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def histogrammer(theData):

    '''
    This function creates a histogram.
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

###
### MAIN
###

print('\nwelcome to marginTracker...\n')

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

margins={}

co2levels=[300,1000]
epochs=[0,1,2]
growths=['exp','sta']
diurnals=['AM','PM']

for co2level in co2levels:
    margins[co2level]=[]
    for epoch in epochs:
        for growth in growths:
            for diurnal in diurnals:
                print(co2level,epoch,growth,diurnal)
                # obtain the replicate sample IDs
                sampleIDs=[]
                for sampleID in metadata.keys():
                    if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:
                        sampleIDs.append(sampleID)
                # fill the margin array
                localMargins=[]
                for geneName in sortedGeneNames:
                    m=[]
                    for sampleID in sampleIDs:
                        value=expression[geneName][sampleID]
                        logValue=numpy.log10(value+1)
                        m.append(logValue)
                    if m != []:
                        margin=max(m)-min(m)
                        # removing low expressed transcripts
                        if min(m) >= 1:
                            localMargins.append(margin)
                if len(localMargins) != 0:
                    # remove 2.5% at each side
                    localMargins.sort()
                    extreme=int(len(localMargins)*0.025)-1
                    trimmedLocalMargins=localMargins[extreme:-extreme]
                    # compute PDF
                    x,y=histogrammer(trimmedLocalMargins)
                    margins[co2level].append([x,y])
    # 3. plotting figures
    for i in range(len(margins[co2level])):
        curve=margins[co2level][i]
        matplotlib.pyplot.plot(curve[0],curve[1],'-',color=matplotlib.cm.tab20(i),label=i)
        # fit a log-normal distribution
        
    # close figure
    figureName='figure.{}.pdf'.format(co2level)
    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=1,fontsize=12)
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()
    

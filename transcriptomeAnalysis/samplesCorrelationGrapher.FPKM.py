### this script computes a pairwise Spearman correlation based on the intersect of genes that have more than tau FPKM.

import os,sys,scipy,matplotlib,numpy,random
from scipy import stats
from matplotlib import pyplot
matplotlib.rcParams['pdf.fonttype']=42 # this cryptical line is necessary for Illustrator compatibility of text saved as pdf
matplotlib.rcParams.update({'xtick.labelsize':20,'ytick.labelsize':20})

def correlationQuantifier(working_co2level,working_epoch,biorepSCCs):

    '''
    this function computes the Spearman correlation coefficient between samples
    it builds a heatmap and a boxplot figure per epoch
    '''

    epochs=[0,1,2]
    growths=['exp','sta']
    diurnals=['AM','PM']
    replicates=['A','B','C']

    orderedSamples=[]
    
    
    for growth in growths:
        for diurnal in diurnals:
            for replicate in replicates:
                for sampleID in metaData.keys():
                    if metaData[sampleID]['co2'] == working_co2level and metaData[sampleID]['epoch'] == working_epoch and metaData[sampleID]['growth'] == growth and metaData[sampleID]['diurnal'] == diurnal and metaData[sampleID]['replicate'] == replicate:
                        orderedSamples.append(sampleID)
                                
    # computing the matrix of correlations
    z=len(orderedSamples)
    M=numpy.ones((z,z))
    conditionNames=[]
    lowDomain=10000
    
    for i in range(z):
        if metaData[orderedSamples[i]]['growth'] == 'exp':
            word1='ear'
        else:
            word1='lat'
        if  metaData[orderedSamples[i]]['diurnal'] == 'AM':
            word2='li'
        else:
            word2='da'
            
        conditionName='{}.{}.{}'.format(word1,word2,metaData[orderedSamples[i]]['replicate'])
        conditionNames.append(conditionName)
        
        for j in range(z):

            if i < j:
                expressionA=[]
                expressionB=[]

                for gene in sortedGenes:
                    a=expression[gene][orderedSamples[i]]
                    b=expression[gene][orderedSamples[j]]

                    if a >= tau and b >= tau:
                        expressionA.append(a)
                        expressionB.append(b)

                # computing correlation
                rho,pval=scipy.stats.spearmanr(expressionA,expressionB)
                M[i,j]=rho

                #if len(expressionA) < lowDomain:
                #    lowDomain=len(expressionA)
                #    print(lowDomain)

                # selecting values for biorepSCCs
                growthA=metaData[orderedSamples[i]]['growth']
                growthB=metaData[orderedSamples[j]]['growth']
                diurnalA=metaData[orderedSamples[i]]['diurnal']
                diurnalB=metaData[orderedSamples[j]]['diurnal']

                if working_co2level not in biorepSCCs.keys():
                    biorepSCCs[working_co2level]={}
                if working_epoch not in biorepSCCs[working_co2level].keys():
                    biorepSCCs[working_co2level][working_epoch]=[]

                if growthA == growthB and diurnalA == diurnalB:
                    biorepSCCs[working_co2level][working_epoch].append(rho)

    # mirroring the values of the matrix
    for i in range(z):
        for j in range(z):
            if i > j:
                M[i,j]=M[j,i]

    # building the heatmap of correlations
    matplotlib.pyplot.imshow(M,interpolation='none',cmap='viridis',vmin=0.53,vmax=1.)
    cbar=matplotlib.pyplot.colorbar()
    cbar.set_label('SCC',size=24)
    matplotlib.pyplot.grid(False)
    matplotlib.pyplot.xticks(range(len(conditionNames)),conditionNames,rotation=90)
    matplotlib.pyplot.yticks(range(len(conditionNames)),conditionNames)
    matplotlib.pyplot.tight_layout()
    figureFileName='correlations.FPKM.{}.{}.png'.format(working_co2level,working_epoch)
    matplotlib.pyplot.savefig(figureFileName)
    matplotlib.pyplot.clf()
    
    return biorepSCCs

def expressionReader():

    '''
    this function reads the matrix of expression in FPKM into the format of a dictionary
    '''

    expression={}
    with open(expressionFile,'r') as f:
        header=f.readline()
        prelabels=header.split('\t')[1:]
        labels=[element.split('_')[0] for element in prelabels]
        next(f)
        for line in f:
            vector=line.split('\t')

            geneName=vector[0]
            expression[geneName]={}
            
            preValues=vector[1:]
            values=[float(element) for element in preValues]
            for i in range(len(values)):
                expression[geneName][labels[i]]=values[i]

    return expression

def metadataReader():

    '''
    this function returns a dictionary with the metadata of the expression data
    '''
    
    metaData={}

    with open(metaDataFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            if vector[4] != '':

                sampleID=vector[6]
                if vector[0] != '':
                    epoch=int(vector[0])
                growth=vector[1]
                diurnal=vector[2]
                co2=int(vector[3].replace(',',''))
                replicate=vector[4]
                preCollapse=bool(int(vector[5]))

                metaData[sampleID]={}
                metaData[sampleID]['growth']=growth
                metaData[sampleID]['epoch']=epoch
                metaData[sampleID]['diurnal']=diurnal
                metaData[sampleID]['co2']=co2
                metaData[sampleID]['replicate']=replicate
                metaData[sampleID]['pre-collapse']=preCollapse

    return metaData

def overalCorrelationPlotter(biorepSCCs):

    '''
    this function plots the overall correlation between biological replicates
    '''

    blur=0.05
    
    labelsDone=[]
    
    for level in biorepSCCs.keys():
        
        if level == 300:
            theColor='blue'
            theLabel='LC'
            sep=-0.1
        else:
            theColor='red'
            theLabel='HC'
            sep=0
            
        for epoch in biorepSCCs[level].keys():
            
            x=epoch
            y=[element for element in biorepSCCs[level][epoch]]

            # plotting values
            for i in range(len(y)):
                py=y[i]
                r=blur/2
                #r=blur*random.random()
                px=epoch+sep+r
                matplotlib.pyplot.plot(px,py,'o',color=theColor,alpha=0.33,mew=0.,ms=8)

            # plotting mean
            meanValue=numpy.mean(y)
            if theLabel not in labelsDone:
                matplotlib.pyplot.plot([x+sep,x+sep+blur],[meanValue,meanValue],'-',lw=4,color=theColor,label=theLabel)
                labelsDone.append(theLabel)
            else:
                matplotlib.pyplot.plot([x+sep,x+sep+blur],[meanValue,meanValue],'-',lw=4,color=theColor)

    # closing figure
    matplotlib.pyplot.legend(loc='lower left',fontsize=20)
    matplotlib.pyplot.ylim([0.53,1])
    matplotlib.pyplot.ylabel('Spearman rank CC',size=24)
    matplotlib.pyplot.xlabel('epoch',size=24)
    matplotlib.pyplot.tight_layout()
    figureFileName='correlationTrend.png'
    matplotlib.pyplot.savefig(figureFileName)
    matplotlib.pyplot.clf()

    return None

# 0. defining user variables
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
tau=5

# 1. reading files

# 1.1. reading metadata
print('reading metadata...')
metaData=metadataReader()

# 1.2. reading expression
print('reading expression...')
expression=expressionReader()
sortedGenes=list(expression.keys())
sortedGenes.sort()

# 2. computing correlation between biological replicates
print('computing matrix of correlations and building graphs...')

biorepSCCs={} # biorepSCCs[300,1000][0,1,2]=[]

biorepSCCs=correlationQuantifier(300,0,biorepSCCs)
biorepSCCs=correlationQuantifier(300,1,biorepSCCs)
biorepSCCs=correlationQuantifier(1000,0,biorepSCCs)
biorepSCCs=correlationQuantifier(1000,1,biorepSCCs)
biorepSCCs=correlationQuantifier(1000,2,biorepSCCs)

# 3. plotting overall correlations per co2 level and epoch
print('building figure of correlation trends...')
overalCorrelationPlotter(biorepSCCs)

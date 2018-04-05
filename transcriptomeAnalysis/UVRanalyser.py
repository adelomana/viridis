###
### This script compares the UVR-response genes in HC vs LC.
### It compares if there is consistent expression difference in epoch 0 vs 1 with respect to HC and LC on UVR-response genes.
###

import sys,numpy,os,pickle
import multiprocessing,multiprocessing.pool
import matplotlib,matplotlib.pyplot
import library

matplotlib.rcParams.update({'font.size':18,'font.family':'Arial','xtick.labelsize':14,'ytick.labelsize':14})
matplotlib.rcParams['pdf.fonttype']=42

def analysis(uvrGenes,verbose,plotting):

    '''
    This function computes the difference between Stage 2 and Stage 1 with respect to HC/LC fold change for UVR-response genes and others.
    '''

    # f.1. define variables to return
    effectsUVR=[]; effectsNOUVR=[]; effectsUVRdict={}

    # f.2. iterate over genes
    for geneName in sortedGeneNames:

        # f.2.1. find trajectory for LC
        co2level=300
        orderedExpressionLC=[]

        for epoch in epochs:
            for growth in growths:
                for diurnal in diurnals:

                    localExpression=[]
                    for sampleID in metadata.keys():
                        if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:
                            value=expression[geneName][sampleID]
                            localExpression.append(value)
                    orderedExpressionLC.append(localExpression)

        # f.2.2. find trajectory for HC
        co2level=1000
        orderedExpressionHC=[]

        for epoch in epochs:
            for growth in growths:
                for diurnal in diurnals:

                    localExpression=[]
                    for sampleID in metadata.keys():
                        if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:
                            value=expression[geneName][sampleID]
                            localExpression.append(value)
                    orderedExpressionHC.append(localExpression)

        # f.2.3. convert into array and compute means
        LC=numpy.array(orderedExpressionLC)+1
        HC=numpy.array(orderedExpressionHC)+1

        mLC=numpy.mean(LC,axis=1)
        mHC=numpy.mean(HC,axis=1)

        # f.2.4. discard any gene whose max expression does not reach 10 FPKM
        top=max([max(mLC),max(mHC)])
        if top < 10+1:
            pass
            #print('\t skipping transcript {} for low expression {:.4f}...'.format(geneName,top))
        else:

            # compute Stage 2 delta HC/LC
            deltaStage2=numpy.log2(mHC[4:]/mLC[4:])

            # compute Stage 1 delta HC/LC
            deltaStage1=numpy.log2(mHC[:4]/mLC[:4])

            # compute UVR effect
            effect=(sum(deltaStage2)/4)-(sum(deltaStage1)/4)

            # f.2.5. sort effects based on being UVR response genes or not
            if geneName in uvrGenes:
                effectsUVR.append(effect)
                effectsUVRdict[geneName]=effect
                # plotting if required
                if plotting == True:
                    figureFileName='figures/uvr.trajectory.{}.pdf'.format(geneName)
                    xTickLabels=['exp.AM','exp.PM','sta.AM','sta.PM','exp.AM','exp.PM','sta.AM','sta.PM']
                    x=[i for i in range(len(xTickLabels))]
                    matplotlib.pyplot.plot(x[:4],mLC[:4],':',color='blue',lw=2,label='LC')
                    matplotlib.pyplot.plot(x[:4],mHC[:4],':',color='red',lw=2,label='HC')
                    matplotlib.pyplot.plot(x[4:],mLC[4:],':',color='blue',lw=2)
                    matplotlib.pyplot.plot(x[4:],mHC[4:],':',color='red',lw=2)
                    matplotlib.pyplot.errorbar(x,mLC,yerr=stdDevLC,fmt='o',color='blue')
                    matplotlib.pyplot.errorbar(x,mHC,yerr=stdDevHC,fmt='o',color='red')
                    matplotlib.pyplot.fill_between(x[:4],mHC[:4],mLC[:4],facecolor='green',alpha=0.2,edgecolor='None')
                    matplotlib.pyplot.fill_between(x[4:],mHC[4:],mLC[4:],facecolor='orange',alpha=0.2,edgecolor='None')
                    matplotlib.pyplot.xticks(x,xTickLabels)
                    matplotlib.pyplot.ylabel('Expression (FPKM)')
                    matplotlib.pyplot.legend(markerscale=1.5,framealpha=1,loc=2,fontsize=12)
                    matplotlib.pyplot.tight_layout()
                    matplotlib.pyplot.savefig(figureFileName)
                    matplotlib.pyplot.clf()
            else:
               effectsNOUVR.append(effect)

    # f.3. print top 5 UVR effects
    if verbose == True:
        for element in sorted(effectsUVRdict,key=effectsUVRdict.get,reverse=True)[:5]:
            print('\t{}\t{}'.format(element,effectsUVRdict[element]))
    
    return effectsUVR,effectsNOUVR,effectsUVRdict

def crossValidation():

    '''
    This function performs the crossvalidation of number of outliers.
    '''

    # f.1. check for right jar, if not found, run a crossvalidation
    rightJar='uvr.crossvalidation.{}.iterations.pckl'.format(crossValidationIterations)
    if os.path.exists(rightJar) == True:
        # f.1.1. open the jar
        f=open(rightJar,'rb')
        hits=pickle.load(f)
        f.close()
        print('\t a {} crossvalidation recovered.'.format(len(hits)))
    else:
        # f.1.2. run the cross-validation
        hydra=multiprocessing.pool.Pool(numberOfThreads)
        tasks=[i for i in range(crossValidationIterations)]
        hits=hydra.map(outlierFinder,tasks)

        # save a jar
        f=open(rightJar,'wb')
        pickle.dump(hits,f)
        f.close()
        print('\t a {} crossvalidation pickled.'.format(len(hits)))

    # f.2. plot a histogram
    H=numpy.array(hits)
    pvalue=len(H[H>3])/len(H)
    print('\t crossvalidation p-value: {}'.format(pvalue))
    
    x,y=histogrammer(hits)
    cdfx=[]; cdfy=[]
    for i in range(len(x)):
        if y[i] != 0:
            cdfx.append(x[i])
            cdfy.append(y[i])
    print(numpy.cumsum(cdfy))
    matplotlib.pyplot.plot(cdfx,numpy.cumsum(cdfy),'o',mew=0,color='black')
    matplotlib.pyplot.axvline(x=3,color='red',ls=':',lw=2)
    
    matplotlib.pyplot.xlabel('Outlier rank')
    matplotlib.pyplot.ylabel('Cumulative probability')
    
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.savefig('uvr.crossvalidation.pdf')
    matplotlib.pyplot.clf()          
    return

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

def outlierFinder(tempo):

    '''
    This function is called by crossValidation and calls analysis to detect the number of outliers.
    '''

    # f.1. define a set of random elements from all genes the same size of UVR set
    randomSet=numpy.random.choice(sortedGeneNames,size=len(uvrGenes))

    # f.2. define how many of them pass a threshold
    outlierRank=0
    tempo1,tempo2,effectsRandomSetdict=analysis(randomSet,verbose=False)
    for geneName in effectsRandomSetdict:
        effect=effectsRandomSetdict[geneName]
        if effect > topThreshold:
            outlierRank=outlierRank+1

    return outlierRank

def plotter():

    '''
    This function plots two histograms, one for each set of gene values.
    '''

    fig, (axTop, axBottom) = matplotlib.pyplot.subplots(2,1, gridspec_kw = {'height_ratios':[9, 1]})
    #fig.subplots_adjust(hspace=0.03)
    
    # f.1. plot distribution of all genes
    x,y=histogrammer(effectsNOUVR)
    
    left=[[],[]]; center=[[],[]]; right=[[],[]]
    centers=[]
    halfDelta=(x[1]-x[0])/2

    for i in range(len(x)):
        acc=numpy.sum(y[:i+1])
        if acc < 0.025:
            left[0].append(x[i]); left[1].append(y[i])
        elif acc > 0.975:
            right[0].append(x[i]); right[1].append(y[i])
        else:
            center[0].append(x[i]); center[1].append(y[i])
            centers.append(x[i])

    top=max(centers)+halfDelta
    bottom=min(centers)-halfDelta

    print('\t thresholds: {} {}'.format(bottom,top))

    axTop.plot([top,top],[0,0.08],ls=':',color='red',lw=1)
    axTop.plot([bottom,bottom],[0,0.08],ls=':',color='red',lw=1)

    axTop.plot(left[0],left[1],'-',color='red',lw=1)
    axTop.fill_between(left[0],left[1],numpy.zeros(len(left[1])),facecolor='red',alpha=0.2,edgecolor='None')
    
    axTop.plot(center[0],center[1],'-',color='black',lw=1)
    axTop.fill_between(center[0],center[1],numpy.zeros(len(center[1])),facecolor='black',alpha=0.2,edgecolor='None')
    
    axTop.plot(right[0],right[1],'-',color='red',lw=1)
    axTop.fill_between(right[0],right[1],numpy.zeros(len(right[1])),facecolor='red',alpha=0.2,edgecolor='None')
    
    axTop.set_xticks([])
    axTop.set_xlim([-3.5,3.5])
    axTop.set_ylabel('Probability density function')

    # f.2. plotting UVR genes
    x,y=histogrammer(effectsUVR)
    for element in effectsUVR:
        if element > top:
            axBottom.plot([element],[-0.01],'o',color='red',alpha=0.5,mew=0)
        else:
            axBottom.plot([element],[-0.01],'o',color='black',alpha=0.5,mew=0)
            
    matplotlib.pyplot.xlabel('log$_2$ FC(HC/LC)$_{\mathrm{Stage \, 2}}$ - log$_2$ FC(HC/LC)$_{\mathrm{Stage \, 1}}$')
    axBottom.set_xticks([-3,-2,-1,0,1,2,3])
    axBottom.set_yticks([])
    axBottom.set_xlim([-3.5,3.5])
    
    # f.3. closing the figure
    matplotlib.pyplot.tight_layout(pad=0.4,w_pad=0.5,h_pad=0)
    matplotlib.pyplot.savefig('uvr.figure.pdf')
    matplotlib.pyplot.clf()

    return top,bottom

def uvrRespondingGenesReader():

    '''
    This function retrieves the names of UVR-responding genes.
    '''

    names=[]
    with open(uvrRespondingGenesFile,'r') as f:
        next(f)
        for line in f:
            vector=line.split(',')
            name='gene:Thaps{}'.format(vector[1])
            names.append(name)

    names.sort()

    return names

###
### MAIN
###

# 0. user defined variables
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
uvrRespondingGenesFile='/Volumes/omics4tb/alomana/projects/dtp/data/functionalAnnotation/UVR_Responsive_Genes.csv'

epochs=[0,1]
growths=['exp','sta']
diurnals=['AM','PM']

numberOfThreads=4
crossValidationIterations=10000

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
print('\t found {} UVR response genes...'.format(len(uvrGenes)))

# 2. perform analysis
print('performing analysis...')
effectsUVR,effectsNOUVR,effectsUVRdict=analysis(uvrGenes,verbose=True,plotting=False)

# 3. plot one histogram for each set of genes, one for UVR-response genes, another for the rest
print('plotting...')
topThreshold,bottomThreshold=plotter()

# 4. perform random sets of 101 genes and define the number of
print('running cross-validation...')
crossValidation()

# 5. final message
print('... done.')

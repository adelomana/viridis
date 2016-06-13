### this script locates all epochs samples into a diurnal vs. growth space characterized in epoch=0

import sys,numpy,matplotlib,random,pickle,math
import scipy.stats,scipy.spatial
import matplotlib.pyplot
from matplotlib.patches import Ellipse
matplotlib.rcParams['pdf.fonttype']=42 # this cryptical line is necessary for Illustrator compatibility of text saved as pdf

def addArrows():


    begin=(-0.4,1.5)
    end=(0.4,1.5)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    begin=(-0.4,-1.5)
    end=(0.4,-1.5)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    begin=(0.4,0.4)
    end=(-0.4,-0.4)
    matplotlib.pyplot.annotate("", xy=end, xytext=begin,arrowprops=dict(arrowstyle="->",alpha=0.2))

    return None

def bordersCalculator(descriptors,borders,flag):

    '''
    this function computes the borders of the descriptors in log10 scale
    '''

    # f.0. defining a variable for recapitulating the data borders necessary for mapping samples into new space
    borders[flag]={}
    
    # f.1. defining the expression values
    if flag == 'diurnal':
        expressionA=expressionRetriever(descriptors,flag,'AM')
        expressionB=expressionRetriever(descriptors,flag,'PM')        
    elif flag == 'growth':
        expressionA=expressionRetriever(descriptors,flag,'exp')
        expressionB=expressionRetriever(descriptors,flag,'sta')
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # f.2. recovering the info for the mapping samples
    listOfDescriptors=sorted(descriptors,key=descriptors.__getitem__,reverse=True) # ranking the descriptors
    for geneID in listOfDescriptors:
        borders[flag][geneID]={}
        
        # f.2.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        # f.2.2. saving the info for the mapping samples
        xa=numpy.min(logx); xb=numpy.median(logx); xc=numpy.max(logx); sdx=numpy.std(logx)
        ya=numpy.min(logy); yb=numpy.median(logy); yc=numpy.max(logy); sdy=numpy.std(logy)
        if xb > yb:
            center=((xa-yc)/2.)+yc
        else:
            center=((ya-xc)/2.)+xc
        # f.2.3. incorporating the data
        #! w=weightRankCalculator(geneID,listOfDescriptors)
        w=weightProportionalCalculator(geneID,flag)
        borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        
    return borders

def boxPlotGrapher(descriptors,borders,flag):

    '''
    this function plots the expression of the descriptors in log10 scale using boxplots
    '''

    # f.0. retrieving the expression data
    if flag == 'diurnal':
        expressionA=expressionRetriever(descriptors,flag,'AM')
        expressionB=expressionRetriever(descriptors,flag,'PM')        
    elif flag == 'growth':
        expressionA=expressionRetriever(descriptors,flag,'exp')
        expressionB=expressionRetriever(descriptors,flag,'sta')
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    listOfDescriptors=sorted(descriptors,key=descriptors.__getitem__,reverse=True) # ranking the descriptors
    
    # f.1. computing the figures for specific AM/exp or PM/sta descriptors
    # AM/exp
    trend='AM-exp'
    boxPlotMaker(listOfDescriptors,expressionA,expressionB,flag,trend)

    # PM/sta
    trend='PM-sta'
    boxPlotMaker(listOfDescriptors,expressionA,expressionB,flag,trend)

    return None

def boxPlotMaker(listOfDescriptors,expressionA,expressionB,flag,trend):

    '''
    this function builds the figure for the distribution of the descriptors
    '''

    boxPlotPosition=0
    names=[]
    figureCount=0
    borA=[]
    borB=[]

    if flag == 'diurnal':
        optimalRange=36
        optimalFontSize=15
    else:
        optimalRange=38
        optimalFontSize=13
    
    for geneID in listOfDescriptors:
        
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        xa=numpy.min(logx); xb=numpy.median(logx); xc=numpy.max(logx); sdx=numpy.std(logx)
        ya=numpy.min(logy); yb=numpy.median(logy); yc=numpy.max(logy); sdy=numpy.std(logy)

        if trend == 'AM-exp':
             if xb > yb:
                boxPlotPosition=boxPlotPosition+1
                borA.append(xa); borB.append(yc)

                bp=matplotlib.pyplot.boxplot([logx],positions=[boxPlotPosition],patch_artist=True)
                if flag == 'diurnal':
                    setBoxColors(bp,'orange')
                else:
                    setBoxColors(bp,'#0571b0')
                bp=matplotlib.pyplot.boxplot([logy],positions=[boxPlotPosition],patch_artist=True)
                if flag == 'diurnal':
                    setBoxColors(bp,'darkgreen')
                else:
                    setBoxColors(bp,'#ca0020')

                name=geneID.split('Thaps')[1]
                names.append(name)
                
        elif trend == 'PM-sta':
            if xb < yb:
                boxPlotPosition=boxPlotPosition+1
                borA.append(xc); borB.append(ya)

                bp=matplotlib.pyplot.boxplot([logx],positions=[boxPlotPosition],patch_artist=True)
                if flag == 'diurnal':
                    setBoxColors(bp,'orange')
                else:
                    setBoxColors(bp,'#0571b0')
                bp=matplotlib.pyplot.boxplot([logy],positions=[boxPlotPosition],patch_artist=True)
                if flag == 'diurnal':
                    setBoxColors(bp,'darkgreen')
                else:
                    setBoxColors(bp,'#ca0020')

                name=geneID.split('Thaps')[1]
                names.append(name)
        else:
            print 'error from boxPlotMaker about trend selection. Exiting...'
            sys.exit()

        # closing figure because of number of descriptors
        if boxPlotPosition >= optimalRange:
            matplotlib.pyplot.fill_between(range(1,len(names)+1),borA,borB,facecolor='black',alpha=0.2,edgecolor='None')
            
            matplotlib.pyplot.xlim([0,boxPlotPosition+1])
            matplotlib.pyplot.ylim([-0.2,5.])
            theXticks=range(boxPlotPosition)
            theXticksPosition=[element+1 for element in theXticks]
            matplotlib.pyplot.xticks(theXticksPosition,names,rotation=90,fontsize=optimalFontSize)
            matplotlib.pyplot.yticks(fontsize=optimalFontSize)
            matplotlib.pyplot.ylabel('log$_{10}$ FPKM',fontsize=optimalFontSize)
            matplotlib.pyplot.tight_layout(pad=2.5)
            matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
            matplotlib.pyplot.tick_params(axis='y',which='both',right='off')
            
            if flag == 'diurnal':
                matplotlib.pyplot.plot([-1],[-1],color='orange',lw=2,label='light')
                matplotlib.pyplot.plot([-1],[-1],color='darkgreen',lw=2,label='dark')
            else:
                matplotlib.pyplot.plot([-1],[-1],color='#0571b0',lw=2,label='early')
                matplotlib.pyplot.plot([-1],[-1],color='#ca0020',lw=2,label='late')
            matplotlib.pyplot.legend(fontsize=optimalFontSize)

            figureName='figures/boxplots.%s.%s.%s.pdf'%(flag,trend,int(figureCount))
            matplotlib.pyplot.savefig(figureName)
            
            matplotlib.pyplot.clf()

            boxPlotPosition=0
            names=[]
            figureCount=figureCount+1
            borA=[]
            borB=[]

    # final closing of figure
    matplotlib.pyplot.fill_between(range(1,len(names)+1),borA,borB,facecolor='black',alpha=0.2,edgecolor='None')
    
    matplotlib.pyplot.xlim([0,boxPlotPosition+1])
    matplotlib.pyplot.ylim([-0.2,5.])
    theXticks=range(boxPlotPosition)
    theXticksPosition=[element+1 for element in theXticks]
    matplotlib.pyplot.xticks(theXticksPosition,names,rotation=90,fontsize=optimalFontSize)
    matplotlib.pyplot.yticks(fontsize=optimalFontSize)
    matplotlib.pyplot.ylabel('log$_{10}$ FPKM',fontsize=optimalFontSize)
    matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off')

    if flag == 'diurnal':
        matplotlib.pyplot.plot([-1],[-1],color='orange',lw=2,label='light')
        matplotlib.pyplot.plot([-1],[-1],color='darkgreen',lw=2,label='dark')
    else:
        matplotlib.pyplot.plot([-1],[-1],color='#0571b0',lw=2,label='early')
        matplotlib.pyplot.plot([-1],[-1],color='#ca0020',lw=2,label='late')
    matplotlib.pyplot.legend(fontsize=optimalFontSize)
            
    figureName='figures/boxplots.%s.%s.%s.pdf'%(flag,trend,int(figureCount))
    matplotlib.pyplot.savefig(figureName)
    matplotlib.pyplot.clf()

    return None

def coefficientOfVariationDistributionCalculator(workingSamples):

    '''
    this function returns the distribution of coefficient of variations for each set of descriptors, normalized by their weight
    '''

    distribution=[]
    for dimension in borders.keys():
        for descriptor in borders[dimension].keys():
            w=borders[dimension][descriptor][7] # weight
            v=[]
            for sampleID in workingSamples:
                s=numpy.log10(expression[descriptor][sampleID]+1.)
                v.append(s)
            if max(v) == 0.:
                cv=0.
            else:
                cv=numpy.std(v)/numpy.mean(v)
            distribution.append(cv)

    # dealing with zero values and log transforming
    lowerValue=min([element for element in distribution if element != 0.])
    formattedDistribution=[lowerValue if (element == 0.) else element for element in distribution]
    logDistribution=numpy.log10(formattedDistribution)

    return logDistribution

def descriptorsFilter(descriptors,flag):

    '''
    this function selects the descriptors that have at least separation between max and min values into the two cases in log space
    '''
    
    # 1. defining the expression values
    if flag == 'diurnal':
        selectedExpressionA=expressionRetriever(descriptors,flag,'AM')
        selectedExpressionB=expressionRetriever(descriptors,flag,'PM')
    elif flag == 'growth':
        selectedExpressionA=expressionRetriever(descriptors,flag,'exp')
        selectedExpressionB=expressionRetriever(descriptors,flag,'sta')
    else:
        print 'error handling flags from boxPlotGrapher'
        sys.exit()

    # 2. computing the distances
    separations={}
    for geneID in descriptors: 
        
        # 2.1 transforming into log10 scale
        x=numpy.array(selectedExpressionA[geneID])
        y=numpy.array(selectedExpressionB[geneID])

        x=x+1.
        y=y+1.

        logx=numpy.log10(x)
        logy=numpy.log10(y)

        # computing a separation threshold: it should be larger than the average of the standard deviations of the two distributions
        averageSD=0.5*numpy.std(logx) + 0.5*numpy.std(logy)
        
        if numpy.mean(logx) > numpy.mean(logy):
            separation=numpy.min(logx)-numpy.max(logy)
        elif numpy.mean(logy) > numpy.mean(logx):
            separation=numpy.min(logy)-numpy.max(logx)
        else:
            print 'error computing the separation between the two distributions from descriptorsFilter'
            sys.exit()
        if separation > averageSD:
            separations[geneID]=separation

    return separations

def descriptorsRetriever(flag):

    '''
    this function reads the output of cuffdiff and select the genes that pass the following rule: 1. log2 fold > 1., and 2. statistical significance
    '''

    descriptorsFolds={}
    descriptorsSignificances={}
    
    inputFile=cuffdiffDir+'%s/gene_exp.diff'%flag
    with open(inputFile,'r') as f:
        f.next()
        for line in f:
            vector=line.split('\t')
            geneName=vector[0]
            significance=vector[-1].replace('\n','')
            foldChange=abs(float(vector[-5]))
            qValue=float(vector[-2])
            if significance == 'yes' and foldChange > 1.: # the selection rule is log2 fold > 1.0 and statistically different
                descriptorsFolds[geneName]=foldChange
                descriptorsSignificances[geneName]=qValue
    # sorting descriptors
    listOfDescriptors=sorted(descriptorsFolds,key=descriptorsFolds.__getitem__,reverse=True)

    return listOfDescriptors

def descriptorsWriter(selected,flag):

    '''
    this function writes the gene descriptors into a file
    '''

    # f.1. differentiating phases of descriptors
    if flag == 'diurnal':
        f=open('descriptors/epoch0_diurnalDescriptors_AM_ranked.txt','w')
        g=open('descriptors/epoch0_diurnalDescriptors_PM_ranked.txt','w')
    elif flag == 'growth':
        f=open('descriptors/epoch0_growthDescriptors_exp_ranked.txt','w')
        g=open('descriptors/epoch0_growthDescriptors_sta_ranked.txt','w')

    sortedList=sorted(selected,key=selected.__getitem__,reverse=True)
    for element in sortedList:
        m=borders[flag][element][1]
        n=borders[flag][element][4]
        if m > n:
            f.write('%s\n'%element)
        else:
            g.write('%s\n'%element)

    f.close()
    g.close()

    # f.2. adding the relative importance of each descriptor
    if flag == 'diurnal':
        f=open('descriptors/epoch0_diurnalDescriptors_ranked.txt','w')
        header1='median expression light (FPKM)'
        header2='median expression dark (FPKM)'
    elif flag == 'growth':
        f=open('descriptors/epoch0_growthDescriptors_ranked.txt','w')
        header1='median expression early (FPKM)'
        header2='median expression late (FPKM)'

    accW=0.
    rank=0
    allW=[]
    
    f.write('#rank\tgeneID\t%s\t%s\tweight\taccumulated weight\n'%(header1,header2))
    for element in sortedList:
        rank=rank+1
        localW=borders[flag][element][-3]
        accW=accW+localW
        allW.append(localW)
        expressionA=(10**(borders[flag][element][1]))-1.
        expressionB=(10**(borders[flag][element][4]))-1.
        
        f.write('%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\n'%(rank,element,expressionA,expressionB,localW,accW))
    f.close()

    # f.3. making a plot of the distribution of weights
    figureFile='figures/weightDistribution.%s.pdf'%flag
    x=range(1,len(allW)+1)
    matplotlib.pyplot.plot(x,allW,'o-',color='black')
    matplotlib.pyplot.xlabel('rank')
    matplotlib.pyplot.ylabel('relative weight')
    matplotlib.pyplot.savefig(figureFile)
    matplotlib.pyplot.clf()

    # f.4 pickle the descriptors for other tools, like GSE_Mapper.py
    jar=flag+'.pckl'
    f=open(jar,'w')
    pickle.dump(sortedList,f)
    f.close()

    return None

def distanceQuantifier(flag):

    '''
    this function computes the distance between replicates per epoch at each sample point
    '''
    
    co2level=int(flag)
    orderedSamples=sampleOrderer(co2level)

    distances={}
    epochLabel=0
    for timepoint in orderedSamples.keys():
        workingSamples=orderedSamples[timepoint]

        if (timepoint-1)%4 == 0:
            epochLabel=epochLabel+1
            distances[epochLabel]=[]

        # computing positions
        positions=[]
        for sampleID in workingSamples:
            x,y=newCoordinateCalculator(sampleID)
            positions.append([x,y])

        # computing distances
        for i in range(len(positions)):
            for j in range(len(positions)):
                if i < j:
                    d=scipy.spatial.distance.euclidean(positions[i],positions[j])
                    distances[epochLabel].append(d)

    # plotting distances
    x = distances.keys()
    y = distances.values()

    violinParts=matplotlib.pyplot.violinplot(y,x,showmeans=True,showextrema=False)

    timeLabels=['epoch.'+str(element) for element in x]
    colors=[]
    for timeLabel in timeLabels:
        if '.1' in timeLabel:
            colors.append('blue')
        elif '.2' in timeLabel:
            colors.append('green')
        elif '.3' in timeLabel:
            colors.append('red')
        else:
            print 'error assigning violin colors. exiting...'
            sys.exit()

    for i in range(len(timeLabels)):
        violinParts['bodies'][i].set_linewidths(None)
        violinParts['bodies'][i].set_facecolor(colors[i])
    violinParts['cmeans'].set_edgecolor('black')

    matplotlib.pyplot.ylabel('replicate distance')
    matplotlib.pyplot.xticks(range(1,len(timeLabels)+1),timeLabels,rotation=-45)
    matplotlib.pyplot.ylim([0,1.2])
    matplotlib.pyplot.tight_layout()

    fileName='figures/figure.distances.%s.pdf'%(flag)
    matplotlib.pyplot.savefig(fileName)
    matplotlib.pyplot.clf()

    # checking for significance
    for i in range(len(x)):
        for j in range(len(x)):
            if i < j:
                a=distances[i+1]
                b=distances[j+1]
                statistic,pvalue=scipy.stats.mannwhitneyu(a,b)
                print flag,i+1,j+1,statistic,pvalue


    return None

def ellipseSizeCalculator(flag1,flag2):

    '''
    this function calculates size of ellipsoids
    '''

    if flag1 == 'diurnal':
        rankedDimensions=sorted(diurnalFilteredDescriptors,key=diurnalFilteredDescriptors.__getitem__,reverse=True)
    elif flag1 == 'growth':
        rankedDimensions=sorted(growthFilteredDescriptors,key=growthFilteredDescriptors.__getitem__,reverse=True)

    averageDispersion=0.
    for transcript in rankedDimensions:

        if flag2 == 'AM' or flag2 == 'exp':
        
            mean=borders[flag1][transcript][1]
            deviation=borders[flag1][transcript][-2]
            top=borders[flag1][transcript][2]
            weight=borders[flag1][transcript][-3]
            
        elif flag2 == 'PM' or flag2 == 'sta':
            
            mean=borders[flag1][transcript][4]
            deviation=borders[flag1][transcript][-1]
            top=borders[flag1][transcript][5]
            weight=borders[flag1][transcript][-3]
     
        
        stretch=top-mean
        if stretch != 0.:
            position=0.5*(deviation*1.96)/stretch
        else:
            position=0.
        value=position*weight
        averageDispersion=averageDispersion+value

    print 'dispersion found for',flag1,flag2,'conditions: ',averageDispersion

    return averageDispersion

def expressionReader():

    '''
    this function reads the matrix of expression in FPKM into the format of a dictionary
    '''

    expression={}
    with open(expressionFile,'r') as f:
        header=f.readline()
        prelabels=header.split('\t')[1:]
        labels=[element.split('_')[0] for element in prelabels]
        f.next()
        for line in f:
            vector=line.split('\t')

            geneName=vector[0]
            expression[geneName]={}
            
            preValues=vector[1:]
            values=[float(element) for element in preValues]
            for i in range(len(values)):
                expression[geneName][labels[i]]=values[i]

    return expression

def expressionRetriever(descriptors,flag,condition):

    '''
    this function returns the expression values for a set of genes under a specific condition
    '''

    # 1. selecting the acceptable samples
    selectedSamples=[]
    for sampleID in metaData.keys():
        if metaData[sampleID]['epoch'] == 0:
            if metaData[sampleID][flag] == condition:
                selectedSamples.append(sampleID)
                
    # 2. defining the expression values
    selectedExpression={}
    for geneName in descriptors:
        selectedExpression[geneName]=[]
        for sampleID in selectedSamples:
            value=expression[geneName][sampleID]
            selectedExpression[geneName].append(value)

    return selectedExpression

def loadCalculator(sampleID,flag):

    '''
    this function computes the value of the sample in the new axis 
    '''

    averageLoad=0.
    A=[]
    W=[]
    
    for descriptor in borders[flag].keys():
        s=numpy.log10(expression[descriptor][sampleID]+1.) # sample expression level

        a=borders[flag][descriptor][0] # min of positive
        b=borders[flag][descriptor][1] # median of positive
        c=borders[flag][descriptor][2] # max of positive

        d=borders[flag][descriptor][3] # min of negative
        e=borders[flag][descriptor][4] # median of negative
        f=borders[flag][descriptor][5] # max of negative

        NML=borders[flag][descriptor][6] # No Man's Land
        w=borders[flag][descriptor][7] # weight
        
        # differentiate if it's positive (AM,exp) or negative (PM,sta) or miss-regulated
        positive=None
        missR=None
        if b > e:
            if s > NML:
                positive=True
            else:
                positive=False
            if s < a and s > f:
                missR=True
            else:
                missR=False
        if e > b:
            if s > NML:
                positive=False
            else:
                positive=True
            if s < d and s > c:
                missR=True
            else:
                missR=False
        
        # working with the misregulated samples
        if missR == True:
            if b > e:
                if s > NML:
                    stretch=abs(a-NML)
                    value=(s-NML)/stretch
                else:
                    stretch=abs(f-NML)
                    value=-abs(s-NML)/stretch
            else:
                if s > NML:
                    stretch=abs(d-NML)
                    value=-abs(s-NML)/stretch
                else:
                    stretch=abs(c-NML)
                    value=abs(s-NML)/stretch
                        
        # dealing with values within previously observed
        else:
            if positive == True: # it is a light sample
                # assuming light boxplot is above
                if b > e:
                    if s > b:
                        stretch=abs(c-b)
                        deviation=abs(s-b)
                        value=1.5+0.5*(deviation/stretch)
                    else:
                        stretch=abs(b-a)
                        deviation=abs(s-b)
                        value=1.5-0.5*(deviation/stretch)
                # assuming light boxplot is below
                else:
                    if s > b:
                        stretch=abs(c-b)
                        deviation=abs(s-b)
                        value=1.5-0.5*(deviation/stretch)
                    else:
                        stretch=abs(b-a)
                        deviation=abs(s-b)
                        if stretch != 0.:
                            value=1.5+0.5*(deviation/stretch)
                        else:
                            value=2.
            else: # it is a dark sample
                # assuming light boxtplot is above
                if b > e:
                    if s > e:
                        stretch=abs(f-e)
                        deviation=abs(s-e)
                        value=-1.5+0.5*(deviation/stretch)
                    else:
                        stretch=abs(e-d)
                        deviation=abs(s-e)
                        value=-1.5-0.5*(deviation/stretch)
                # assuming light boxplot is below
                else:
                    if s > e:
                        stretch=abs(f-e)
                        deviation=abs(s-e)
                        value=-1.5-0.5*(deviation/stretch)
                    else:
                        stretch=abs(e-d)
                        deviation=abs(s-e)
                        value=-1.5+0.5*(deviation/stretch)

        # weighting the value
        term=w*value
        averageLoad=averageLoad+term

        # returning the full list of descriptors
        A.append(value)
        W.append(w)

    return averageLoad,A,W

def metadataReader():

    '''
    this function returns a dictionary with the metadata of the expression data
    '''
    
    metaData={}

    with open(metaDataFile,'r') as f:
        f.next()
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

def newCoordinateCalculator(sampleID):

    '''
    this function computes the new coordinates of a sample given the borders calculated from the filtered descriptors
    '''

    x,xA,xW=loadCalculator(sampleID,'diurnal')
    y,yA,yW=loadCalculator(sampleID,'growth')

    # converting to inverse values for figure representation
    x=-x

    return x,y

def newSpaceMapper(flag):

    '''
    this function plots the samples into a new space
    '''

    preselectedSamples=[sampleID for sampleID in metaData.keys() if metaData[sampleID]['co2'] == int(flag)]
    print 'selected ', len(preselectedSamples), 'samples for plotting on ',flag, 'condition.'

    # starting the figure
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(111)
    
    for sampleID in preselectedSamples:
        
        x,y=newCoordinateCalculator(sampleID)
    
        theSize=10
        theAlpha=.85

        # defining the color depending on the light/dark
        if metaData[sampleID]['diurnal'] == 'AM':
            if metaData[sampleID]['replicate'] == 'A':
                theColor='chocolate'
            elif metaData[sampleID]['replicate'] == 'B':
                theColor='orange'
            elif metaData[sampleID]['replicate'] == 'C':
                theColor='orangered'
            else:
                print 'error defining AM replicate at newSpaceMapper. exiting...'
                sys.exit()
        elif metaData[sampleID]['diurnal'] == 'PM':
            if metaData[sampleID]['replicate'] == 'A':
                theColor='olive'
            elif metaData[sampleID]['replicate'] == 'B':
                theColor='green'
            elif metaData[sampleID]['replicate'] == 'C':
                theColor='lightseagreen'
            else:
                print 'error defining PM replicate at newSpaceMapper. exiting...'
                sys.exit()
        else:
            print 'error while defining the color from main'
            sys.exit()

        # defining the marker type depending on epoch
        if metaData[sampleID]['epoch'] == 0:
            theMarker='o'
        elif metaData[sampleID]['epoch'] == 1:
            theMarker='s'
        elif metaData[sampleID]['epoch'] == 2:
            theMarker='^'
        else:
            print 'error while defining the marker from main'
            sys.exit()

        # defining the facecolor and markeredgecolor depending on exp/sta
        if metaData[sampleID]['growth'] == 'exp':
            theMFC='None'; theMEC=theColor
        elif metaData[sampleID]['growth'] == 'sta':
            theMFC=theColor; theMEC='None'
        else:
            print 'error while defining the marker from main'
            sys.exit()
      
        ax.plot(x,y,marker=theMarker,mew=1,color=theColor,ms=theSize,alpha=theAlpha,mfc=theMFC,mec=theColor,zorder=10)

        # marking the pre-collapse samples
        if metaData[sampleID]['pre-collapse'] == True and metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['replicate'] == 'B':
            ax.plot(x,y,marker='*',color='black',ms=3,zorder=10)

    # plotting the ellipses
    a=ellipseSizeCalculator('diurnal','AM')
    b=ellipseSizeCalculator('diurnal','PM')
    c=ellipseSizeCalculator('growth','exp')
    d=ellipseSizeCalculator('growth','sta')

    e=matplotlib.patches.Ellipse(xy=(1.5,1.5),width=b,height=c,edgecolor='darkgreen',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    e=matplotlib.patches.Ellipse(xy=(1.5,-1.5),width=b,height=d,edgecolor='darkgreen',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    
    e=matplotlib.patches.Ellipse(xy=(-1.5,1.5),width=a,height=c,edgecolor='orange',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    e=matplotlib.patches.Ellipse(xy=(-1.5,-1.5),width=a,height=d,edgecolor='orange',fc='None',lw=1,alpha=1.,ls='--')
    ax.add_patch(e)

    addArrows()

    # setting ranges
    matplotlib.pyplot.xlim([-2.,2.])
    matplotlib.pyplot.ylim([-2.,2.])
    matplotlib.pyplot.xticks([-1.5,1.5],['light','dark'],fontsize=24)
    matplotlib.pyplot.yticks([-1.5,1.5],['late','early'],fontsize=24)
    matplotlib.pyplot.xlabel('diurnal cycle',fontsize=24)
    matplotlib.pyplot.ylabel('growth phase',fontsize=24)
    matplotlib.pyplot.title(flag+' ppm',fontsize=28)

    # defining health zones
    matplotlib.pyplot.plot([-2,-1],[1,1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([-1,-1],[1,2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([2,1],[1,1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([1,1],[1,2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([2,1],[-1,-1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([1,1],[-1,-2],color='black',alpha=0.2)

    matplotlib.pyplot.plot([-2,-1],[-1,-1],color='black',alpha=0.2)
    matplotlib.pyplot.plot([-1,-1],[-1,-2],color='black',alpha=0.2)

    # defining the grid lines
    matplotlib.pyplot.plot([-1.5,-1.5],[1,2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([-2,-1],[1.5,1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([-1.5,-1.5],[-1,-2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([-2,-1],[-1.5,-1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([1.5,1.5],[-1,-2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([2,1],[-1.5,-1.5],color='black',ls=':',alpha=0.4)

    matplotlib.pyplot.plot([1.5,1.5],[1,2],color='black',ls=':',alpha=0.4)
    matplotlib.pyplot.plot([2,1],[1.5,1.5],color='black',ls=':',alpha=0.4)

    # defining the misregulation zone
    matplotlib.pyplot.plot([-1,-1],[-1,1],color='black',alpha=0.7,lw=2.)
    matplotlib.pyplot.plot([-1,1],[1,1],color='black',alpha=0.7,lw=2.)
    matplotlib.pyplot.plot([1,1],[1,-1],color='black',alpha=0.7,lw=2.)
    matplotlib.pyplot.plot([-1,1],[-1,-1],color='black',alpha=0.7,lw=2.)

    # aspect
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    
    matplotlib.pyplot.savefig('figures/sampleLocation.%s.pdf'%(str(int(flag))))
    matplotlib.pyplot.clf()
    
    return None

def newSpaceProbabilisticMapper(flag):

    '''
    this function plots the samples into a new space considering in a probabilistic way each of the descriptors
    '''

    preselectedSamples=[sampleID for sampleID in metaData.keys() if metaData[sampleID]['co2'] == int(flag)]
    print 'selected ', len(preselectedSamples), 'samples for plotting on ',flag, 'condition.'

    # starting the figure
    epochT={}
    epochGP={}
    for sampleID in preselectedSamples:
        for i in range(3):
            
            epochT[i]=[]
            epochGP[i]=[]
                

        x,y=newCoordinateCalculator(sampleID)
        P,Q,grid=probabilisticCoordinateCalculator(sampleID)

        m=[]
        for i in range(len(grid)):
            v=[]
            for j in range(len(grid)):
                product=P[i]*Q[j]
                #if product == 0.:
                #    product=float('nan')

                v.append(product)
            m.append(v)
        m=numpy.array(m)
        t=numpy.transpose(m)

        print numpy.max(t)

        # plotting the figure
        gridPosition=[20+x*10,20+y*10]
        matplotlib.pyplot.imshow(t,interpolation='none',cmap='Blues',origin='lower',extent=(0,40,0,40)) # kaiser is a good interpolation
        matplotlib.pyplot.plot(gridPosition[0],gridPosition[1],'o',color='red',mew=0.)
        matplotlib.pyplot.hold(True)

    # setting ranges and labels
    matplotlib.pyplot.xlim([0,40.])
    matplotlib.pyplot.ylim([0,40.])
    matplotlib.pyplot.xlabel('diurnal cycle',fontsize=24)
    matplotlib.pyplot.ylabel('growth phase',fontsize=24)

    # aspect
    matplotlib.pyplot.tight_layout()
    matplotlib.pyplot.axes().set_aspect('equal')
    
    
    matplotlib.pyplot.savefig('figures/sampleLocation.prob.%s.png'%(str(int(flag))))
    matplotlib.pyplot.clf()

    sys.exit()

    return None

def oneDimensionTimeMapper():

    '''
    this function maps samples into 1D space and plots them as a time series
    '''

    epsilon=0.1
    theSize=15
    theAlpha=.6
    
    # f.1. building a structure holding samples as a time series
    co2levels=[300,1000]
    epochs=[0,1,2]
    growths=['exp','sta']
    diurnals=['AM','PM']

    for co2level in co2levels:
        orderedSamples={}
        time=0
        for epoch in epochs:
            for growth in growths:
                for diurnal in diurnals:
                    time=time+1
                    for sampleID in metaData.keys():
                        if metaData[sampleID]['co2'] == co2level and metaData[sampleID]['epoch'] == epoch and metaData[sampleID]['growth'] == growth and metaData[sampleID]['diurnal'] == diurnal:
                            if time in orderedSamples:
                                orderedSamples[time].append(sampleID)
                            else:
                                orderedSamples[time]=[sampleID]
            
        # f.2. computing distances and building trajectories for each co2 condition
        timePoints=orderedSamples.keys()
        timePoints.sort()
        for timePoint in timePoints:
            samples=orderedSamples[timePoint]

            yPoints=[]
            for sampleID in samples:
                
                yPos=distanceCalculator(sampleID)

                # selecting position and color
                if metaData[sampleID]['replicate'] == 'A':
                    xPos=timePoint-epsilon
                    if co2level == 300:
                        theColor='lightblue'
                    else:
                        theColor='salmon'
                elif metaData[sampleID]['replicate'] == 'B':
                    xPos=timePoint
                    if co2level == 300:
                        theColor='blue'
                    else:
                        theColor='red'
                elif metaData[sampleID]['replicate'] == 'C':
                    xPos=timePoint+epsilon
                    if co2level == 300:
                        theColor='darkblue'
                    else:
                        theColor='darkred'
                else:
                    print 'error choosing the replicates from oneDimensionTimeMapper. exiting...'
                    sys.exit()
            
                # actual plotting the dots
                matplotlib.pyplot.plot(xPos,yPos,marker='.',color=theColor,ms=theSize,alpha=theAlpha,mew=0.,zorder=1)
                yPoints.append(yPos)
            # plotting the bars
            average=numpy.mean(yPoints)
            if co2level == 300:
                theColor='blue'
            else:
                theColor='red'
            matplotlib.pyplot.plot([timePoint-epsilon,timePoint+epsilon],[average,average],'-',color=theColor,lw=2,zorder=2)
             
    # f.3. finishing up the figure
    matplotlib.pyplot.plot([1,12],[1.5,1.5],color='black',ls=':')
    matplotlib.pyplot.plot([1,12],[1.,1.],color='magenta',ls='--')
    matplotlib.pyplot.plot([1,12],[-1.,-1.],color='magenta',ls='--')
    matplotlib.pyplot.plot([4.5,4.5],[-1.5,2],color='black',alpha=0.2,ls='--')
    matplotlib.pyplot.plot([8.5,8.5],[-1.5,2],color='black',alpha=0.2,ls='--')
    
    matplotlib.pyplot.xlabel('epoch',fontsize=24)
    matplotlib.pyplot.ylabel('state',fontsize=24)
    matplotlib.pyplot.subplots_adjust(left=0.2, right=0.95, top=0.9, bottom=0.11)

    matplotlib.pyplot.xticks((2.5,6.5,10.5),('0','1','2'))
    matplotlib.pyplot.yticks((-1.5,0.,1.5),('inverse','misregulation','expected'))

    matplotlib.pyplot.xlim([0,13])
    matplotlib.pyplot.ylim([-1.5,2])
    
    matplotlib.pyplot.savefig('figures/distance.pdf')
    matplotlib.pyplot.clf()

    return None

def setBoxColors(bp,theColor):

    '''
    this function access the elements of a boxplot and colors them appropriately
    '''

    matplotlib.pyplot.setp(bp['boxes'],color=theColor)
    matplotlib.pyplot.setp(bp['caps'],color='None')
    matplotlib.pyplot.setp(bp['whiskers'],color=theColor,ls='-')
    matplotlib.pyplot.setp(bp['fliers'],markeredgecolor=theColor,marker='+')
    matplotlib.pyplot.setp(bp['medians'],color=theColor)    

    return None

def sampleOrderer(co2level):

    '''
    this function order samples timewise
    '''

    epochs=[0,1,2]
    growths=['exp','sta']
    diurnals=['AM','PM']

    orderedSamples={}
    time=0
    timeLabels=[]
    for epoch in epochs:
        for growth in growths:
            for diurnal in diurnals:
                time=time+1
                for sampleID in metaData.keys():
                    if metaData[sampleID]['co2'] == co2level and metaData[sampleID]['epoch'] == epoch and metaData[sampleID]['growth'] == growth and metaData[sampleID]['diurnal'] == diurnal:

                        timeLabel=diurnal+'.'+growth+'.'+str(epoch+1)
                        if timeLabel not in timeLabels:
                            timeLabels.append(timeLabel)
                            
                        if time in orderedSamples:
                            orderedSamples[time].append(sampleID)
                        else:
                            orderedSamples[time]=[sampleID]

    return orderedSamples

def weightedHistogrammer(scores,scoresWeights):

    '''
    this function computes a weighted histogram of a probability density distribution
    '''

    n,bins=numpy.histogram(scores,weights=scoresWeights,bins=20,range=(min(scores),max(scores)))

    z=[]
    halfBin=(bins[1]-bins[0])/2.
    for bin in bins:
        center=bin+halfBin
        z.append(center)
    z.pop()

    p=[]
    p=numpy.array(n)
    p=p/float(sum(p))

    return p,z

def weightProportionalCalculator(geneID,flag):

    '''
    this function computes the weight of the descriptor based on the empty space between the descriptors
    '''

    if flag == 'diurnal':
        sumSpaces=sum(diurnalFilteredDescriptors.values())
        value=diurnalFilteredDescriptors[geneID]
    elif flag == 'growth':
        sumSpaces=sum(growthFilteredDescriptors.values())
        value=growthFilteredDescriptors[geneID]

    weight=value/sumSpaces

    return weight

def weightRankCalculator(geneID,descriptors):

    '''
    this function computes the weight of the descriptor based on its rank
    '''

    inverseRank=len(descriptors)-descriptors.index(geneID)
    sumOfRanks=sum(numpy.arange(1.,len(descriptors)+1.))
    weight=float(inverseRank)/sumOfRanks

    return weight

### MAIN

# 0. preliminaries
print
print 'welcome to sampleMapper'
print
print 'initializing variables...'
# 0.1. user defined variables and paths
cuffdiffDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
expressionFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/allSamples/genes.fpkm_table.v2.txt'
metaDataFile='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'

boxplotPlotting=False
time300=numpy.array([1.375,1.625,3.375,3.708333333,5.291666667,5.708333333,7.333333333,7.75])
time1000=numpy.array([1.375,1.625,3.375,3.708333333,5.291666667,5.708333333,7.333333333,7.75,15.45833333,15.79166667,17.41666667,17.79166667])

# 0.2. reading metadata
metaData=metadataReader()

# 0.3. reading expression
expression=expressionReader()

# 1. selecting descriptors
print
print 'selecting state descriptors...'

# 1.1. recover the DET genes
print 'recovering DETs...'

diurnalDescriptors=descriptorsRetriever('light_epoch0')
print len(diurnalDescriptors),'diurnal descriptors detected.'

growthDescriptors=descriptorsRetriever('growth_epoch0')
print len(growthDescriptors),'growth descriptors detected.'

# 1.2. filtering the descriptors based on separation
print
print 'computing descriptors based on separation rules...'
diurnalFilteredDescriptors=descriptorsFilter(diurnalDescriptors,'diurnal')
growthFilteredDescriptors=descriptorsFilter(growthDescriptors,'growth')
print len(diurnalFilteredDescriptors),'filtered diurnal descriptors.'
print len(growthFilteredDescriptors),'filtered growth descriptors.'
print

# 1.3. plotting a boxplots of the best descriptors
print 'computing the border values for the descriptors...'
borders={}
borders=bordersCalculator(diurnalFilteredDescriptors,borders,'diurnal')
borders=bordersCalculator(growthFilteredDescriptors,borders,'growth')

# 1.4. saving the descriptors
print 'writing the descriptors...'
descriptorsWriter(diurnalFilteredDescriptors,'diurnal')
descriptorsWriter(growthFilteredDescriptors,'growth')

# 1.5. plotting the boxplots
if boxplotPlotting == True:
    print 'plotting expression graphs for descriptors...'
    boxPlotGrapher(diurnalFilteredDescriptors,borders,'diurnal')
    boxPlotGrapher(growthFilteredDescriptors,borders,'growth')

# 2. map samples into a new space of dark/light distributed in x:-2:-1/1:2 and stationary/exponential y:-2:-1/1:2
print
print 'mapping samples into new space...'
newSpaceMapper('300')
print
newSpaceMapper('1000')

# 3. generating plots of condition dispersion
print
print 'computing sample distances...'
distanceQuantifier('300')
distanceQuantifier('1000')


# 3. map samples into a new space in a probability manner
#print
#print 'mapping samples into a probabilistic space...'
#newSpaceProbabilisticMapper('300')
#sys.exit()

# 3. map samples into a 1D variable
#print 'mapping samples into 1D space...'
#oneDimensionTimeMapper()

# 4. final message
print
print '... analysis completed.'
print

### this script locates all epochs samples into a diurnal vs. growth space characterized in epoch=0, using microarray data

import sys,numpy,matplotlib,random,pickle
from matplotlib import pyplot
from matplotlib.patches import Ellipse

def bordersCalculator(descriptors,borders,flag):

    '''
    this function computes the borders of the descriptors in log10 scale
    '''

    # 0. defining a variable for recapitulating the data borders necessary for mapping samples into new space
    borders[flag]={}
    
    # 1. defining the expression values
    if flag == 'diurnal':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,'light')
        expressionB=descriptorsExpressionRetriever(descriptors,flag,'dark')           
    elif flag == 'growth':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,[1,2])
        expressionB=descriptorsExpressionRetriever(descriptors,flag,[4,5])
    else:
        print(flag)
        print('error handling flags from boxPlotGrapher')
        sys.exit()

    # 3. plotting and recovering the info for the mapping samples
    listOfDescriptors=sorted(descriptors,key=descriptors.__getitem__,reverse=True) # ranking the descriptors
    for geneID in listOfDescriptors:
        borders[flag][geneID]={}
        
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        xa=numpy.min(x); xb=numpy.median(x); xc=numpy.max(x); sdx=numpy.std(x)
        ya=numpy.min(y); yb=numpy.median(y); yc=numpy.max(y); sdy=numpy.std(y)
        if xb > yb:
            center=((xa-yc)/2.)+yc
        else:
            center=((ya-xc)/2.)+xc
        # incorporating the data
        w=weightNMLCalculator(geneID,flag)
        if flag == 'diurnal':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        elif flag == 'growth':
            borders[flag][geneID]=[xa,xb,xc,ya,yb,yc,center,w,sdx,sdy]
        else:
            print('error handling flags from boxPlotGrapher (bis)')
            sys.exit()
        
    return borders

def boxPlotGrapher(descriptors,borders,flag):

    '''
    this function plots the expression of the descriptors in log10 scale using boxplots
    '''

    # 1. defining the expression values
    if flag == 'diurnal':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,'light')
        expressionB=descriptorsExpressionRetriever(descriptors,flag,'dark')           
    elif flag == 'growth':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,[1,2])
        expressionB=descriptorsExpressionRetriever(descriptors,flag,[4,5])
    else:
        print(flag)
        print('error handling flags from boxPlotGrapher')
        sys.exit()

    # 3. plotting and recovering the info for the mapping samples
    boxPlotPosition=0
    names=[]
    listOfDescriptors=sorted(descriptors,key=descriptors.__getitem__,reverse=True) # ranking the descriptors
    for geneID in listOfDescriptors:
        boxPlotPosition=boxPlotPosition+1
        
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        # 3.2. actual plotting
        bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
        if flag == 'diurnal':
            setBoxColors(bp,'orange')
        else:
            setBoxColors(bp,'#0571b0')
        bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
        if flag == 'diurnal':
            setBoxColors(bp,'darkgreen')
        else:
            setBoxColors(bp,'#ca0020')

        name=geneID.split('_')[1]
        names.append(name)
        
    # 3.3. closing the figure
    matplotlib.pyplot.xlim([0,boxPlotPosition+1])
    matplotlib.pyplot.ylim([-5.,5.])
    theXticks=range(boxPlotPosition)
    theXticksPosition=[element+1 for element in theXticks]
    theFontSize=int(600./len(listOfDescriptors))
    matplotlib.pyplot.xticks(theXticksPosition,names,rotation=90,fontsize=theFontSize)
    matplotlib.pyplot.ylabel('log2 fold')
    matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off')
    matplotlib.pyplot.savefig('figuresGSE/boxplots_%s.pdf'%flag)
    matplotlib.pyplot.clf()

    # 4. computing the figures for specific AM/exp or PM/sta descriptors

    # AM/exp
    boxPlotPosition=0
    names=[]
    borA=[]
    borB=[]
    for geneID in listOfDescriptors:
        
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        xa=numpy.min(x); xb=numpy.median(x); xc=numpy.max(x); sdx=numpy.std(x)
        ya=numpy.min(y); yb=numpy.median(y); yc=numpy.max(y); sdy=numpy.std(y)

        if xb > yb:
            boxPlotPosition=boxPlotPosition+1
            borA.append(xa); borB.append(yc)
            
            bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
            if flag == 'diurnal':
                setBoxColors(bp,'orange')
            else:
                setBoxColors(bp,'#0571b0')
            bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
            if flag == 'diurnal':
                setBoxColors(bp,'darkgreen')
            else:
                setBoxColors(bp,'#ca0020')

            name=geneID.split('_')[1]
            names.append(name)

    # closing the figure
    matplotlib.pyplot.fill_between(range(1,len(names)+1),borA,borB,facecolor='black',alpha=0.2,edgecolor='None')
    
    matplotlib.pyplot.xlim([0,boxPlotPosition+1])
    matplotlib.pyplot.ylim([-5.,5.])
    theXticks=range(boxPlotPosition)
    theXticksPosition=[element+1 for element in theXticks]
    theFontSize=int(600./len(listOfDescriptors))
    if flag == 'diurnal':
        matplotlib.pyplot.plot([-1],[-1],color='orange',lw=2,label='light')
        matplotlib.pyplot.plot([-1],[-1],color='darkgreen',lw=2,label='dark')
        matplotlib.pyplot.legend()
    else:
        matplotlib.pyplot.plot([-1],[-1],color='#0571b0',lw=2,label='early')
        matplotlib.pyplot.plot([-1],[-1],color='#ca0020',lw=2,label='late')
        matplotlib.pyplot.legend()
    matplotlib.pyplot.xticks(theXticksPosition,names,rotation=90,fontsize=14)
    matplotlib.pyplot.ylabel('log2 fold')
    matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off')
    matplotlib.pyplot.savefig('figuresGSE/boxplots_%s.trend.up.pdf'%flag)
    matplotlib.pyplot.clf()

    # PM/sta
    boxPlotPosition=0
    names=[]
    borA=[]
    borB=[]
    for geneID in listOfDescriptors:
        
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        xa=numpy.min(x); xb=numpy.median(x); xc=numpy.max(x); sdx=numpy.std(x)
        ya=numpy.min(y); yb=numpy.median(y); yc=numpy.max(y); sdy=numpy.std(y)

        if xb < yb:
            boxPlotPosition=boxPlotPosition+1
            borA.append(xc); borB.append(ya)
            
            bp=matplotlib.pyplot.boxplot([x],positions=[boxPlotPosition],patch_artist=True)
            if flag == 'diurnal':
                setBoxColors(bp,'orange')
            else:
                setBoxColors(bp,'#0571b0')
            bp=matplotlib.pyplot.boxplot([y],positions=[boxPlotPosition],patch_artist=True)
            if flag == 'diurnal':
                setBoxColors(bp,'darkgreen')
            else:
                setBoxColors(bp,'#ca0020')

            name=geneID.split('_')[1]
            names.append(name)

    # closing the figure
    matplotlib.pyplot.fill_between(range(1,len(names)+1),borA,borB,facecolor='gray',alpha=0.2,edgecolor='None')
    
    matplotlib.pyplot.xlim([0,boxPlotPosition+1])
    matplotlib.pyplot.ylim([-5.,5.])
    theXticks=range(boxPlotPosition)
    theXticksPosition=[element+1 for element in theXticks]
    theFontSize=int(600./len(listOfDescriptors))
    if flag == 'diurnal':
        matplotlib.pyplot.plot([-1],[-1],color='orange',lw=2,label='light')
        matplotlib.pyplot.plot([-1],[-1],color='darkgreen',lw=2,label='dark')
        matplotlib.pyplot.legend()
    else:
        matplotlib.pyplot.plot([-1],[-1],color='#0571b0',lw=2,label='early')
        matplotlib.pyplot.plot([-1],[-1],color='#ca0020',lw=2,label='late')
        matplotlib.pyplot.legend()
    matplotlib.pyplot.xticks(theXticksPosition,names,rotation=90,fontsize=14)
    matplotlib.pyplot.ylabel('log2 fold')
    matplotlib.pyplot.tight_layout(pad=2.5)
    matplotlib.pyplot.tick_params(axis='x',which='both',top='off')
    matplotlib.pyplot.tick_params(axis='y',which='both',right='off')
    matplotlib.pyplot.savefig('figuresGSE/boxplots_%s.trend.down.pdf'%flag)
    matplotlib.pyplot.clf()

    return None

def descriptorsFilter(descriptors,flag):

    '''
    this function filters out descriptors
    '''

    # 1. defining the expression values
    if flag == 'diurnal':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,'light')
        expressionB=descriptorsExpressionRetriever(descriptors,flag,'dark')        
    elif flag == 'growth':
        expressionA=descriptorsExpressionRetriever(descriptors,flag,[1,2])
        expressionB=descriptorsExpressionRetriever(descriptors,flag,[4,5])
    else:
        print('error handling flags from boxPlotGrapher')
        sys.exit()

    # 2. computing the separations
    separations={}
    for geneID in descriptors: 
        
        # 2.1 transforming into log10 scale
        x=numpy.array(expressionA[geneID])
        y=numpy.array(expressionB[geneID])

        # computing a relative separation: it should be larger than the average of the standard deviations of the two distributions
        averageSD=0.5*numpy.std(x) + 0.5*numpy.std(y)

        # defining the separation
        if numpy.mean(x) > numpy.mean(y):
            separation=numpy.min(x)-numpy.max(y)
        elif numpy.mean(y) > numpy.mean(x):
            separation=numpy.min(y)-numpy.max(x)
        else:
            print('error computing the separation between the two distributions from descriptorsFilter')
            sys.exit()
        if separation > averageSD:
            separations[geneID]=separation

    return separations

def descriptorsRetriever(flag):

    '''
    this function matches the previously found descriptors with the array data
    '''

    jar=flag+'.pckl'
    f=open(jar,'r')
    fpkmDescriptors=pickle.load(f)
    f.close()
    print(len(fpkmDescriptors),flag,'FPKM descriptors recovered.')

    listOfDescriptors=[]
    for descriptor in fpkmDescriptors:
        for putative in expression.keys():
            if 'draft' in descriptor:
                a=int(descriptor.split('gene:Thapsdraft')[1])
            else:
                a=int(descriptor.split('gene:Thaps')[1])
            label='THAPSDRAFT_'
            if label in putative:
                b=int(putative.split('THAPSDRAFT_')[1])
            else:
                b=''
            if a == b:
                listOfDescriptors.append(putative)

    return listOfDescriptors

def expressionReader():

    '''
    this function reads the matrix of expression in FPKM into the format of a dictionary
    '''

    expression={}
    with open(expressionFile,'r') as f:
        header=f.readline()
        prelabels=header.split('\t')
        labels=[element.split('_')[0] for element in prelabels[11:30]]
        for line in f:
            vector=line.split('\t')
            geneName=vector[1]
            expression[geneName]={}
            
            preValues=vector[11:30]
            values=[float(element) for element in preValues]
            for i in range(len(values)):
                expression[geneName][labels[i]]=values[i]

    return expression

def descriptorsExpressionRetriever(descriptors,flag,condition):

    '''
    this function returns the expression values for a set of genes under a specific condition
    '''

    # 1. selecting the acceptable samples
    selectedSamples=[]
    for sampleID in metaData.keys():
        if metaData[sampleID]['co2'] == 400:
            if condition == 'light' or condition == 'dark':
                if metaData[sampleID][flag] == condition:
                    selectedSamples.append(sampleID)
            if len(condition) == 2:
                if metaData[sampleID][flag] in condition:
                    selectedSamples.append(sampleID)

    print(flag,condition,selectedSamples)
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
    
    for descriptor in borders[flag].keys():
        s=expression[descriptor][sampleID] # sample expression level

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

    return averageLoad

def metadataReader():

    '''
    this function returns a dictionary with the metadata of the expression data
    '''
    
    metaData={}

    with open(expressionFile,'r') as f:
        header=f.readline()
        vector=header.split('\t')

        for sampleID in vector:
            if 'X' in sampleID:

                # co2 levels
                co2=int(sampleID.split('.')[0].replace('X',''))

                # diurnal
                if sampleID.split('.')[-1] == 'Dk':
                    diurnal='dark'
                elif sampleID.split('.')[-1] == 'Lt':
                    diurnal='light'

                # day
                growth=int(sampleID.split('.')[1].replace('Day',''))

                # filling up
                metaData[sampleID]={}
                metaData[sampleID]['co2']=co2
                metaData[sampleID]['diurnal']=diurnal
                metaData[sampleID]['growth']=growth

    return metaData

def newCoordinateCalculator(sampleID):

    '''
    this function computes the new coordinates of a sample given the borders calculated from the filtered descriptors
    '''

    x=-loadCalculator(sampleID,'diurnal')
    y=loadCalculator(sampleID,'growth')
    
    return x,y

def newSpaceMapper(flag):

    '''
    this function plots the samples into a new space
    '''

    preselectedSamples=[sampleID for sampleID in metaData.keys() if metaData[sampleID]['co2'] == flag]
    print('selected ', len(preselectedSamples), 'samples for plotting.')

    # starting the figure
    fig=matplotlib.pyplot.figure()
    ax=fig.add_subplot(111)

    for sampleID in preselectedSamples:
        
        x,y=newCoordinateCalculator(sampleID)
        print(sampleID,x,y)
        
        theSize=16
        theAlpha=.85
        theMarker='o'

        # defining the color
        if metaData[sampleID]['diurnal'] == 'light':
            theColor='orange'
        elif metaData[sampleID]['diurnal'] == 'dark':
            theColor='green'
        else:
            print('error selecting the color from newSpaceMapper.')
            sys.exit()

        theDay=sampleID.split('.')[1].replace('Day','')

        ax.plot(x,y,marker=theMarker,mew=2.,mec=theColor,mfc='None',ms=theSize,alpha=theAlpha,zorder=10)
        ax.plot(x,y,marker="$%s$"%theDay,zorder=20)

    # setting ranges
    matplotlib.pyplot.xlim([-2.,2.])
    matplotlib.pyplot.ylim([-2.,2.])
    matplotlib.pyplot.xticks([-1.5,1.5],['light','dark'],fontsize=24)
    matplotlib.pyplot.yticks([-1.5,1.5],['late','early'],fontsize=24)
    matplotlib.pyplot.xlabel('diurnal cycle',fontsize=24)
    matplotlib.pyplot.ylabel('growth phase',fontsize=24)
    matplotlib.pyplot.title(str(flag)+' ppm',fontsize=28)

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

    matplotlib.pyplot.savefig('figuresGSE/sampleLocation.%s.pdf'%str(flag))
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

def weightNMLCalculator(geneID,flag):

    '''
    this function computes the weight of the descriptor based on the empty space between the descriptors
    '''

    if flag == 'diurnal':
        sumSpaces=sum(diurnalFilteredDescriptors.values())
        value=diurnalFilteredDescriptors[geneID]
    elif flag == 'growth':
        sumSpaces=sum(growthFilteredDescriptors.values())
        value=growthFilteredDescriptors[geneID]
    else:
        print('error from weightNMLCalculator')
        sys.exit()

    weight=value/sumSpaces

    return weight

### MAIN

# use the points defined in day 1 and use the gap as unit of measurement. rank equal weights.

# 0. preliminaries
print ('')
print('welcome to GSE_Mapper')
print('')
print('initializing variables...')

# 0.1. user defined variables and paths
expressionFile='../data/GSE45252_value_Matrix.txt'
boxplotPlotting=True

# 0.2. reading metadata
metaData=metadataReader()

# 0.3. reading expression
expression=expressionReader()

# 1. define descriptors

# 1.1. recover descriptors
print('')
print('recovering descriptors...')

diurnalDescriptors=descriptorsRetriever('diurnal')
growthDescriptors=descriptorsRetriever('growth')
print('')
print(len(diurnalDescriptors),'diurnal descriptors found after annotation check.')
print(len(growthDescriptors),'growth descriptors found after annotation check.')

# 1.2. filtering the descriptors based on separation
print('')
print('filtering descriptors based on separation...')

diurnalFilteredDescriptors=descriptorsFilter(diurnalDescriptors,'diurnal')
growthFilteredDescriptors=descriptorsFilter(growthDescriptors,'growth')
print(len(diurnalFilteredDescriptors),'filtered diurnal descriptors.')
print(len(growthFilteredDescriptors),'filtered growth descriptors.')
print('')

# 1.3. plotting a boxplots of the best descriptors
print('plotting expression graphs for descriptors...')
borders={}
borders=bordersCalculator(diurnalFilteredDescriptors,borders,'diurnal')
borders=bordersCalculator(growthFilteredDescriptors,borders,'growth')

# 1.4. plotting the boxplots
if boxplotPlotting == True:
    print('plotting expression graphs for descriptors...')
    boxPlotGrapher(diurnalFilteredDescriptors,borders,'diurnal')
    boxPlotGrapher(growthFilteredDescriptors,borders,'growth')

# 2. map samples into a new space of dark/light distributed in x:-2:-1/1:2 and stationary/exponential y:-2:-1/1:2
print('')
print('mapping samples into new space...')
newSpaceMapper(400)
print('')
newSpaceMapper(800)

# 4. final message
print('')
print('... analysis completed.')
print('')

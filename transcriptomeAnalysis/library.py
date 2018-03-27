def expressionReader(expressionFile):

    '''
    This function reads the matrix of expression in FPKM into the format of a dictionary.
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

def metadataReader(metaDataFile):

    '''
    This function returns a dictionary with the metadata of the expression data.
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

def sampleOrderer(co2level,metadata):

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
                for sampleID in metadata.keys():
                    if metadata[sampleID]['co2'] == co2level and metadata[sampleID]['epoch'] == epoch and metadata[sampleID]['growth'] == growth and metadata[sampleID]['diurnal'] == diurnal:

                        timeLabel=diurnal+'.'+growth+'.'+str(epoch+1)
                        if timeLabel not in timeLabels:
                            timeLabels.append(timeLabel)
                            
                        if time in orderedSamples:
                            orderedSamples[time].append(sampleID)
                        else:
                            orderedSamples[time]=[sampleID]

    return orderedSamples

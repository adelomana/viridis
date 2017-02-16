import sys,os

def cuffdiffCaller(bamFilesA,bamFilesB,label):

    '''
    this function calls cuffdiff to compute statistical test about expression differences
    '''
    
    outputDir=cuffdiffDir+label+'/'
        
    term1='time cuffdiff %s '%(gtfFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='--library-type fr-firststrand '
    term5='--multi-read-correct '
    term6='-b %s '%fastaFile
    
    term8a=','.join(bamFilesA)
    term8b=','.join(bamFilesB)
    term8=term8a+' '+term8b

    cmd=term1+term2+term3+term4+term5+term6+term8

    print()
    print(cmd)
    print()

    os.system(cmd)

    return None

def cuffnormCaller():

    '''
    this function calls cuffnorm using all generated abundance binary files by cuffquant
    '''

    outputDir=cufflinksDir+'allSamples'

    term1='cuffnorm %s '%(gtfFile)
    term2=' '.join(abundanceFiles)+' '
    term3='-o %s '%outputDir
    term4='-p %s '%numberOfThreads
    term5='--library-type fr-firststrand '
    term6='-L '+','.join(labels)+' '
    cmd=term1+term2+term3+term4+term5+term6
    
    print()
    print(cmd)
    print()

    os.system(cmd)

    return None

def cuffquantCaller(inputFile):

    '''
    this function calls the different steps of the cufflinks pipeline.
    '''

    label=inputFile.split('/')[-2]
    outputDir=cufflinksDir+label
    
    # cuffquant
    term1='cuffquant %s '%(gtfFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='-M %s '%maskFile
    term5='--library-type fr-firststrand '
    term6='--multi-read-correct '
    cmd=term1+term2+term3+term4+term5+term6+inputFile

    print()
    print(cmd)
    print()

    os.system(cmd)

    return None

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

def samples2bamfiles(samplesA,samplesB):

    '''
    this function builds the full path of the BAM files
    '''

    bamFilesA=[]
    for sample in samplesA:
        matching = [s for s in fullPathBamFiles if sample in s]
        fullName=matching[0]
        bamFilesA.append(fullName)

    bamFilesB=[]
    for sample in samplesB:
        matching = [s for s in fullPathBamFiles if sample in s]
        fullName=matching[0]
        bamFilesB.append(fullName)

    return bamFilesA,bamFilesB

# 0. defining input files
metaDataFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cufflinksDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/'
cuffdiffDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
fastaFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.dna.genome.fa'
maskFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/mask.gff3'
numberOfThreads=12

# 1. reading the metadata
metaData=metadataReader()

# 2. defining the BAM and abundance files
roots=os.listdir(bamFilesDir)
roots.remove('secondRun')

fullPathBamFiles=[bamFilesDir+element+'/Aligned.sortedByCoord.out.bam' for element in roots]
abundanceFiles=[cufflinksDir+element+'/abundances.cxb' for element in roots]
labels=[element.split('_')[-1] for element in roots]

# 2. calling cuffquantCaller 
print 'calling cuffquant...'
for inputFile in bamFiles:
    cuffquantCaller(inputFile)

# 3. calling cuffnorm
print 'calling cuffnorm...'
cuffnormCaller()

# 4. running cuffdiff
print('calling cuffdiff...')

# 4.1. running cuffdiff late versus early samples (only in the light) (1 test)
testLabel='early.vs.late-AM'
samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['growth'] == 'exp' and metaData[sampleID]['diurnal'] == 'AM':
        samplesA.append(sampleID)
    elif metaData[sampleID]['growth'] == 'sta' and metaData[sampleID]['diurnal'] == 'AM':
        samplesB.append(sampleID)

print(samplesA,len(samplesA))
print(samplesB,len(samplesB))
    
bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print(len(bamFilesA),len(bamFilesB))

cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

# 4.2. running cuffdiff late versus early samples (only the light), in two separate carbon levels (2 tests)
testingConditions=[300,1000]
for testingCondition in testingConditions:

    testLabel='early.vs.late-AM.{}'.format(testingCondition)
    samplesA=[]
    samplesB=[]

    for sampleID in metaData.keys():
        if metaData[sampleID]['growth'] == 'exp' and metaData[sampleID]['diurnal'] == 'AM' and metaData[sampleID]['co2'] == testingCondition:
            samplesA.append(sampleID)
        elif metaData[sampleID]['growth'] == 'sta' and metaData[sampleID]['diurnal'] == 'AM' and metaData[sampleID]['co2'] == testingCondition:
            samplesB.append(sampleID)

    print(samplesA,len(samplesA))
    print(samplesB,len(samplesB))

    bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

    print(len(bamFilesA),len(bamFilesB))

    cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

# 4.3. running cuffdiff late versus early samples (both in the light), in separate carbon and epochs (10 tests)
testingConditionsA=[300,1000]
testingConditionsB=[0,1,2]
testingConditionsC=['AM','PM']
for testingConditionA in testingConditionsA:
    for testingConditionB in testingConditionsB:
        for testingConditionC in testingConditionsC:

            print(testingConditionA,testingConditionB,testingConditionC)
            
            testLabel='early.vs.late-{}.{}.{}'.format(testingConditionA,testingConditionB,testingConditionC)
            print('### {}'.format(testLabel))
            
            samplesA=[]
            samplesB=[]

            for sampleID in metaData.keys():
                if metaData[sampleID]['growth'] == 'exp' and metaData[sampleID]['co2'] == testingConditionA and metaData[sampleID]['epoch'] == testingConditionB and metaData[sampleID]['diurnal'] == testingConditionC:
                    samplesA.append(sampleID)
                elif metaData[sampleID]['growth'] == 'sta' and metaData[sampleID]['co2'] == testingConditionA and metaData[sampleID]['epoch'] == testingConditionB and metaData[sampleID]['diurnal'] == testingConditionC:
                    samplesB.append(sampleID)

            print(samplesA,len(samplesA))
            print(samplesB,len(samplesB))

            bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

            print(len(bamFilesA),len(bamFilesB))

            cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

# 4.4. running cuffdiff HC vs LC, excluding last epoch (8 tests)
testingConditionsA=[0,1,2]
testingConditionsB=['exp','sta']
testingConditionsC=['AM','PM']
for testingConditionA in testingConditionsA:
    for testingConditionB in testingConditionsB:
        for testingConditionC in testingConditionsC:

            print(testingConditionA,testingConditionB,testingConditionC)
            
            testLabel='LC.vs.HC-{}.{}.{}'.format(testingConditionA,testingConditionB,testingConditionC)
            print('### {}'.format(testLabel))
            
            samplesA=[]
            samplesB=[]

            for sampleID in metaData.keys():
                if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] == testingConditionA and metaData[sampleID]['growth'] == testingConditionB and metaData[sampleID]['diurnal'] == testingConditionC:
                    samplesA.append(sampleID)
                elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] == testingConditionA and metaData[sampleID]['growth'] == testingConditionB and metaData[sampleID]['diurnal'] == testingConditionC:
                    samplesB.append(sampleID)

            print(samplesA,len(samplesA))
            print(samplesB,len(samplesB))

            bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

            print(len(bamFilesA),len(bamFilesB))

            cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

# 4.5. running cuffdiff HC vs LC, all combined, or first epoch
testLabel='LC.vs.HC'
samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300:
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000:
        samplesB.append(sampleID)

print(samplesA,len(samplesA))
print(samplesB,len(samplesB))

bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print(len(bamFilesA),len(bamFilesB))

cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

###

testLabel='LC.vs.HC.only.epoch.0'
samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] == 0:
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] == 0:
        samplesB.append(sampleID)

print(samplesA,len(samplesA))
print(samplesB,len(samplesB))

bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print(len(bamFilesA),len(bamFilesB))

cuffdiffCaller(bamFilesA,bamFilesB,testLabel)

# 5. final message
print('... all done.')

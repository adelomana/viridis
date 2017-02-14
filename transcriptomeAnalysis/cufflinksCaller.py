import sys,os

def cuffdiffCaller(bamFilesA,bamFilesB,label):

    '''
    this function calls cuffdiff to compute statistical test about expression differences
    '''
    
    outputDir=cuffdiffDir+label+'/'
        
    term1='cuffdiff %s '%(gtfFile)
    term2='-o %s '%outputDir
    term3='-p %s '%numberOfThreads
    term4='--library-type fr-firststrand '
    term5='--multi-read-correct '
    term6='-b %s '%fastaFile
    
    term8a=','.join(bamFilesA)
    term8b=','.join(bamFilesB)
    term8=term8a+' '+term8b

    cmd=term1+term2+term3+term4+term5+term6+term8

    print
    print cmd
    print

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

# 0. defining input files
metaDataFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cufflinksDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/'
cuffdiffDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
maskFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/mask.gff3'
numberOfThreads=12

# 1. reading the metadata
metaData=metadataReader()

# 2. defining the BAM and abundance files
roots=os.listdir(bamFilesDir)
roots.remove('secondRun')

bamFiles=[bamFilesDir+element+'/Aligned.sortedByCoord.out.bam' for element in roots]
abundanceFiles=[cufflinksDir+element+'/abundances.cxb' for element in roots]
labels=[element.split('_')[-1] for element in roots]

# 2. calling cuffquantCaller 
#print 'calling cuffquant...'
#for inputFile in bamFiles:
#    cuffquantCaller(inputFile)

# 3. calling cuffnorm
#print 'calling cuffnorm...'
#cuffnormCaller()

# 4. running cuffdiff
print('calling cuffdiff...')

# 4.1. running cuffdiff late versus early samples (only in the light) (1 test)
testLabel='early.vs.late.inlight'
samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['growth'] == 'early':
        samplesA.append(sampleID)
    elif metaData[sampleID]['growth'] == 'late':
        samplesB.append(sampleID)
    else:
        print('error selecting samples. exiting...')
        sys.exit()

cuffdiffCaller(samplesA,samplesB)

# 4.2. running cuffdiff late versus early samples (only the light), in two separate carbon levels (2 tests)

# # 4.3. running cuffdiff late versus early samples (only in the light), in separate carbon and epochs (10 tests)

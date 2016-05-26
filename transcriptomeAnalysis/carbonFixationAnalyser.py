### this script tests for the effect of C levels to specific genes associated with CCM and CA (carbonic anhydrase)

import os,sys

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

def samples2bamfiles(samplesA,samplesB):

    '''
    this function builds the full path of the BAM files
    '''

    bamFilesA=[]
    for sample in samplesA:
        matching = [s for s in verboseBamDirs if sample in s]
        fullName=matching[0]
        bamFilesA.append(bamFilesDir+fullName+'/Aligned.sortedByCoord.out.bam')

    bamFilesB=[]
    for sample in samplesB:
        matching = [s for s in verboseBamDirs if sample in s]
        fullName=matching[0]
        bamFilesB.append(bamFilesDir+fullName+'/Aligned.sortedByCoord.out.bam')

    return bamFilesA,bamFilesB

def tableFormatting(targetTranscripts,hypoLabel,functionLabel,n):

    '''
    this function builds a table with the formatted results from the hypothesis tests
    '''

    accCount=0

    outputFile=cuffdiffDir+hypoLabel+'/formattedInfo.%s.txt'%(functionLabel)
    g=open(outputFile,'w')

    header=['GeneID','LC (FPKM), n = %s'%n,'HC (FPKM), n = %s'%n,'log2 FC','p-value','significance','accCount']
    headerString='\t'.join(header)
    g.write(headerString)
    g.write('\n')

    inputFile=cuffdiffDir+'%s/isoform_exp.diff'%hypoLabel
    with open(inputFile,'r') as f:
        f.next()
        for line in f:
            vector=line.split('\t')

            geneName=vector[0].split(':')[1]
            if geneName in targetTranscripts:
                
                expressionLC=vector[7]
                expressionHC=vector[8]
                foldChange=vector[9]
                pValue=vector[11]
                significance=vector[13].replace('\n','')
                
                if significance == 'yes':
                    accCount=accCount+1

                accCountString='%s/%s'%(accCount,len(targetTranscripts))

                info=[geneName,expressionLC,expressionHC,foldChange,pValue,significance,accCountString]
                infoLine='\t'.join(info)
                g.write(infoLine)
                g.write('\n')

    g.close()

    return None
    

# 0. user defined variables
metaDataFile='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/metadata/metadata.v2.tsv'
bamFilesDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/bamFiles/'
cufflinksDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cufflinks/'
cuffdiffDir='/proj/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/cuffdiff/'
gtfFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.gff3'
fastaFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/Thalassiosira_pseudonana.ASM14940v1.29.dna.genome.fa'
maskFile='/proj/omics4tb/alomana/projects/dtp/data/ensembl/mask.gff3'
numberOfThreads=4 

ccmTranscripts=['Thaps10360','Thaps14147','Thaps2078','Thaps22208','Thaps233','Thaps25840','Thaps260953','Thaps263924','Thaps264181','Thaps269942','Thaps3353','Thaps34030','Thaps34125','Thaps34585','Thaps36208','Thaps39799','Thaps4819','Thaps4820','Thaps6528','Thaps6529','Thaps9903']
caTranscripts=['Thaps814','Thaps34125','Thaps233','Thaps25840','Thaps262009','Thaps34094','Thaps22391','Thaps41566','Thaps31046','Thaps22596','Thaps31215','Thaps5284','Thaps262006']

# 1. reading the metadata
metaData=metadataReader()
verboseBamDirs=os.listdir(bamFilesDir)

# 2. hypo1: is CCM expression different from LC to HC samples?
print
print 'about to call cuffdiff for all conditions, (24/24) conditions expected...'

samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] != 2:
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] != 2:
        samplesB.append(sampleID)

print samplesA
print samplesB

bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print len(bamFilesA),len(bamFilesB)

cuffdiffCaller(bamFilesA,bamFilesB,'hypo1')

tableFormatting(ccmTranscripts,'hypo1','ccm',len(bamFilesA))
tableFormatting(caTranscripts,'hypo1','ca',len(bamFilesA))

# 3. hypo2: is CCM expression different from LC to HC specifically at early light conditions?
print
print 'about to call cuffdiff for AM early conditions, (6/6) conditions expected...'

samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['growth'] == 'exp' and metaData[sampleID]['diurnal'] == 'AM':
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['growth'] == 'exp' and metaData[sampleID]['diurnal'] == 'AM':
        samplesB.append(sampleID)

print samplesA
print samplesB
    
bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print len(bamFilesA),len(bamFilesB)

cuffdiffCaller(bamFilesA,bamFilesB,'hypo2')

tableFormatting(ccmTranscripts,'hypo2','ccm',len(bamFilesA))
tableFormatting(caTranscripts,'hypo2','ca',len(bamFilesA))

# 4. hypo3: is CCM expression different from LC to HC specifically light conditions?
print
print 'about to call cuffdiff for AM conditions, (12/12) conditions expected...'

samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['diurnal'] == 'AM':
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['diurnal'] == 'AM':
        samplesB.append(sampleID)

print samplesA
print samplesB
    
bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)

print len(bamFilesA),len(bamFilesB)

cuffdiffCaller(bamFilesA,bamFilesB,'hypo3')

tableFormatting(ccmTranscripts,'hypo3','ccm',len(bamFilesA))
tableFormatting(caTranscripts,'hypo3','ca',len(bamFilesA))

# 5. hypo4: is CCM expression different from LC to HC specifically light late conditions?
print
print 'about to call cuffdiff for AM late conditions, (6/6) conditions expected...'

samplesA=[]
samplesB=[]

for sampleID in metaData.keys():
    if metaData[sampleID]['co2'] == 300 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['diurnal'] == 'AM' and metaData[sampleID]['growth'] == 'sta':
        samplesA.append(sampleID)
    elif metaData[sampleID]['co2'] == 1000 and metaData[sampleID]['epoch'] != 2 and metaData[sampleID]['diurnal'] == 'AM' and metaData[sampleID]['growth'] == 'sta':
        samplesB.append(sampleID)

print samplesA
print samplesB

bamFilesA,bamFilesB=samples2bamfiles(samplesA,samplesB)
    
print len(bamFilesA),len(bamFilesB)

cuffdiffCaller(bamFilesA,bamFilesB,'hypo4')

tableFormatting(ccmTranscripts,'hypo4','ccm',len(bamFilesA))
tableFormatting(caTranscripts,'hypo4','ca',len(bamFilesA))

# 6. final message
print '... all done!'

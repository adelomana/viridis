### this script computes a pairwise Spearman correlation based on the intersect of genes that have more than tau counts.

import os,sys,scipy,matplotlib
from scipy import stats
from matplotlib import pyplot

matplotlib.rcParams.update({'xtick.labelsize':8,'ytick.labelsize':8})

def expressionReader(sample):

    expression={}
    inputFile=countsDir+sample
    with open(inputFile,'r') as f:
        for line in f:
            vector=line.split()
            if '__' not in vector[0]:
                expression[vector[0]]=int(vector[1])

    return expression

# 0. defining user variables
countsDir='/Volumes/omics4tb/alomana/projects/dtp/data/expression/tippingPoints/counts/'
tau=10

# 1. reading files
print('reading files...')
sampleLabels=[
    'A1000_6_27_2015A',
    'B1000_6_27_2015A',
    'C1000_6_27_2015A',
    'A1000_6_27_2015P',
    'B1000_6_27_2015P',	
    'C1000_6_27_2015P',
    'A1000_6_29_2015A',
    'B1000_6_29_2015A',
    'C1000_6_29_2015A',
    'A1000_6_29_2015P',
    'B1000_6_29_2015P',	
    'C1000_6_29_2015P',
    'A1000_7_01_2015A',	
    'B1000_7_01_2015A',	
    'C1000_7_01_2015A',	
    'A1000_7_01_2015P',	
    'B1000_7_01_2015P',	
    'C1000_7_01_2015P',		
    'A1000_7_03_2015A',	
    'B1000_7_03_2015A',	
    'C1000_7_03_2015A',		
    'A1000_7_03_2015P',	
    'B1000_7_03_2015P',	
    'C1000_7_03_2015P',		
    'A1000_7_11_2015A',	
    'C1000_7_11_2015A',		
    'A1000_7_11_2015P',	
    'C1000_7_11_2015P',	
    'A1000_7_13_2015A',	
    'C1000_7_13_2015A',		
    'A1000_7_13_2015P',
    'C1000_7_13_2015P',
    'D300_6_27_2015A',
    'E300_6_27_2015A',			
    'F300_6_27_2015A',	
    'D300_6_27_2015P',
    'E300_6_27_2015P',	
    'F300_6_27_2015P',		
    'D300_6_29_2015A',	
    'E300_6_29_2015A',
    'F300_6_29_2015A',		
    'D300_6_29_2015P',
    'E300_6_29_2015P',
    'F300_6_29_2015P',
    'D300_7_01_2015A',	
    'E300_7_01_2015A',
    'F300_7_01_2015A',		
    'D300_7_01_2015P',	
    'E300_7_01_2015P',
    'F300_7_01_2015P',	
    'D300_7_03_2015A',	
    'E300_7_03_2015A',
    'F300_7_03_2015A',	
    'D300_7_03_2015P',	
    'E300_7_03_2015P',			
    'F300_7_03_2015P'
    ]

"""
###
sampleLabels=[
   
    'A1000_7_11_2015A',	
    'C1000_7_11_2015A',		

    'A1000_7_11_2015P',	
    'C1000_7_11_2015P',	

    'A1000_7_13_2015A',	
    'C1000_7_13_2015A',
    	
    'A1000_7_13_2015P',
    'C1000_7_13_2015P'
]
###

###
sampleLabels=[
   
    'A1000_7_11_2015A',	
    'C1000_7_11_2015A',		

    'A1000_7_13_2015A',	
    'C1000_7_11_2015P',	

    'A1000_7_11_2015P',	
    'C1000_7_13_2015A',
    	
    'A1000_7_13_2015P',
    'C1000_7_13_2015P'
]
###
"""
    
samples=[]
allFiles=os.listdir(countsDir)
for label in sampleLabels:
    finalLabel=label.replace('_','-')
    for file in allFiles:
        if file.find(finalLabel) != -1 and file.find('1000') != -1:
            samples.append(file)

#samples=samples[:10]

print(samples)

# 2. computing the correlation coefficient
print('computing the correlation coefficients...')
M=[]
labels=[]
controlSamples=[]
for sample1 in samples:
    label=sample1.split('_')[1].split('.')[0]
    labels.append(label)
    controlSamples.append(sample1)
    expression1=expressionReader(sample1)
    trimmedExpression1={key:expression1[key] for key in expression1 if expression1[key]>tau}
    V=[]
    for sample2 in samples:
        expression2=expressionReader(sample2)
        trimmedExpression2={key:expression2[key] for key in expression2 if expression2[key]>tau}

        # finding the overlap
        genes1=trimmedExpression1.keys()
        genes2=trimmedExpression2.keys()
        intersect=list(set(genes1) & set(genes2))

        # computing the correlation
        x=[]
        y=[]
        for element in intersect:
            x.append(expression1[element])
            y.append(expression2[element])
        rho,pval=scipy.stats.spearmanr(x,y)
        print(sample1,sample2,rho)
        V.append(rho)
    M.append(V)
        
# 3. plotting the matrix of correlations
print('plotting the matrix...')
matplotlib.pyplot.imshow(M,interpolation='none',cmap='jet')
matplotlib.pyplot.colorbar()
matplotlib.pyplot.grid(False)
matplotlib.pyplot.xticks(range(len(labels)),labels,rotation=90)
matplotlib.pyplot.yticks(range(len(labels)),labels)
matplotlib.pyplot.tight_layout()
matplotlib.pyplot.savefig('correlations.counts.1000.png')

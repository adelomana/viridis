# viridis
Tools for quantitative analyses of diatom resilience.  

## Growth Analysis

library.py: library containing the common functions, including the fitting functions.  
timeSeriesFitter.py: script to fit each time series individually with the purpose of manually checking the goodness of fit.  
epochGrapher.py: script to generate the epoch figures for each condition.  
nicheBreadthGrapher.py: script to generate the niche breadth increase based on maximum growth.  
tippingPointTransitionGrapher.py: script to generate the tipping point transition based on growth lag and carrying capacity. 

## Transcriptome Analysis

The natural order for the analysis would follow as:   

#### cleaning reads
readsCleaner.py: script to call Trimmomatic and clean the reads.

#### mapping reads
readsMapper.py: script to call STAR and map the reads to the genome.

#### quantifying reads (counts)
readsCounter.py: script to call HTSeq to counts the reads per transcript.  

#### quantifying reads (FPKM)
cufflinksCaller.py: script to call cufflinks and quantify the mapped reads.  

#### finding differentially expressed transcripts  
classifiersFinder.py: script that calls cuffdiff to find DETs between two conditions.  
carbonFixationAnalyser.py: script that calls cuffdiff to find DETs between LC and HC for specific CCM transcripts.  

#### visualizing samples
samplesCorrelationGrapher.py: script to compute and plot the correlation among samples. Available both on counts and FPKM values.  
pcaGrapher.py: script to compute the PCA on expression values in FPKM.  
sampleMapper.py: script to map samples into diurnal and growth coordinates.  
GSE_Mapper.py: script to map samples into diurnal and growth space from microarray data.  

#### misc
logValuesMatrixCreator.py: script to transform absolute (FPKM) into relative expression (log2 fold change).
cartoonFigureGenerator.py: script to generate a cartoon about the state space calculation.

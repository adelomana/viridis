# viridis
Tools for quantifying gene expression in diatoms.  
    
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

#### visualizing samples
samplesCorrelationGrapher.py: script to compute and plot the correlation among samples.  
pcaGrapher.py: script to compute the PCA on expression values in FPKM.  
sampleMapper.py: script to map samples into diurnal and growth coordinates.  
GSE_Mapper.py: script to map samples into diurnal and growth space from microarray data.  

#### misc
logValuesMatrixCreator.py: script to transform absolute (FPKM) into relative expression (log2 fold change).


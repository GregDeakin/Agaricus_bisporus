#<i>Agaricus bisporus</i> microarray analysis

## Version 5 microarrays

###	Data handling
The data were imported into ‘R’ where normalisation, filtering and statistical analysis were undertaken using the Bioconductor Limma package and DESeq2. 
####	Normalisation
The data were background corrected using normal exponential convolution with an offset of 50. Between array normalisation was performed using the quantile method.
####	Data filtering
To reduce false discovery rate of differentially expressed genes, the data were filtered to remove features that showed no expression (using the Agilent Well Above Background Boolean flag). 
Each biological replicate set was filtered using an aggressive filter (≥ 50% present) due to the small number of biological replicates (<=four)  for each set. Features failing the filter across all sets were removed before calculating significance. Control probes were removed prior to normalisation. 
Outliers within a technical replicate set were replaced with the median value of that set. Outliers were defined by two criteria:

(1) log expression value greater/less than 1.0 from set median value and 
<br>(2) replicate set range greater than 2.0. 

Technical replicates were then combined and log expression value replaced with the set median value.
###	Statistical analysis
Comparisons were made between treated and control samples, with each time-point taken separately. Differences in transcript levels between the comparisons were calculated using Empirical Bayes Statistics for Differential Expression with Benjamini and Hochberg False Discovery Rate correction for multiple testing and DESeq2. Statistical significance was taken as p <= 0.05.
###	Go functional analysis
Go functional analysis was performed using the BiNGO plugin (2.44) for Cytoscape 2.8.3. Gene Ontology was downloaded from the Gene Ontology Consortium (geneontology.org/page/downloads). An <i>A. bisporus</i> gene to GO annotation file was created manually ((gene_association.GO_Agaricus) from publically available <i>A. bisporus</i> functional annotation 


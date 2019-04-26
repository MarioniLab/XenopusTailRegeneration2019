Instructions:

1. Download (and unzip) counts matrix data from array express https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7716/E-MTAB-7716.processed.1.zip 

2. Download gene lists from the additionalFiles github folder

3. Download code from github:
- scripts.R                       : collection of useful scripts 
- projectAndCluster.Rmd           : initial UMAP projection and clustering
- xenotailFigures.Rmd            : code to produce figures

See additionFiles/rpackages.txt for a list of installed R packages. UMAP must also be installed for code/projectAndCluster.Rmd. We use the python installation (https://umap-learn.readthedocs.io/en/latest/), and run this in R using reticulate (see code/scripts.R, also additionalFiles/environment.yaml).  



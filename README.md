# perturb.met
perturb-Met can identify transcriptional perturbations centered at hundreds of metabolites in a genome-scale metabolic network, using single cell or single nucleus RNA-seq data.  
perturb-Met can be used to identify metabolic changes between any two conditions, including using bulk transcriptomic data (RNA-seq or microarray).  
perturb-Met adapts the Kullback-Leibler (KL) divergence concept from information theory to quantify relative changes in the distribution of transcript abundances for all genes connected to a metabolite between two conditions (e.g., AD vs. control). KL divergence quantifies how a statistical distribution differs from a second, reference distribution.  

perturb-Met offers several advantages over existing methods. 
* First, perturb-Met addresses both problems by considering relative changes of transcript abundances at individual metabolites, thus it has higher resolution than pathway-level analysis.   
* Second, unlike the Reproter Metabolite analysis method, perturb-Met does not require a a minimum and a maximum number of genes connected to a metabolite (e.g., minimum 3 genes, maximum 50 genes). This avoid the problem of missing important metabolites connected by only two genes. 
* Third, even if more than 3 genes were connected to a metabolite, Reporter Metabolite analysis would still require enough of them to be significantly differentially expressed. A metabolite may not be selected if only one of its connected genes is differentially expressed. However, using  perturb-Met, it is still considered noteworthy if the magnitude of expression change for one gene significantly alters the relative distribution of transcript abundances of all genes around a metabolite.  


To install perturb-Met R package, you need to first install the _devltools_ R package, then run:
```r
library(devtools)
install_github("yuliangwang/perturb.met")
library(perturb.met)
```
You can check out how to use perturb-Met in the [demo on Rpubs](https://rpubs.com/wang341/perturb-met).   
Learn more about perturb-Met and its application on single cell RNA-seq of Alzheimer's disease and metabolic differences across neuronal subtypes in AD by reading my preprint: [Identifying neuron subtype-specific metabolic network changes in single cell transcriptomics of Alzheimer's Disease using perturb-Met](https://www.biorxiv.org/content/10.1101/2021.01.18.427154v1).


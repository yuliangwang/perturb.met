#'Identify perturbed metabolites based on changes in the relative expression of enzymes around a metabolite
#'This function does NOT perform metabolic network topology permutation. This is intended to be used when performing sample label permutation.
#'
#' @param expression_data data frame of the input expression data.First column must be gene symbol; second column is the expression level in reference condition; third comun is the expression level in treated/perturbed condition; expression in TPM or FPKM, not log transformed.
#' @param mets2genes data frame that describe gene-metabolite connections. It includes two columns: "met", which lists metabolite symbols; and "genes" which include gene symbols for genes that produce or consume the corresponding metabolite
#' @param expressed_threshold a value above which a gene is considered expressed
#' @return a list that has two components. 1. perturbed_mets: list of perturbed metabolites sorted by p-value; 2. met_gene_pairs: metabolite-gene pairs after filtering down to expressed genes and their connected metabolites;
#' @export
findPerturbMet_noPerm<-function(expression_data,mets2genes,expressed_threshold){
  expressed_genes<-as.character(expression_data[,1][expression_data[,2]>expressed_threshold|expression_data[,3]>expressed_threshold])

  met_gene_pairs<-mets2genes
  met_gene_pairs<-met_gene_pairs[met_gene_pairs$genes %in% expression_data[,1],]
  met_gene_pairs<-met_gene_pairs[-grep("^SLC",met_gene_pairs$genes),]

  gene_freq<-table(met_gene_pairs$genes)
  met_gene_pairs<-met_gene_pairs[met_gene_pairs$genes %in% names(gene_freq)[gene_freq<=15],] #because LDH is connected to 15 metabolites. Consider that connectivity only includes carbon substrates, no co-factors or water

  met_gene_pairs<-dplyr::arrange(met_gene_pairs,mets)
  temp<-met_gene_pairs[met_gene_pairs$genes %in% expressed_genes,]
  met_expressed_freq<-table(temp$names)

  met_freq<-table(met_gene_pairs$names)

  met_gene_pairs<-met_gene_pairs[(met_gene_pairs$names %in% names(met_freq)[met_freq>=2]) & (met_gene_pairs$names %in% names(met_expressed_freq)[met_expressed_freq>=1]),]
  met_sets<-unique(met_gene_pairs$mets)


  ind<-match(met_gene_pairs$genes,expression_data[,1])
  met_gene_pairs$exprs_ref<-expression_data[,2][ind]
  met_gene_pairs$exprs_treat<-expression_data[,3][ind]

  KL_met_set<-rep(0,length(met_sets))
  Dirichlet_met_set<-KL_met_set
  for (i in seq_along(met_sets)){
    sub_met_genes<-met_gene_pairs[met_gene_pairs$mets==met_sets[i],]
    KL_met_set[i]<-calculateKL_abs(sub_met_genes$exprs_ref,sub_met_genes$exprs_treat)
    Dirichlet_met_set[i]<-calculateKL_Dirichlet(sub_met_genes$exprs_ref,sub_met_genes$exprs_treat)
  }

  perturbed_mets<-data.frame(mets=met_sets,KL=KL_met_set,Dirichlet=Dirichlet_met_set)

}

#' Make a bar plot comparing the percentage changes of all enzymes connected to a metabolite between two conditions.
#' @param met_gene_exprs a dataframe that includes metabolite-gene connections and gene expression values.
#' @param met2look which metabolite to plot, using metabolite symbols
#' @param legend_coord where to place the legend. c(0,0) is bottom left; c(1,1) is top right.
#' @param n2plot how many genes to plot
#' @param ylimit y axis limit
#' @return Return a plot of bar plots

plot_pct_change<-function(met_gene_exprs,met2look,legend_coord,n2plot,ylimit){
  require(ggplot2)
  require(dplyr)
  sub_met_genes<-met_gene_exprs[met_gene_exprs$mets==met2look,]
  pseudo_count<-1e-1
  x<-sub_met_genes$exprs_ref+pseudo_count
  y<-sub_met_genes$exprs_treat+pseudo_count
  x<-x/sum(x)
  y<-y/sum(y)
  df<-data.frame(genes=rep(sub_met_genes$genes,2),pct=c(x,y),Conditions=rep(c("ref","treat"),each=length(x)))
  df$Conditions<-factor(df$Conditions,levels=c("ref","treat"),ordered=T)
  df<-arrange(df,Conditions)
  #plot only genes that changed much#
  sub_met_genes$exprs_ref<-round(sub_met_genes$exprs_ref,digits = 1)
  sub_met_genes$exprs_treat<-round(sub_met_genes$exprs_treat,digits = 1)
  sub_met_genes$pct_ref<-(sub_met_genes$exprs_ref+1e-6)/sum(sub_met_genes$exprs_ref+1e-6)
  sub_met_genes$pct_treat<-(sub_met_genes$exprs_treat+1e-6)/sum(sub_met_genes$exprs_treat+1e-6)
  sub_met_genes$contribution<- sub_met_genes$pct_treat*log2(sub_met_genes$pct_treat/sub_met_genes$pct_ref)
  sub_met_genes<-arrange(sub_met_genes,desc(abs(contribution)))

  df<-df[df$genes %in% sub_met_genes$genes[1:n2plot],]
  ggplot(data=df,aes(x=genes,y=pct)) + geom_bar(aes(fill=Conditions),stat='identity',position='dodge')+ scale_fill_manual(values=c('darkolivegreen3','darkorange')) + xlab("Genes") + ylab("Percentage") + scale_y_continuous(labels = scales::percent,limits=c(0,ylimit/100)) + theme(panel.background = element_rect(fill = "white",color='black'),panel.grid.major = element_line(color='grey'),legend.position = legend_coord,legend.justification = legend_coord,legend.title = element_blank())

}




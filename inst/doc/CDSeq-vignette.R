## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10, 
  fig.height = 7,
  fig.align = 'center'
)

## ----quick start, echo=TRUE, eval=FALSE---------------------------------------
#  result<-CDSeq(bulk_data =  mixtureGEP,
#                cell_type_number = 6,
#                mcmc_iterations = 1000,
#                cpu_number=1)

## ----installation, eval=FALSE-------------------------------------------------
#  install_github("kkang7/CDSeq_R_Package")

## ----library, results='hide', message=FALSE-----------------------------------
library(CDSeq)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(dplyr)

## ----example 1, echo=TRUE ,eval=FALSE, cache=FALSE, results='hide'------------
#  
#  result1<-CDSeq(bulk_data =  mixtureGEP,
#                 cell_type_number = 6,
#                 mcmc_iterations = 2000,
#                 gene_length = as.vector(gene_length),
#                 reference_gep = refGEP,  # gene expression profile of pure cell lines
#                 cpu_number = 1)

## ----output, eval=TRUE, echo=TRUE---------------------------------------------
ls(result1)

## ----fig1, echo=TRUE,cache=FALSE,eval=TRUE, fig.width=7,fig.height=5,fig.align = 'center'----
trueGEP <- true_GEP_rpkm[,result1$cell_type_assignment]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1), pty="s")
for (i in 1:6) {
  max_v <- max(c(result1$estGEP[,i], trueGEP[,i]))
  plot(result1$estGEP[,i],
       trueGEP[,i],
       xlab = "CDSeq-estimated GEPs", 
       ylab = "True GEPs",
       pch = 21,
       bg = 'green',
       cex=2, 
       xlim = c(0, max_v), 
       ylim = c(0,max_v),
       main = paste0("cell type ",i))
  lines(c(0,max_v),c(0,max_v), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----fig2, echo=TRUE,cache=FALSE,eval=TRUE,fig.width=7,fig.height=5,fig.align = 'center'----
trueProp <- true_prop_cell[result1$cell_type_assignment,]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1), pty="s")
for (i in 1:6) {
    plot(result1$estProp[i,],
         trueProp[i,],
         pch = 21,
         cex = 3,
         bg='red',
         xlab = "CDSeq-estimated proportions", 
         ylab = "True proportions", 
         xlim = c(0,1),
         ylim = c(0,1),
         main = paste0("cell type ",i))
    lines(c(0,1),c(0,1), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----time, echo=TRUE,cache=FALSE,eval=TRUE------------------------------------
result1$gibbsRunningTime

## ----example 2, echo=TRUE ,eval=FALSE, cache=FALSE, results='hide'------------
#  result2<-CDSeq(bulk_data =  mixtureGEP,
#                cell_type_number = 2:10,
#                mcmc_iterations = 2000,
#                dilution_factor = 1,
#                block_number = 1,
#                gene_length = as.vector(gene_length),
#                reference_gep = refGEP, # gene expression profile of pure cell lines
#                cpu_number = 9, # use multiple cores to save time. Set the cpu_number = length(cell_type_number) if there is enough cores.
#                print_progress_msg_to_file = 0)
#  

## ----fig3, echo=TRUE, eval=TRUE, fig.width=7,fig.height=5,fig.align = 'center'----
trueGEP <- true_GEP_gene[,result2$cell_type_assignment]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1),  pty="s")
for (i in 1:6) {
  max_v <- max(c(result2$estGEP[,i], trueGEP[,i]))
  plot(result2$estGEP[,i],
       trueGEP[,i],
       xlab = "CDSeq-estimated GEPs", 
       ylab = "True GEPs",
       pch = 21,
       bg = 'green',
       cex=2, 
       xlim = c(0, max_v), 
       ylim = c(0,max_v),
       main = paste0("cell type ",i))
  lines(c(0,max_v),c(0,max_v), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----fig4, echo=TRUE, eval=TRUE, fig.width=7,fig.height=5,fig.align = 'center'----
trueProp <- true_prop_cell[result2$cell_type_assignment,]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1),  pty="s")
for (i in 1:6) {
    plot(result2$estProp[i,],
         trueProp[i,],
         pch = 21,
         cex = 3,
         bg='red',
         xlab = "CDSeq-estimated proportions", 
         ylab = "True proportions", 
         xlim = c(0,1),
         ylim = c(0,1),
         main = paste0("cell type ",i))
    lines(c(0,1),c(0,1), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----fig5,echo=TRUE, eval=TRUE,  fig.width=7,fig.height=6,fig.align = 'center'----
lgpst <- rep(0,9)
for (i in 1:9) {
  lgpst[i] <- result2$est_all[[i]]$lgpst
}
plot(2:10,lgpst,xlab = "number of cell typess", ylab = "log posterior")
points(6,lgpst[5],pch=16,col="red")
lines(2:10, lgpst)
lines(c(6,6),c(0,lgpst[5]),lty=2)

## ----example 3, eval=FALSE, cache=FALSE, results='hide'-----------------------
#  result3<-CDSeq(bulk_data =  mixtureGEP,
#                cell_type_number = 6,
#                mcmc_iterations = 2000,
#                dilution_factor = 1,
#                block_number = 10,
#                gene_subset_size = 100,
#                gene_length = as.vector(gene_length),
#                reference_gep = refGEP,
#                cpu_number = 10,
#                print_progress_msg_to_file = 0)

## ----fig6, echo=TRUE, eval=TRUE, fig.width=7,fig.height=5,fig.align = 'center'----

trueGEP <- true_GEP_rpkm[,result3$cell_type_assignment]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1),  pty="s")
for (i in 1:6) {
  max_v <- max(c(result3$estGEP[,i], trueGEP[,i]))
  plot(result3$estGEP[,i],
       trueGEP[,i],
       xlab = "CDSeq-estimated GEPs", 
       ylab = "True GEPs",
       pch = 21,
       bg = 'green',
       cex=2, 
       xlim = c(0, max_v), 
       ylim = c(0,max_v),
       main = paste0("cell type ",i))
  lines(c(0,max_v),c(0,max_v), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----fig7, echo=TRUE, eval=TRUE, fig.width=7,fig.height=5,fig.align = 'center'----

trueProp <- true_prop_cell[result3$cell_type_assignment,]
.pardefault <- par()
par(mar=c(4,3,1,1), mfrow=c(2,3),mai = c(.1, .01, 0.1, 0.1),  pty="s")
for (i in 1:6) {
    plot(result3$estProp[i,],
         trueProp[i,],
         pch = 21,
         cex = 3,
         bg='red',
         xlab = "CDSeq-estimated proportions", 
         ylab = "True proportions", 
         xlim = c(0,1),
         ylim = c(0,1),
         main = paste0("cell type ",i))
    lines(c(0,1),c(0,1), type = "l")
}
par(mar = .pardefault$mar, mfrow = .pardefault$mfrow, mai=.pardefault$mai, pty=.pardefault$pty)


## ----example 4, echo=TRUE, eval=FALSE, cache=FALSE, results='hide'------------
#  cdseq.result <- CDSeq::CDSeq(bulk_data = pbmc_mix,
#                               cell_type_number = seq(3,12,3),
#                               beta = 0.5,
#                               alpha = 5,
#                               mcmc_iterations = 700,
#                               cpu_number = 4,
#                               dilution_factor = 10)

## ----fig8, echo=TRUE, eval=TRUE, results='hide',message=FALSE, warning=FALSE----
cdseq.result.celltypeassign <- vector(mode = "list", length = 4)

for (i in 1:4) {
  cdseq_gep <- cdseq.result$est_all[[i]]$estGEP
  cdseq_prop <- cdseq.result$est_all[[i]]$estProp
  cdseq.result.celltypeassign[[i]] <- cellTypeAssignSCRNA(cdseq_gep = cdseq_gep, # CDSeq-estimated cell-type-specific GEPs
                                                     cdseq_prop = cdseq_prop, # CDSeq-estimated cell type proportions
                                                     sc_gep = sc_gep,         # PBMC single cell data
                                                     sc_annotation = sc_annotation,# PBMC single data annotations
                                                     sc_pt_size = 3,
                                                     cdseq_pt_size = 6,
                                                     seurat_nfeatures = 100,
                                                     seurat_npcs = 50,
                                                     seurat_dims=1:5,
                                                     plot_umap = 1,
                                                     plot_tsne = 0)
  
  cdseq.result.celltypeassign[[i]]$cdseq_scRNA_umap
}

## ----show GEP vs scRNAseq, echo=TRUE, eval=TRUE, fig.width = 20, fig.height = 70, fig.align = 'center',out.width = "95%"----
ggarrange(cdseq.result.celltypeassign[[1]]$cdseq_scRNA_umap, 
          cdseq.result.celltypeassign[[2]]$cdseq_scRNA_umap,
          cdseq.result.celltypeassign[[3]]$cdseq_scRNA_umap,
          cdseq.result.celltypeassign[[4]]$cdseq_scRNA_umap,
          nrow = 4,ncol = 1, common.legend = TRUE) 

## ----prepare df for ggplot, echo=TRUE, eval=FALSE-----------------------------
#  for (i in 1:4) {
#    if(i==1){
#      rownames(cdseq.result.celltypeassign[[i]]$cdseq_prop_merged) <- rownames(true_prop)
#      cdseq_prop_df <- melt(cdseq.result.celltypeassign[[i]]$cdseq_prop_merged)
#      names(cdseq_prop_df)[1:3] <- c("celltype", "sampleID",paste0("estprop_",i))
#      true_prop_df <- melt(true_prop)
#      names(true_prop_df)[1:3] <- c("celltype", "sampleID","trueprop")
#  
#      merge_df <- merge(cdseq_prop_df,true_prop_df, by=c("celltype","sampleID"))
#    }else{
#      rownames(cdseq.result.celltypeassign[[i]]$cdseq_prop_merged) <- rownames(true_prop)
#      cdseq_prop_df <- melt(cdseq.result.celltypeassign[[i]]$cdseq_prop_merged)
#      names(cdseq_prop_df)[1:3] <- c("celltype", "sampleID",paste0("estprop_",i))
#  
#      merge_df <- merge(merge_df, cdseq_prop_df, by=c("celltype","sampleID"))
#    }
#  }
#  }

## ----compare GEP, echo=TRUE, eval=TRUE, results='hide', message=FALSE---------
pbmc_ggplot <- vector(mode = "list", length = 4)
for (i in 1:4) {
  B_plot<- merge_df %>% filter(celltype=="Bcell") %>% ggscatter(x = paste0("estprop_",i), y = "trueprop",
                                                                add = "reg.line",  
                                                                size = 3,
                                                                title = "B cell",
                                                                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                                                conf.int = TRUE 
  ) + stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01)
  
  T_plot<- merge_df %>% filter(celltype=="Tcell") %>% ggscatter(x = paste0("estprop_",i), y = "trueprop",
                                                                add = "reg.line", 
                                                                size = 3,
                                                                title = "T cell",
                                                                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                                                conf.int = TRUE # Add confidence interval
  ) + stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01)
  
  M_plot<- merge_df %>% filter(celltype=="Monocyte") %>% ggscatter(x = paste0("estprop_",i), y = "trueprop",
                                                                   add = "reg.line",  
                                                                   size = 3,
                                                                   title = "Monocyte",
                                                                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                                                                   conf.int = TRUE 
  ) + stat_cor(method = "pearson", p.accuracy = 0.001, r.accuracy = 0.01)
  
  
  pbmc_ggplot[[i]] <- ggarrange(B_plot, T_plot, M_plot, ncol = 3, nrow = 1)
}

## ----show plots, echo=TRUE, eval=TRUE, fig.width=8, fig.height=12, fig.align = 'center', out.width ="95%"----
ggarrange(pbmc_ggplot[[1]], 
          pbmc_ggplot[[2]], 
          pbmc_ggplot[[3]], 
          pbmc_ggplot[[4]],
          nrow = 4, ncol = 1)


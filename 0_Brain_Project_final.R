setwd("/Users/gabriel/Downloads/tsv/")
#Format got MetaXcan
dir.create("formatted")
files <- list.files("/Users/gabriel/Downloads/tsv/GWAS/")
for (i in 1:length(files)){
  print(i)
  token <- files[i]
  output <- read.table(paste("/Users/gabriel/Downloads/tsv/GWAS/",token,sep=""),
                       sep = "\t",
                       header = T)
  output_final <- output
  
  output_final$snpid <- output_final$rsids
  output_final$chr <- output_final$chrom
  output_final$a1 <- output_final$ref
  output_final$a2 <- output_final$alt
  output_final$freq <- output_final$maf
  output_final$beta <- output_final$beta
  output_final$se <- output_final$sebeta
  output_final$pval <- output_final$pval
  output_MTAG <- as.data.frame(cbind(output_final$snpid,
                                     output_final$chr,
                                     output_final$a1,
                                     output_final$a2,
                                     output_final$freq,
                                     output_final$beta,
                                     output_final$se,
                                     output_final$pval))
  colnames(output_MTAG) <- c("snpid","chr","a1","a2","freq","beta",
                             "se","pval")
  output_MTAG <- output_MTAG[output_MTAG$snpid != ".",]
  dir.create(paste("formatted/",token,sep=""))
  write.table(output_MTAG,file = paste("formatted/",token,"/output.txt",sep = ""),
              sep = " ",row.names = F,col.names = T, quote = F)
}
write.table(list.files("GWAS/"), file = "list.txt",sep= " ",
            quote = F,row.names = F,col.names = F)

#Read Correlation
setwd("/Volumes/Year4PhD/Brain_Project/")
files_tissue <- list.files("results/")
files_tissue <- substring(files_tissue[grepl("Brain",files_tissue)],1,
                          nchar(files_tissue[grepl("Brain",files_tissue)])-4)
files <- list.files("tsv/formatted/")

list_tissue_expression <- list()
for (j in 1:length(files_tissue)){
  tissue <- files_tissue[j]
  list_expression <- list()
  print(j)
  for (i in 1:length(files)){
    token <- read.csv(paste("tsv/formatted/",files[i],"/",tissue,".csv",sep=""))
    token$image <- files[i]
    list_expression[[i]] <- token
  }
  head(token)
  total <- do.call("rbind",list_expression)
  all_genes <- unique(total$gene_name)
  
  library(plyr)
  token_df <- as.data.frame(cbind(all_genes,0))
  colnames(token_df) <- c("gene_name","zscore")
  list_expression_new <- vector("list", length = 185)
  list_expression_new[[1]] <- token_df
  for (i in 1:length(list_expression)){
    token_new <- list_expression[[i]]
    token_new <- as.data.frame(cbind(token_new$gene_name,token_new$zscore))
    colnames(token_new) <- c("gene_name","zscore")
    list_expression_new[[i+1]] <- token_new
  }
  expression_matrix <- join_all(dfs = list_expression_new, by='gene_name', type='left',)
  gene_names_placeholder <- expression_matrix$gene_name
  expression_matrix <- expression_matrix[,-c(1:2)]
  expression_matrix[is.na(expression_matrix)] <- 0
  expression_matrix <- sapply( expression_matrix, as.numeric )
  colnames(expression_matrix) <- files
  rownames(expression_matrix) <- gene_names_placeholder
  list_tissue_expression[[j]] <- expression_matrix
  names(list_tissue_expression)[j] <- files_tissue[j]
}
save(list_tissue_expression, file = "list_tissue_expression.Rdata")

#PCA
load("list_tissue_expression.Rdata")
library("glmpca")
library(ggrepel)
pca_plot <- prcomp(list_tissue_expression[[5]],scale. = TRUE,center = TRUE)
pca_plot$terms <- rownames(pca_plot)
pca_dim <- as.data.frame(pca_plot$rotation)
pca_dim$terms <- rownames(pca_plot$rotation)

if(all.equal(pca_dim$terms,idp$phenotype)){
  pca_dim$region <- idp$region_category
  pca_dim$imaging_measure <- idp$`imaging measure`
  pca_dim$white_grey <- idp$white_grey_combined
}
library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
colors <- getPalette(12)
plot_1 <- ggplot(pca_dim, aes(x = PC1, y = PC2)) +
  geom_point(size = 3, aes(shape = region, col = imaging_measure)) +
  coord_fixed() + ggtitle("PCA of Cerebellar Imputation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = colors)
plot_2 <- ggplot(pca_dim, aes(x = PC1, y = PC3)) +
  geom_point(size = 3, aes(shape = region, col = imaging_measure)) +
  coord_fixed() + ggtitle("PCA of Cerebellar Imputation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = colors)
plot_3 <- ggplot(pca_dim, aes(x = PC2, y = PC3)) +
  geom_point(size = 3, aes(shape = region, col = imaging_measure)) +
  coord_fixed() + ggtitle("PCA of Cerebellar Imputation") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_color_manual(values = colors)

plot_total <- ggarrange(plot_1,plot_2,plot_3,
                        ncol = 3, nrow = 1,common.legend = TRUE, legend = "right")
pdf("pca_total.pdf",width = 15, height = 15)
plot_total
dev.off()
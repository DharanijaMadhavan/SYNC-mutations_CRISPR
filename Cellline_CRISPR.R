library(ggplot2)
library(data.table)
library(stringr)
library(dplyr)
library(tidyverse)
library(tibble)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(gridExtra)
##############################
#Reading in data
##############################
depmap_meta <- read.csv("DepMap-2019q1-celllines_v2.csv", header=T) 

data.mut.depmap <- read.csv("depmap_19Q1_mutation_calls.csv", header = T)

depmap_ceres <- read.csv("gene_effect_corrected.csv", header =T)
row.names(depmap_ceres) <- depmap_ceres[,1]
depmap_ceres <- depmap_ceres[,-1]


##############################
#Identifying cell lines with mutations
##############################
Gene <- "PCBP1"
mutations <- c("L100Q", "L102Q", "P99S") 


mut.depmap <- subset(data.mut.depmap, Hugo_Symbol == Gene) 
mut.depmap <- mut.depmap %>% 
            mutate("CCLE_Name"= slice(depmap_meta, match(mut.depmap$DepMap_ID, depmap_meta$DepMap_ID))$CCLE_Name) %>% 
            mutate("Primary.Disease"= slice(depmap_meta, match(mut.depmap$DepMap_ID, depmap_meta$DepMap_ID))$Primary.Disease)
                  
mut.depmap$Primary.Disease <-   gsub(" Cancer", "", mut.depmap$Primary.Disease) 



##For checking if rows matched correctly
#mut.meta <- slice(depmap_meta, match(mut.depmap$DepMap_ID, depmap_meta$DepMap_ID))
#as.vector(mut.depmap$DepMap_ID) == as.vector(mut.meta$DepMap_ID)

#Plotting all variants
varA <- ggplot(mut.depmap,aes(x=Protein_Change, y=Primary.Disease, colour=Variant_Classification))+
          geom_point()+theme(axis.text.x = element_text(angle = 90, hjust = 1),
                     axis.title = element_text(size=10),  axis.text=element_text(size=8), 
                     legend.position="bottom",legend.title = element_blank(),
                     legend.text = element_text(size=8))



#Plotting non-silent
varNS <- mut.depmap %>% 
  filter(Variant_Classification != "Silent") %>% 
  mutate(xy=ifelse(grepl(paste(mutations,collapse="|"),Protein_Change) == TRUE,  "SYNC Mutations","")) %>% 
  ggplot(aes(x=Protein_Change, y=Primary.Disease))+
  geom_point(aes(col=xy, shape= xy))+ scale_color_manual(values=c("black","red")) +
  theme(legend.position="bottom",legend.title = element_blank(),legend.text = element_text(size=8),
        axis.title = element_text(size=10),axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text=element_text(size=8))

  


#Cell Line information
Cell.lines <- mut.depmap %>% filter(grepl(paste(mutations,collapse="|"),Protein_Change)) %>%  
  select(c(CCLE_Name, DepMap_ID, Protein_Change)) %>% 
  distinct(DepMap_ID,.keep_all = T) #Removing duplicate rows wherein the only difference is codon change


print(Cell.lines)


mut.depmap2 <- mut.depmap %>% distinct(DepMap_ID,.keep_all = T) 


##########################################
#CERES distribution of Gene across all cell lines
###########################################
#tib of CERES score for gene of interest and corresponding DepmapID
gene_ceres <- depmap_ceres %>%  tibble::rownames_to_column("aa") %>% 
    select(c(aa,contains(Gene))) %>% rename( Gene = contains(Gene))  
 
 
 
 #Ceres score of gene of interest in only the non-silent mutant celll lines 
gene_mut_ceres <-  gene_ceres %>% 
                  slice(match(filter(mut.depmap2, Variant_Classification != "Silent")$DepMap_ID,
                              gene_ceres$aa)) %>% 
                             column_to_rownames("aa") 

colour <- ifelse(gene_mut_ceres$Gene < -0.5, "red", "black")  
labels <- depmap_meta[match(row.names(gene_mut_ceres), depmap_meta$DepMap_ID),"CCLE_Name"]
labels <- paste(labels,mut.depmap2[match(row.names(gene_mut_ceres), mut.depmap2$DepMap_ID),"Protein_Change"],
                      sep=":")
ceres.dist <-  ggplot(gene_ceres, aes(x = Gene)) + geom_histogram(colour= "black", fill="white") + 
        labs(x= " CERES Distribution")+theme(axis.title = element_text(size=10))+
  theme(axis.text = element_text(size=8))+
  geom_vline(aes(xintercept = Gene), gene_mut_ceres, col=colour)+
   annotate("text", x=gene_mut_ceres$Gene, y=70, label=labels, angle=90, vjust=-0.5, size=3) #only if legible
   
   title <- ggdraw() +  draw_label(Gene, fontface='bold')
plot_grid(title, plot_grid(varA, varNS, ncol=2),NULL, ceres.dist, nrow = 4, rel_heights = c(0.2, 1,0.02, 0.8))
 
########################################
#Getting relevant CRISPR data of cell line with SYNC mutations
########################################
  
mut.ceres <- slice(depmap_ceres, match(Cell.lines$DepMap_ID, row.names(depmap_ceres))) 


rownames(mut.ceres) <- row.names(depmap_ceres)[na.omit(match( Cell.lines$DepMap_ID,
                                                         row.names(depmap_ceres)))]
Cell.lines_crispr2 <- Cell.lines[na.omit(match(row.names(depmap_ceres),Cell.lines$DepMap_ID)),
                                c("CCLE_Name","Protein_Change")] 
                        
  
Cell.lines_crispr <- paste(Cell.lines_crispr2$CCLE_Name, Cell.lines_crispr2$Protein_Change, sep=":")



############
#Statistical test
############
depmap_ceres2 <- depmap_ceres
#Matching with CERES_null...rda
colnames(depmap_ceres2) <- gsub("\\..*","", colnames(depmap_ceres2))
mm_set2 <- mm_set[match(colnames(depmap_ceres2), mm_set$symbol) , ] 

which((colnames(depmap_ceres2) == mm_set2$symbol) == FALSE)


##################################################
#Comparing distributions of genes and plotting
##################################################
mut.ceres.list <- as.list(as.data.frame(t(mut.ceres)))


#Calculate p value for each gene/CERES in a cell line against distribution of that gene's
#CERES across all other cell lines
pval <- lapply(mut.ceres.list, function(x)
  2*pnorm(x, mean= mm_set2$mu_null, sd= mm_set2$sigma_null ) )

###Adjusted pValue
adj.pval <- lapply(pval, function(x) data.frame("FDR" = p.adjust(x, "fdr")) %>% 
          mutate("sig" = ifelse(FDR < 0.1, "sig", "ns")))  
        
#Make list of data frames for plotting
mut_ceres <- lapply(mut.ceres.list, function(x) 
  data.frame("dCERES" = x-mm_set$mu_null, "Gene" = gsub("\\..*","", colnames(mut.ceres))) )

mut_ceres_pval <- Map(cbind, mut_ceres, adj.pval)

##Plotting 

plot_dceresDist <- list()
for (i in 1:length(mut_ceres_pval)){
  local({
    i <- i
p1 <- ggplot(mut_ceres_pval[[i]],aes(x=dCERES))+ geom_histogram()+
  ggtitle(Cell.lines_crispr[i])+
       xlab("\u0394CERES distribution")+
      theme(axis.title.x  = element_blank(), 
            axis.title.y = element_text(size=10), axis.text = element_text(size=10))
    print(i)
    print(p1)
    plot_dceresDist[[i]] <<- p1  # add each plot into plot list
  })}


plot_volcano <- list()
for (i in 1:length(mut_ceres_pval)){
  local({
    i <- i
    p1 <- ggplot(mut_ceres_pval[[i]], aes(x= dCERES, y= -log10(FDR), colour = sig))+
  geom_point(alpha=0.4,size=1.75)+xlab("CERES") + ylab("-log10 FDR")+
  geom_text_repel(data= subset(mut_ceres_pval[[i]], mut_ceres_pval[[i]]$FDR < 0.1),
                  aes(dCERES,y=-log10(FDR),label=Gene), size=3) +
      theme(legend.position="none",axis.title.x  = element_blank(), 
            axis.title.y = element_text(size=10), axis.text = element_text(size=10))

    print(i)
    print(p1)
    plot_volcano[[i]] <<- p1 
  })}
  
 plot_combine <- rbind(plot_dceresDist, plot_volcano)
 

nrows <- ifelse(length(mut_ceres_pval) == 1, 1,  2 )


marrangeGrob(grobs = plot_combine, ncol=2, nrow=nrows, bottom ="\u0394CERES (cell line/\u03bc)")


##############################
#Results : Significant genes
##############################

sl <- lapply(mut_ceres_pval, function(x) subset(x, x$FDR < 0.1) %>% select(c(Gene,dCERES, FDR))) 
names(sl) <- Cell.lines_crispr
             
 sl <- lapply(sl, function(x) {
           cbind(x, "dCERES" = round(x$dCERES,2), "FDR"= scientific(x$FDR, digits=3))}) %>% 
      lapply(function(x) x[,-c(2,3)])


list(as.factor(Reduce(intersect, sapply(sl, function(x) x$Gene))))


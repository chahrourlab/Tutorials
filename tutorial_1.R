################################################
#####Tutorial 1- 8/24/23
#####Volcano plot
################################################

#Download libraries
install.packages("ggplot2")
install.packages("ggrepel")

#Load libraries
library(ggplot2)
library(ggrepel)

#Set working directory
setwd("/data")

#Read the data we will use for plotting
df<- read.table(file="pseudobulk_allcells_kovswt.txt", header = TRUE, sep = "\t")

#Lets see what the data looks like
View(df)

#Basic scatter plot
ggplot(data=df,aes(x=log2FoldChange,y=padj))+geom_point()

#Doesn't look like a volcano plot. What is different? 
#Change y axis scale
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()

#Lets get rid of the background grid lines 
#Add vertical lines that indicate log2FoldChange threshold and horizontal lines for padj threshold
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj))) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) 

#Which genes are differentially expressed? The ones in the upper left and upper right corner
#Add a column to the table to specify if a gene is differentially expressed
df$Differential_Expression<-"No"
df$Differential_Expression[df$log2FoldChange>0.3 & df$padj<0.05 & df$pvalue<0.05]<-"Upregulated"
df$Differential_Expression[df$log2FoldChange< -0.3 & df$padj<0.05 & df$pvalue<0.05]<-"Downregulated"

#Lets see the updated table
View(df)

#Re-plot with "colour" parameter- this will colour the DEGs
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=Differential_Expression)) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) 

#Lets change these default colours and re-plot
#Get rid of the "No" in the legend
mycolors <- c("black", "darkred","darkblue")
names(mycolors) <- c("No", "Upregulated","Downregulated")
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=Differential_Expression)) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) +
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated"))

#Subset upregulated genes and order their log2FoldChange in descending order
up<-df[df$Differential_Expression=="Upregulated",]
up<-up[order(up$log2FoldChange,decreasing = TRUE),]

#Subset downregulated genes and order their log2FoldChange in ascending order
down<-df[df$Differential_Expression=="Downregulated",]
down<-down[order(down$log2FoldChange),]

#Lets make a list of the top 10 upregulated genes and top 10 downregulated genes
top_genes<-append((up[1:10,1]),(down[1:10,1]))
#instead of this, you can also label specific genes of your choice by making a list of gene names like you did for colors. 

#Write the names of these top DEGs next to their point
df$delabel<-NA
#check if the value in the gene column exists in the top_genes list
#If it does, change the NA label to the gene name
df$delabel[df$gene %in% top_genes]<-df$gene[df$gene %in% top_genes]

#Lets look at the updated table
View(df)

#Plot with labels
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=Differential_Expression, label=delabel)) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) +
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated")) + 
  geom_text_repel(arrow = arrow(length = unit(0.02, "npc")),box.padding = 1,max.overlaps=Inf)

#Change titles and clean up
ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=Differential_Expression, label=delabel)) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) +
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated")) + 
  geom_text_repel(arrow = arrow(length = unit(0.02, "npc")),box.padding = 1,max.overlaps=Inf) +
  xlim(-2.5,2.5)+
  xlab("Log2 fold change")+
  ylab("-log10(adjusted p-value)")+
  ggtitle("Differentially expressed genes KO v/s WT")+ 
  theme(legend.title=element_blank())

#save the plot
p<- ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=Differential_Expression, label=delabel)) + geom_point()+
  theme_classic()+
  geom_vline(xintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  geom_hline(yintercept=-log10(0.05), col="darkgreen",linetype=2) +
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated")) + 
  geom_text_repel(arrow = arrow(length = unit(0.02, "npc")),box.padding = 1,max.overlaps=Inf) +
  xlim(-2.5,2.5)+
  xlab("Log2 fold change")+
  ylab("-log10(adjusted p-value)")+
  ggtitle("Differentially expressed genes KO v/s WT")+ 
  theme(legend.title=element_blank())

png(file="volcano_plot.png",width=800,height=800,res=100)
print(p)
dev.off()

################################################
#####Challenge 1- 7/27/23
#####MA plot
################################################
top_genes<-c("Sox18","Foxq1","Isl1","Drd2")
df$delabel<-NA
df$delabel[df$gene %in% top_genes]<-df$gene[df$gene %in% top_genes]

q<- ggplot(data=df, aes(x=log2(baseMean), y=log2FoldChange, col=Differential_Expression, label=delabel)) + geom_point()+
  theme_classic()+
  geom_hline(yintercept=c(-0.3, 0.3), col="darkgreen",linetype=2) +
  scale_colour_manual(values = mycolors,breaks = c("Upregulated","Downregulated")) + 
  geom_text_repel(arrow = arrow(length = unit(0.02, "npc")),box.padding = 1,max.overlaps=Inf)+
  xlab("Log2 base mean")+
  ylab("Log2 fold change")+
  ggtitle("Differentially expressed genes KO v/s WT")+ 
  theme(legend.title=element_blank()) 

png(file="ma_plot.png",width=800,height=800,res=100)
print(q)
dev.off()

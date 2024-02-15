

#library(phyloseq)
library(MicrobiotaProcess) # an R package for analysis, visualization and biomarker discovery of Microbiome.

library(phyloseq) # Handling and analysis of high-throughput microbiome census data.
library(ggplot2) # Create Elegant Data Visualisations Using the Grammar of Graphics.
library(tidyverse) # Easily Install and Load the 'Tidyverse'.
library(vegan) # Community Ecology Package.
library(coin) # Conditional Inference Procedures in a Permutation Test Framework.
library(reshape2) # Flexibly Reshape Data: A Reboot of the Reshape Package.
library(ggnewscale) # Multiple Fill and Colour Scales in 'ggplot2'.
library(phyloseq)
library(dplyr)
library(UpSetR)
library(microbiome)
library(microViz)
library(microbiome)
library(dplyr)







setwd('D:/base/cystcase/cyststoolfastqs/analysis/case144')
otu <- read.csv('tax_case144.csv',row.names = 1,check.names = FALSE)
meta <- read.csv('meta_case144.csv',row.names = 1)
taxa <- read.csv('linm_case144.csv',row.names = 1)


dim(otu)
dim(taxa)
OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(as.matrix(taxa))
samples = sample_data(meta)
head(OTU)

ps=phyloseq::phyloseq(OTU, TAX, samples)

classtaxa <- get_taxadf(obj=ps, taxlevel=6)




dev.new()

alphaobj <- get_alphaindex(classtaxa)
head(as.data.frame(alphaobj))
write.csv(as.data.frame(alphaobj),'alphaobj_enterotype.csv')


p_alpha <- ggbox(alphaobj, geom="violin", factorNames="Enterotype") +
  scale_fill_manual(values=c("#00AED7", "#FD9347"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p_alpha

distme <- get_dist(classtaxa, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(classtaxa), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$Enterotype <- factor(sampleda$Enterotype)
set.seed(1024)
adores=0
library(vegan)

adores <- adonis2(distme ~ Enterotype, data=sampleda, permutation=9999)
data.frame(adores)


pcares <- get_pca(obj=classtaxa, method="hellinger")
# Visulizing the result
pcaplot1 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
                       factorNames=c("Enterotype"), ellipse=TRUE) +
  scale_color_manual(values=c("#00AED7", "#FD9347")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347"))
pcaplot1

head(pcoares@sampleda)
pcoares <- get_pcoa(obj=classtaxa, distmethod="bray", method="hellinger")
# Visualizing the result
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("Enterotype"), ellipse=TRUE) +
  scale_color_manual(values=c("#00AED7", "#FD9347")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347"))

pcoaplot1


library(ggplot2)
library(MicrobiotaProcess)
library(VennDiagram)
library(UpSetR)

vennlist <- get_vennlist(obj=classtaxa, factorNames="Enterotype")
upsetda <- get_upset(obj=classtaxa, factorNames="Enterotype")

vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "#FD9347"),
                      cat.col=c("#00AED7", "#FD9347"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c("#00AED7", "chartreuse", "blue4"),
                      cat.col=c("#00AED7", "chartreuse", "blue4"),
                      alpha = 0.85, 
                      fontfamily = "serif",
                      fontface = "bold",
                      cex = 1.2,
                      cat.cex = 1.3,
                      cat.default.pos = "outer",
                      cat.dist=0.1,
                      margin = 0.1, 
                      lwd = 3,
                      lty ='dotted',
                      imagetype = "svg")
dev.new()
grid::grid.draw(vennp)

upset(upsetda, sets=unique(as.vector(sample_data(classtaxa)$Enterotype)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")

write.csv(upsetda, file = "upsett_vennlist_Enterotype.csv")























library(cluster)
library(clusterSim)
list.files()
#Download the example data and set the working directory
#setwd('<path_to_working_directory>')
data=read.table("sampleentero.txt", header=T, row.names=1, dec=".", sep="\t")
data=data[-1,]
dim(data)
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

data.dist=dist.JSD(data)

pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  require(cluster)
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
  return(cluster)
}

data.cluster=pam.clustering(data.dist, k=2)

require(clusterSim)
nclusters = index.G1(t(data), data.cluster, d = data.dist, centrotypes = "medoids")

nclusters=NULL

for (k in 1:20) { 
  if (k==1) {
    nclusters[k]=NA 
  } else {
    data.cluster_temp=pam.clustering(data.dist, k)
    nclusters[k]=index.G1(t(data),data.cluster_temp,  d = data.dist,
                          centrotypes = "medoids")
  }
}
dev.new()
plot(nclusters, type="h", xlab="k clusters", ylab="CH index",main="Optimal number of clusters")

obs.silhouette=mean(silhouette(data.cluster, data.dist)[,3])
cat(obs.silhouette) #0.1899451

cat(obs.silhouette) #0.1899451

#data=noise.removal(data, percent=0.01)

## plot 1
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
dev.new()
dim(obs.bet$ls)
dim(obs.pca$ls)
dim(as.factor(data.cluster))
s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F,sub="Between-class analysis")

#plot 2
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
dev.new()
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")














#data=noise.removal(data, percent=0.01)
library(ade4)
library("factoextra")
## plot 1
obs.pca=dudi.pca(data.frame(t(data)), scannf=F, nf=10)
res.pca=obs.pca
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

fviz_pca_var(res.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)


s.label(res.pca$li, 
        xax = 1,     # Dimension 1
        yax = 2)

fviz_pca_ind(obs.bet, repel = TRUE)
dev.new()
fviz_pca_ind(obs.bet,
             col.ind = as.factor(data.cluster), # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
)
groups <- as.factor(decathlon2$Competition[1:23])
s.class(res.pca$li,
        fac = as.factor(data.cluster),  # color by groups
        col = c("#00AFBB", "#E7B800", "#FC4E07")
)



obs.bet=bca(obs.pca, fac=as.factor(data.cluster), scannf=F, nf=k-1) 
dev.new()

a=s.class(obs.bet$ls, fac=as.factor(data.cluster), clabel = 1, grid=F,sub="Between-class analysis")
a

s.class(obs.bet$li, fac = as.factor(data.cluster), col = rainbow(12), facets = as.factor(data.cluster))

write.csv(obs.bet$ls,'pca.csv')
write.csv(data.cluster,'fact.csv')

plot(table(rpois(100,5)), type = "h", col = "red", lwd=10,
     main="rpois(100,lambda=5)")

x <- obs.bet$ls$CS1
y <- obs.bet$ls$CS2
plot(x, y)
+text(x, y, labels=as.factor(data.cluster))

d=obs.bet$ls
dim(obs.bet$ls)
typeof(data.cluster)
length(data.cluster)
d$row=s.class(obs.bet$ls, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(4,2,3))
#plot 2
obs.pcoa=dudi.pco(data.dist, scannf=F, nf=3)
dev.new()
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F,sub="Principal coordiante analysis")
s.class(obs.pcoa$li, fac=as.factor(data.cluster), grid=F, cell=0, cstar=0, col=c(3,2,4))

d$row=rownames(d)
d$ent=as.factor(data.cluster)

write.csv(d,'ent.pca.csv')

library(ggplot2)
library(ggrepel)
ggplot(d, aes(CS1,CS2)) +
  geom_point() +
  geom_text_repel(aes(label = row))
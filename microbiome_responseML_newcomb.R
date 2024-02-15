remotes::install_github("yiluheihei/microbiomeMarker")
setwd('D:/base/cystcase/responseML/newcomb')
otu <- read.csv('taxgenus.csv',row.names = 1,check.names = FALSE)
meta <- read.csv('metagenus.csv',row.names = 1)
taxa <- read.csv('lingenus.csv',row.names = 1)

OTU = otu_table(as.matrix(otu), taxa_are_rows = TRUE)
TAX = phyloseq::tax_table(as.matrix(taxa))
samples = sample_data(meta)
ps=phyloseq(OTU, TAX, samples)
classtaxa <- get_taxadf(obj=ps, taxlevel=6)


relab_genera = transform_sample_counts(classtaxa, function(x) x / sum(x) * 100) 
relab_genera = transform_sample_counts(classtaxa, function(x) x / sum(x)) 
head(otu_table(relab_genera))
rel_abundance_df <- as.data.frame(t(otu_table(relab_genera)))
rel_abundance_df$sample_names <- sample_names(relab_genera)
rel_abundance_df$NCR <- relab_genera@sam_data$NCR
rel_abundance_by_Type <- aggregate(rel_abundance_df[, -ncol(rel_abundance_df)],
                                   by = list(NCR = rel_abundance_df$NCR),
                                   FUN = mean)

# Print the relative abundance of genera by 'Type'
print(rel_abundance_by_Type)
write.csv(rel_abundance_df, file = "relative_data.csv", row.names = FALSE)

write.csv(rel_abundance_by_Type, file = "relative_abundance_combined_NCR.csv", row.names = FALSE)



dev.new()

pclass <- ggbartax(obj=classtaxa, facetNames="NCR") +
  xlab(NULL) +
  ylab("relative abundance (%)") +
  scale_fill_manual(values=c(colorRampPalette(RColorBrewer::brewer.pal(12,"Set3"))(31))) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5))
pclass


alphaobj <- get_alphaindex(classtaxa)
head(as.data.frame(alphaobj))
write.csv(as.data.frame(alphaobj),'alphaobj_CaseFC_combined_NCR.csv')


p_alpha <- ggbox(alphaobj, geom="violin", factorNames="NCR") +
  scale_fill_manual(values=c("#00AED7", "#FD9347","green"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))
p_alpha



dev.new()
pcares <- MicrobiotaProcess::get_pca(obj=classtaxa, method="hellinger")
# Visulizing the result
pcaplot1 <- ggordpoint(obj=pcares, biplot=TRUE, speciesannot=TRUE,
                       factorNames=c("NCR"), ellipse=TRUE) +
  scale_color_manual(values=c("#00AED7", "#FD9347","green")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347","green"))
pcaplot1

head(pcoares@sampleda)
pcoares <- get_pcoa(obj=classtaxa, distmethod="bray", method="hellinger")
# Visualizing the result
pcoaplot1 <- ggordpoint(obj=pcoares, biplot=TRUE, speciesannot=TRUE,
                        factorNames=c("NCR"), ellipse=TRUE,showsample=FALSE) +
  scale_color_manual(values=c("#00AED7", "#FD9347","green")) +
  scale_fill_manual(values=c("#00AED7", "#FD9347","green"))

pcoaplot1


distme <- get_dist(classtaxa, distmethod ="bray", method="hellinger")
sampleda <- data.frame(sample_data(classtaxa), check.names=FALSE)
sampleda <- sampleda[match(colnames(as.matrix(distme)),rownames(sampleda)),,drop=FALSE]
sampleda$NCR <- factor(sampleda$NCR)
set.seed(1024)
adores=0
library(vegan)

adores <- adonis2(distme ~ NCR, data=sampleda, permutation=9999)
data.frame(adores)


library(VennDiagram)
library(UpSetR)

vennlist <- get_vennlist(obj=classtaxa, factorNames="NCR")
upsetda <- get_upset(obj=classtaxa, factorNames="NCR")
vennp <- venn.diagram(vennlist,
                      height=5,
                      width=5, 
                      filename=NULL, 
                      fill=c(c("#ADD8E6", "#90EE90")),
                      cat.col=c("#000000", "#000000"),
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


upset(upsetda, sets=unique(as.vector(sample_data(classtaxa)$NCR)), 
      sets.bar.color = "#56B4E9",
      order.by = "freq", 
      empty.intersections = "on")

write.csv(upsetda, file = "upsett_vennlist_combined_NCR.csv.csv")


psg <- tax_glom(ps, taxrank="Genus")
relab_genera = transform_sample_counts(psg, function(x) x / sum(x)) 
head(otu_table(relab_genera))

rel_abundance_df <- as.data.frame(t(otu_table(relab_genera)))
#rel_abundance_df$sample_names <- sample_names(relab_genera)
rel_abundance_df$NCR <- relab_genera@sam_data$NCR

write.csv(rel_abundance_df, file = "relative_abundance_combined_genus_lefse_NCR.csv", row.names = FALSE)


classtaxa
library(microbiomeMarker)
psd=ps_filter(ps, Cohort %in% c("Case"))







mm_lr <- run_sl(
  psd,
  group = "NCR",
  nfolds = 2,
  nrepeats = 1,
  taxa_rank = "Genus",
  top_n = 15,
  norm = "TSS",
  method = "LR",
)


marker_table(mm_lr)



set.seed(2021)
#dev.new()
plot_sl_roc(mm_lr, group = "NCR")

pht <- run_posthoc_test(classtaxa, group = "NCR")
pht
plot_postHocTest(pht,feature="")
extract_posthoc_res(pht, "Cutibacterium")[[1]]

extract_posthoc_res(pht)
dim(pht@abundance)

pht@result

show(pht)

dim(as.data.frame(pht))


trace(run_sl, edit = T)


library(pROC)
abundance_matrix <- as.matrix(otu_table(classtaxa))
classes <- sample_data(ps)$NCR
taxon_abundance <- abundance_matrix["Genus", ]
classification <- ifelse(taxon_abundance > 0, "Diseased", "Healthy")
library("data.table")   # Also requires R.utils to read gz and bz2 files
library("phyloseq")
library("ALDEx2")
library("tidyverse")
setwd('D:/base/lauramicrobiomeinsights/picrust_pipeline')
picrust2 = "picrust_microinsight"

list.files(picrust2, recursive = TRUE)


p2_EC = paste0(picrust2, "/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_KO = paste0(picrust2, "/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz")
p2_PW = paste0(picrust2, "/pathways_out/path_abun_unstrat.tsv.gz")

mapfile = "description/"
list.files(mapfile, recursive = TRUE)
mapfile_EC = paste0(mapfile, "/ec_level4_info.tsv.gz")
mapfile_KO = paste0(mapfile, "/ko_info.tsv.gz")
mapfile_PW = paste0(mapfile, "/metacyc_pathways_info.txt.gz")


mapEC = as.data.frame(fread(mapfile_EC, header = FALSE))
colnames(mapEC) = c("function","description")
mapKO = as.data.frame(fread(mapfile_KO, header = FALSE, sep = "\t"))
colnames(mapKO) = c("function","description")
mapPW = as.data.frame(fread(mapfile_PW, header = FALSE))
colnames(mapPW) = c("pathway","description")




p2EC = as.data.frame(fread(p2_EC))
rownames(p2EC) = p2EC$"function"
p2EC = as.matrix(p2EC[,-1])
p2EC = round(p2EC)

p2KO = as.data.frame(fread(p2_KO))
rownames(p2KO) = p2KO$"function"
p2KO = as.matrix(p2KO[,-1])
p2KO = round(p2KO)

p2PW = as.data.frame(fread(p2_PW))
rownames(p2PW) = p2PW$"pathway"
p2PW = as.matrix(p2PW[,-1])
p2PW = round(p2PW)



p2EC1 = p2EC
p2KO1 = p2KO
p2PW1 = p2PW

ps0=read_rds('phyloseq_ps0.rds')


set.seed(12345)
system.time({
  aldex2_EC1 = aldex(p2EC1, sample_data(ps0)$Group, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})



# KO
set.seed(12345)
system.time({
  aldex2_KO1 = aldex(p2KO1, sample_data(ps0)$Group, mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})


# Pathway
set.seed(12345)
system.time({
  aldex2_PW1 = aldex(p2PW1, sample_data(ps0)$Group[-18], mc.samples = 500, test = "t", 
                     effect = TRUE, denom = "iqlr", verbose = TRUE)
})



#Check estimated effect size

png("images/ALDEx2_picrust2_effect_1.png", width = 6, height = 6, units = "in", res = 300)
par(mfrow = c(2,2))
hist(aldex2_EC1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "EC")
hist(aldex2_KO1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "KO")
hist(aldex2_PW1$effect, breaks = 20, xlab = "effect size", col = "yellow", main = "Pathway")
invisible(dev.off())










#Create MW and MA plots

png("images/ALDEx2_picrust2_MW_MA_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
aldex.plot(aldex2_EC1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(EC) MW Plot")

aldex.plot(aldex2_EC1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Log-ratio abundance", ylab = "Difference")
title(main = "(EC) MA Plot")

aldex.plot(aldex2_KO1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(KO) MW Plot")

aldex.plot(aldex2_KO1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(KO) MA Plot")

aldex.plot(aldex2_PW1, type = "MW", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Dispersion", ylab = "Difference")
title(main = "(PW) MW Plot")

aldex.plot(aldex2_PW1, type = "MA", test = "wilcox", all.cex = 0.4, rare.cex = 0.4, 
           called.cex = 0.6, cutoff = 0.05, xlab = "Relative abundance", ylab = "Difference")
title(main = "(PW) MA Plot")
invisible(dev.off())




#Relationship between effect, difference, and P values

png("images/ALDEx2_picrust2_P_adjP_1.png", width = 6, height = 8, units = "in", res = 300)
par(mfrow = c(3,2))
plot(aldex2_EC1$effect, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(EC) Effect size plot")
points(aldex2_EC1$effect, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_EC1$diff.btw, aldex2_EC1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(EC) Volcano plot")
points(aldex2_EC1$diff.btw, aldex2_EC1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_KO1$effect, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(KO) Effect size plot")
points(aldex2_KO1$effect, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_KO1$diff.btw, aldex2_KO1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(KO) Volcano plot")
points(aldex2_KO1$diff.btw, aldex2_KO1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")

plot(aldex2_PW1$effect, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Effect size", ylab = "P value", main = "(PW) Effect size plot")
points(aldex2_PW1$effect, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
legend("bottom", legend = c("P value", "BH-adjusted"), pch = 19, col = c("blue", "red"))

plot(aldex2_PW1$diff.btw, aldex2_PW1$we.ep, log = "y", cex = 0.4, col = "blue", pch = 19, 
     xlab = "Difference", ylab = "P value", main = "(PW) Volcano plot")
points(aldex2_PW1$diff.btw, aldex2_PW1$we.eBH, cex = 0.5, col = "red", pch = 19)
abline(h = 0.05, lty = 2, col = "grey")
invisible(dev.off())







# Merge with map file data

df_EC1 = aldex2_EC1 %>% tibble::rownames_to_column(var = "EC") %>% 
  inner_join(mapEC, by = c("EC" = "function")) %>% arrange(EC)

df_KO1 = aldex2_KO1 %>% tibble::rownames_to_column(var = "KO") %>% 
  inner_join(mapKO, by = c("KO" = "function")) %>% arrange(KO)

df_PW1 = aldex2_PW1 %>% tibble::rownames_to_column(var = "Pathway") %>% 
  inner_join(mapPW, by = c("Pathway" = "pathway")) %>% arrange(Pathway)



#Output to file
# Subset-1
write.table(df_EC1, file = "images/ALDEx2_picrust2_EC_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.csv(df_EC1, file = "images/ALDEx2_picrust2_EC_results_1.csv")
write.table(df_KO1, file = "images/ALDEx2_picrust2_KO_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.csv(df_KO1, file = "images/ALDEx2_picrust2_KO_results_1.csv")
write.table(df_PW1, file = "images/ALDEx2_picrust2_Pathway_results_1.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)
write.csv(df_PW1, file = "images/ALDEx2_picrust2_Pathway_results_1.csv")

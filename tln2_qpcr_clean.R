
cbbFill2 <- scale_fill_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#ABABAB", "#000080" ,"#000000", "#E69F00", 
                                       "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "thistle4", "yellowgreen", "violetred1"))
cbbColour2 <- scale_colour_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#ABABAB",
                                           "#000080", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00"))

## use this script for tln2 visualization and qpcr analysis

## read in qpcr data
dat = read.csv(file = "20191219_085241_CT039093_TIPTON_TLN2 -  Quantification Cq Results_0.csv", header = T)
dat = dat[,c(2:6, 8, 11:12)]

## 2nd round of qpcr - use this as standard. starting amount recorded
dat2 = read.csv(file = "20200116_084916_CT039093_TIPTON_TLN2 -  Quantification Cq Results_0.csv", header = T)
dat2 = dat2[,c(2:6, 8, 11:12)]
# fix dat2 target
dat2$Target = paste("TLN2_", dat2$Target, sep = "")
# adjust starting quantities
Sample = c("S1", "S2", "S3", "S4", "S5", "S6")
sdquant = c(50, 10, 2, 0.4, 0.08, 0.016)
sd10 = log10(sdquant)
sd5 = log(sdquant, base = 5)
sd = data.frame(Sample, sdquant, sd10, sd5)

## split out standards and per primer
std = dat2[which(dat2$Content == "Std"),]
std = merge(std, sd, by = "Sample")

std1 = std[which(std$Target == "TLN2_209"),]
# std1$group = "Single"
# std1$group[grep("P4\\.*", std1$Sample)] = "Multiplex"

ggplot(std1, aes(sd10, Cq)) + geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) + 
  scale_y_continuous(breaks = seq(26, 38, 2))
m = (lm(Cq ~ sd10, std1))
summary(m)
me209 = 10^(-1/m$coefficients[2]) # 99.0%
# se209 = 10^(-1/m$coefficients[4])
# sapply(split(data.frame(std1$Cq, std1$Log.Starting.Quantity), std1$group), function(x) summary(lm(x))$r.sq)

## repeat for 3var
std2 = std[which(std$Target == "TLN2_3var"),]
# std2$group = "Single"
# std2$group[grep("P4\\.*", std2$Sample)] = "Multiplex"

ggplot(std2, aes(sd10, Cq)) + geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) + 
  scale_y_continuous(breaks = seq(26, 38, 2))
m = (lm(Cq ~ sd10, std2))
summary(m)
me3var = 10^(-1/m$coefficients[2])
# se3var = 10^(-1/m$coefficients[4])
me3var # r2 = 99.7%
# se3var
# sapply(split(data.frame(std2$Cq, std2$Log.Starting.Quantity), std2$group), function(x) summary(lm(x))$r.sq)

## repeat for 4var
std3 = std[which(std$Target == "TLN2_4var"),]
# std3$group = "Single"
# std3$group[grep("P4\\.*", std3$Sample)] = "Multiplex"

ggplot(std3, aes(sd10, Cq)) + geom_point(size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", se = F) + 
  scale_y_continuous(breaks = seq(26, 38, 2))
m = (lm(Cq ~ sd10, std3))
summary(m)
me4var = 10^(-1/m$coefficients[2])
# se4var = 10^(-1/m$coefficients[4])
me4var # R2 = 99.5%
# se4var
# sapply(split(data.frame(std3$Cq, std3$Log.Starting.Quantity), std3$group), function(x) summary(lm(x))$r.sq)

## make standard curve plot with all assays
std$assay = std$Target
std$assay = gsub("TLN2_3var", "3 Isoforms", std$assay)
std$assay = gsub("TLN2_4var", "4 Isoforms", std$assay)
std$assay = gsub("TLN2_209", "Canonical\nIsoform", std$assay)

p4 = ggplot(std, aes(sd10, Cq, color = assay)) + 
    geom_smooth(aes(color = assay) ,method = "lm", se = F, alpha = 0.7) +
  geom_point(aes(color = assay, shape = assay), size = 2, alpha = 0.7) +
  scale_y_continuous(breaks = seq(26, 38, 3)) +
  labs(x = "Log10 TLN2 Amplicon", y = "Ct") + cbbColour2 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom", axis.title = element_text(size = 8), axis.text = element_text(size = 8),
        legend.text = element_text(size = 8), legend.title = element_blank()) +
  annotate("text", label = paste("Canonical: E = ", (round(me209/2, 2)*100), "%, R2 = 0.990", sep = ""), 
                   x = 0.59, y = 39, size = 2.5, color = "black") +
  annotate("text", label = paste("3 Isoforms: E = ", (round(me3var/2, 2)*100), "%, R2 = 0.997", sep = ""), 
           x = 0.55, y = 37.5, size = 2.5, color = "black") +
  annotate("text", label = paste("4 Isoforms: E = ", (round(me4var/2, 2)*100), "%, R2 = 0.995", sep = ""), 
           x = 0.55, y = 36, size = 2.5, color = "black")

#####
## calculate logratio (3var/4var - 209)
#####

## cut out standards and blank wells for analysis (9th and 10th columns of 1st round, 4th and 5th for 2nd)
wells = unique(dat$Well)
wells1 = wells[grep("9", wells)]
wells2= wells[grep("10", wells)]
wells = c(as.character(wells1), as.character(wells2))

dat = dat[which(dat$Well %in% wells),]

wells = unique(dat2$Well)
wells1 = wells[grep("4", wells)]
wells2= wells[grep("5", wells)]
wells = c(as.character(wells1), as.character(wells2))

dat2 = dat2[which(dat2$Well %in% wells),]

## match genotypes
dat$Genotype = dat$Sample
dat$Genotype = as.character(dat$Genotype)
# code unique to first qpcr run
dat$Genotype[grep("C.*", dat$Genotype)] = "CC"
dat$Genotype[grep("A.*", dat$Genotype)] = "TT"
dat$Genotype[grep("B.*", dat$Genotype)] = "TC"
dat$Genotype[grep("D.*", dat$Genotype)] = "Unk"

dat2 = dat2[which(dat2$Content == "Unkn"),]
Sample = c("1982", "2139", "2146", "2449", "2597", "2937", "1127", "1451", "1848")
Genotype = c("CC", "TT", "CC", "TC", "CC", "TT", "TT", "TT", "TC")
tmp = data.frame(Sample, Genotype)
dat2 = merge(dat2, tmp, by = "Sample")

## merge both qpcr runs
dat = rbind(dat, dat2)

## drop negatives
dat = dat[which(dat$Content == "Unkn"),]
## how many samples failed?
fail = dat[which(is.na(dat$Cq)),]
table(fail$Target)

## ct efficiency adjustment for calculating relative amount of original template ## formula: P = Po(AE)^raw ct
ae = as.numeric(c(me209, me3var, me4var))
ae = data.frame(ae, c("TLN2_209", "TLN2_3var", "TLN2_4var"))
names(ae) = c("ae", "Target")
dat = merge(dat, ae, by = "Target")
dat = dat[-which(dat$Sample %in% fail$Sample),]
## efficiency adjustment not used given high performance of all three primer sets

dat$adj_ct = 2^-dat$Cq
# dat$adj_ct = dat$ae^-dat$Cq
list = unique(dat$Sample)
dct = data.frame()
for (i in c("TLN2_3var", "TLN2_4var")){
  a = dat[which(dat$Target == i),]
  ref = dat[which(dat$Target == "TLN2_209"),]
  for (j in list){
    b = a[which(a$Sample == j),]
    gen = b$Genotype
    b = b$adj_ct
    b_ref = ref[which(a$Sample == j),]
    b_ref = b_ref$adj_ct
    delta = b/b_ref
    tmp = data.frame(i, j, delta, gen)
    dct = rbind(dct, tmp)
  }
}
names(dct) = c("Target", "Sample", "Delta_Ct", "Genotype")

### qpcr statistical analysis
dct = dct[which(dct$Genotype != "Unk"),]

## test for batch effect
dct2 = dct
dct2$Batch = "A"
Sample = c("1982", "2597", "2937", "1848")
for (i in 1:length(dct2$Batch)){
  if(dct2$Sample[i] %in% Sample){
    dct2$Batch[i] = "B"
  }
}
## batch effect negative

m = aov(log(Delta_Ct) ~ Genotype*Target + Batch + Sample, dct2, na.action = na.omit)
summary(m)
## pairwise tests of each combination of a
a = c("Genotype", "Target", "Genotype:Target")
t = TukeyHSD(m, a)
for (i in a){
  j = gsub(":", "_", i)
  b = paste("tukey_", tolower(j), "_all.csv", sep = "")
  write.csv(file = b, t[i])
}
# contrasts(dct2$Genotype:dct2$Target)
# contrasts(dct2$Genotype:dct2$Target) = cbind(c(0.5, -0.5, 0, 0, 0, 0, 0, 0),
#                                              c(0, 0, 0, 0, 0, 0, 0, 0),
#                                              c(0, 0, 0, 0, 0, 0, 0, 0),
#                                              c(0, 0, 0, 0, 0.5, -0.5, 0, 0),
#                                              c(0, 0, 0, 0, 0.5, -0.5, 0, 0),
#                                              c(0, 0, 0, 0, 0, 0, 0.5, -0.5),
#                                              c(0, 0, 0, 0, 0, 0, 0.5, -0.5))

## manual pairwise comparisons between 3var/4var within each genotype
dct2$GenTar = paste(dct2$Genotype, dct2$Target, sep = "_")
gen1 = dct2[which(dct2$Genotype == "TT"),]
gen2 = dct2[which(dct2$Genotype == "TC"),]
gen3 = dct2[which(dct2$Genotype == "CC"),]
gen1.3 = gen1[which(gen1$Target == "TLN2_3var"),]
gen2.3 = gen2[which(gen2$Target == "TLN2_3var"),]
gen3.3 = gen3[which(gen3$Target == "TLN2_3var"),]
gen1.4 = gen1[which(gen1$Target == "TLN2_4var"),]
gen2.4 = gen2[which(gen2$Target == "TLN2_4var"),]
gen3.4 = gen3[which(gen3$Target == "TLN2_4var"),]
t.test(x = log(gen1.3$Delta_Ct), y = log(gen1.4$Delta_Ct), paired = T, var.equal = T, alternative = "less")
t.test(x = log(gen2.3$Delta_Ct), y = log(gen2.4$Delta_Ct), paired = T, var.equal = T, alternative = "greater")
t.test(x = log(gen3.3$Delta_Ct), y = log(gen3.4$Delta_Ct), paired = T, var.equal = T, alternative = "less")

## https://rpubs.com/aaronsc32/post-hoc-analysis-tukey
##
dct$assay = dct$Target
dct$assay = gsub("TLN2_3var", "3 Isoforms", dct$assay)
dct$assay = gsub("TLN2_4var", "4 Isoforms", dct$assay)

dct$logratio = log(dct$Delta_Ct)
dct$Target = NULL
dct$assay = as.factor(dct$assay)
p = ggplot(dct, aes(Genotype, log(Delta_Ct)+4.6, fill = assay)) + 
  geom_hline(yintercept = 6.75, alpha = 0.4, linetype = "dashed") +
  geom_boxplot(aes(fill = assay) ,outlier.shape = NA, alpha = 0.5) + 
  labs(x = "Genotype", y = "Alternative Isoform logratio", fill = "Assay") + 
  geom_bracket(label = c("**", "*", "ns"), y.position = c(7.5, 8.5, 8), xmin = c("CC", "CC", "TC"), xmax = c("TC", "TT", "TT"),
               inherit.aes = F, label.size = 3.25) +
  geom_bracket(label = c("ns", "ns", "*"), y.position = c(6, 6, 5.75), xmin = c("TT", "TC", "CC"), xmax = c("TT", "TC", "CC"),
               inherit.aes = F, label.size = 3.25) +
  geom_point(size = 2, alpha = 0.75, position = position_jitterdodge(jitter.width = 0.15, jitter.height = 0)) +
  cbbFill2 + cbbColour2 +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "bottom", axis.title = element_text(size = 8), axis.text = element_text(size = 8),
        legend.text = element_text(size = 8), legend.title = element_text(size = 9)) +
  expand_limits(y = c(0, 8.75)) 

p1 = p
# png(file = "qpcr_res_prelim.png", res = 600, width = 3, height = 2, units = 'in')
# print(p)
# dev.off()

#####
# go term figure
#####
go = read.csv(file = 'S3_Table.csv', header = T)
go$name = paste(go$Pathway.Description, " (", go$Observed.Gene.Count, ")", sep = "")
go$name = ordered(go$name, levels = rev(go$name))

go1 = go[which(go$Queried.Database == "KEGG"),]
go1 = go1[1:10,]
go1$tln = "TLN2(-)"
go1$tln[grep("TLN2", go1$Matching.Protein.Labels)] = "TLN2(+)"

p = ggplot(go1, aes(name, -log10(False.Discovery.Rate))) +
  geom_bar(aes(fill = tln) ,stat = "identity") + coord_flip() + 
  labs(x = "", y = "-10log(FDR)", fill = "") + #theme(aspect.ratio = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  cbbFill2 +
  scale_y_continuous(limits = c(0, 3) ,expand = c(0, 0), breaks = seq(0, 3, 1)) +
  ggtitle(label = "KEGG") + theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 6.5),
        legend.position = "none")
p5 = p

go2 = go[which(go$Queried.Database == "Component"),]
# go1 = go1[1:10,]
go2$name = as.character(go2$name)
go2$name[1] = "junctional sarcoplasmic\nreticulum membrane (3)"
go2$name = ordered(go2$name, levels = rev(go2$name))
go2$tln = "TLN2(-)"
go2$tln[grep("TLN2", go2$Matching.Protein.Labels)] = "TLN2(+)"

p = ggplot(go2, aes(name, -log10(False.Discovery.Rate))) +
  geom_bar(aes(fill = tln) ,stat = "identity") + coord_flip() + 
  labs(x = "", y = "-10log(FDR)", fill = "") + #theme(aspect.ratio = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  cbbFill2 +
  scale_y_continuous(limits = c(0, 5) ,expand = c(0, 0), breaks = seq(0, 5, 1)) +
  ggtitle(label = "Cellular Component") + theme(plot.title = element_text(hjust = 0.5, size = 8)) +
  theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7),
        legend.position = "none")
p6 = p

# ## ended up dropping due to space limitations - still reporting in S3_Table
# go3 = go[which(go$Queried.Database == "Process"),]
# go3 = go3[1:10,]
# 
# p = ggplot(go3, aes(name, -log10(False.Discovery.Rate))) +
#   geom_bar(stat = "identity", fill = "black") + coord_flip() + 
#   labs(x = "", y = "-10log(FDR)") + #theme(aspect.ratio = 0.5) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), 
#         axis.line = element_line(colour = "black")) +
#   scale_y_continuous(limits = c(0, 5) ,expand = c(0, 0), breaks = seq(0, 5, 1)) +
#   ggtitle(label = "Biological Process") + theme(plot.title = element_text(hjust = 0.5, size = 8)) +
#   theme(axis.title = element_text(size = 7), axis.text = element_text(size = 7))
# p7 = p

##### 
# test out gvis to make better qpcr diagram
#####
## useful guide: https://bioconductor.org/packages/release/bioc/vignettes/Gviz/inst/doc/Gviz.html#46_biomartgeneregiontrack
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Gviz")
library(Gviz)
library(biomaRt)
library(grid)
library(ggplotify) ## needed to convert gvis object into ggplot compatible object

map = read.csv(file = "talin2_domain_mapping.csv", header = T)

bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", chromosome=15, start=62.39e6, end=62.85e6,
                                    name="ENSEMBL", biomart=bm)
plotTracks(biomTrack)
plotTracks(biomTrack, col.line = NULL, col = NULL)

bm <- useMart(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
biomTrack <- BiomartGeneRegionTrack(genome="hg19", symbol = "TLN2",
                                    name="ENSEMBL", biomart=bm)
plotTracks(biomTrack, transcriptAnnotation="symbol")
## access useful internals: biomTrack@range
biomTrack <- BiomartGeneRegionTrack(genome = "GRCh38", transcript = c("ENST00000561311", "ENST00000636159", "ENST00000494733",
                                                              "ENST00000472902"), 
                                    name="ENSEMBL", biomart=bm)
biomTrack2 = biomTrack
biomTrack@range = biomTrack@range[which(biomTrack@range$feature != 'non_coding'),]
biomTrack@range = biomTrack@range[which(biomTrack@range$feature != 'antisense'),]
biomTrack@range = biomTrack@range[which(biomTrack@range$feature != 'miRNA'),]
biomTrack@range = biomTrack@range[which(biomTrack@range$transcript != 'ENST00000558940'),]
biomTrack@range = biomTrack@range[which(biomTrack@range$transcript != 'ENST00000560347'),]

diff = abs(map$Start[which.min(map$Start)] - 62853564)
map$Start = map$Start + diff
map$End = map$End + diff
map$Domain = gsub("N-terminal Head", "Head Domain", map$Domain)
map = map[which(map$Domain != "PIP1 y90"),]
map = map[which(map$Domain != "PIP2"),]

gtrack = GenomeAxisTrack(biomTrack@range, start = 62750000, end = 63150000, chromosome = 15)
grtrack = GeneRegionTrack(biomTrack@range, transcriptAnnotation="transcript", name = "TLN2\nIsoforms", background.title = 'darkgrey',
                          fill = "grey", cex.title = 0.75)
ht = HighlightTrack(trackList = c(grtrack), start = c(62880000, 63031000, 63070000), width = 4000, chromosome = 15)
anno = AnnotationTrack(start = map$Start, end = map$End, group = map$Domain, chromosome = 15, background.title = "darkgrey",
                       name = "TLN2 Functional\nDomains", lwd = 2, cex.title = 0.75)

pd = as.ggplot(~plotTracks(list(gtrack, ht, anno), groupAnnotation = "group", extend.left = 0.075))

#####
### merge figure using patchwork
#####
library(patchwork)
# https://patchwork.data-imaginist.com/articles/guides/layout.html # helpful guide

## version 2
library(grid)
library(ggplotify)

layout = '
ABBBBBCCCCC
ABBBBBCCCCC
'
patchwork = wrap_plots(A = wrap_elements(grid::textGrob("Enrichment Analysis", rot = 90, gp = gpar(fontsize = 9))),
                       B = p5, C = p6, 
                 design = layout)/wrap_elements(pd)/wrap_plots(p4 + p1) + 
        plot_layout(heights = unit(c(3.80, 6.95, 3.75), c('cm')))
patchwork[[1]] = patchwork[[1]] + plot_layout(tag_level = 'new')
patchwork[[3]] = patchwork[[3]] + theme(panel.spacing = unit(0, "lines"), 
                                        strip.background = element_blank(), strip.placement = "outside")

tiff(file = "Fig4.tiff", res = 600, width = 165, height = 205, units = 'mm', compression = 'lzw')
print(patchwork + plot_annotation(tag_levels = 'a', subtitle = "a"))
dev.off()



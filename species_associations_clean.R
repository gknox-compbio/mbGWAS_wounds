cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbFill2 <- scale_fill_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#ABABAB", "#000080", "#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "thistle4", "yellowgreen", "violetred1"))


#####
# Reading and setting up data
#####
otu <- read.csv("clean_species_table_p1.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu <- otu[complete.cases(otu), ]
otu <- otu / rowSums(otu)

otu2 <- read.csv("clean_species_table_p2.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu2 <- otu2[complete.cases(otu2), ]
## making it so tables can be easily merged. Dropping rare and uncommon species ok for these analyses
otu <- otu[, which(names(otu) %in% names(otu2))]
otu2 <- otu2[, which(names(otu2) %in% names(otu))]
otu <- otu[, order(names(otu))]
otu2 <- otu2[, order(names(otu2))]
otu <- rbind(otu, otu2)
rm(otu2)
## counting occurrences for 10% exclusion
temp <- otu
temp[temp > 0] <- 1
keep <- names(which(colSums(temp) >= (0.1 * nrow(temp)))) # ?keep columns where the sum of a column (number of appearances) is greater than 10% the # of rows (total samples)?
rm(temp)
otu <- otu[, keep]
## exclude individuals who do lost all species
otu <- otu[rowSums(otu) > 0, ]

## read in met and allele calls, merge
gen <- read.csv("allele_calls_p1.csv", stringsAsFactors = FALSE, row.names = 2)
gen <- gen[, -1]
gen <- gen[complete.cases(gen), ]
met <- read.csv("met1.csv")
row.names(met) <- met$SampleCode
gen <- merge(gen, met, by = "row.names")
gen$RacialCategory <- gsub("k ", "k", gen$RacialCategory)

## now plate 2
gen2 <- read.csv("allele_calls_p2.csv", stringsAsFactors = FALSE, row.names = 2)
gen2 <- gen2[, -1]
gen2 <- gen2[complete.cases(gen2), ]

met2 <- read.csv("met2.csv")
row.names(met2) <- met2$SampleCode
gen2 <- merge(gen2, met2, by = "row.names")

met <- rbind(gen, gen2)
temp <- as.data.frame(table(met$chartno))
drop <- temp[temp$Freq > 1, ]
met <- met[-which(met$chartno == drop$Var1), ]
row.names(met) <- met$chartno

otu <- otu[which(row.names(otu) %in% row.names(met)), ]
met <- met[which(row.names(met) %in% row.names(otu)), ]
met$Row.names <- NULL

## bring in eigen vectors from SNP PCA
tt <- read.delim("cov3.txt", sep = " ", stringsAsFactors = FALSE)
tt2 <- read.delim("p2_cov3.txt", sep = " ", stringsAsFactors = FALSE)
tt <- tt[, match(names(tt2), names(tt))]
tt <- tt[-c(grep("Control", tt$FID)), ]
tt <- rbind(tt, tt2)
tt$FID <- NULL
names(tt)[1] <- "SampleCode"
rm(tt2)
met <- merge(met, tt, by = "SampleCode")
row.names(met) <- met$chartno

## drop sample with missing znf521 genotype
met <- met[which(met$rs7236481 != "00"), ]
otu <- otu[which(row.names(otu) %in% row.names(met)), ]

#####
# Alpha Diversity plot
#####
library(reshape2)
library(ggpubr)
mtemp <- met
mtemp <- mtemp[, -c(2:7)]
mtemp <- melt(mtemp, measure.vars = 2:3)
names(mtemp)[24:25] <- c("SNP", "Genotype")
mtemp$Genotype <- ordered(mtemp$Genotype, levels = c("CC", "TC", "TT", "GT", "GG"))
p <- ggplot(mtemp, aes(Genotype, Hill, group = Genotype)) +
  geom_boxplot(aes(fill = Genotype), outlier.shape = NA) +
  facet_grid(~SNP, scales = "free_x") +
  labs(x = "Genotype per SNP", y = expression(Hill[1])) +
  geom_jitter(size = 1.5, alpha = 0.5, height = 0, width = 0.2) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none", strip.text = element_text(face = "bold", size = 10), axis.title = element_text(size = 11), axis.text = element_text(size = 10)) +
  scale_fill_manual(values = cb)

# tiff(file = "boxplot_2snp_alpha.tiff", res = 600, width = 3, height = 1.75, units = 'in')
# print(p)
# dev.off()
p4 <- p
p4

#####
# Species association analysis
#####
dat <- merge(met, otu, by = "row.names")

a <- manova(cbind(Pseudomonas.aeruginosa, Staphylococcus.epidermidis, Anaerococcus.vaginalis, Corynebacterium.striatum, Corynebacterium.tuberculostearicum, Escherichia.coli, Finegoldia.magna, Staphylococcus.aureus, Staphylococcus.lugdunensis, Streptococcus.agalactiae)
~ (rs8031916) + (rs7236481) + (Age) + (Diabetes) + (Sex) + (EV.1) + (EV.2) + (EV.3) + (EV.4) + (EV.5), dat)
summary(a)
b <- summary.aov(a)
# sink("manova_taxa_associations_all.txt")
# print(b)
# sink()
# Sig: tln2 = Pa (p = 0.036), Se (0.034); znf521 = Finegoldia (p = 0.069), Staph lug (p = 0.058)

## repeat with 3 znf521 samples dropped
dat <- dat[-which(dat$rs7236481 == "GG"), ]

a <- manova(cbind(Pseudomonas.aeruginosa, Staphylococcus.epidermidis, Anaerococcus.vaginalis, Corynebacterium.striatum, Corynebacterium.tuberculostearicum, Escherichia.coli, Finegoldia.magna, Staphylococcus.aureus, Staphylococcus.lugdunensis, Streptococcus.agalactiae)
~ (rs8031916) + (rs7236481) + (Age) + (Diabetes) + (Sex) + (EV.1) + (EV.2) + (EV.3) + (EV.4) + (EV.5), dat)
summary(a)
b <- summary.aov(a)
# sink("manova_taxa_associations_cut.txt")
# print(b)
# sink()

# ## plot residual abundance of Se and Pa for Tln2
keep <- c("Staphylococcus.epidermidis", "Pseudomonas.aeruginosa")
temp <- dat
temp <- temp[, -c(33:38, 40, 42, 43)]
temp <- melt(temp, measure.vars = 33:34)
names(temp)[33:34] <- c("Species", "Abundance")

# Extract residual abundance for each identified species
res <- data.frame()
for (i in keep) {
  rtemp <- temp[which(temp$Species == i), ]
  lmres <- lm(Abundance ~ Age + Diabetes + Sex, rtemp)
  rtemp$res <- lmres$residuals
  res <- rbind(res, rtemp)
}

res$Species <- gsub("\\.", " ", res$Species)
p <- ggplot(res, aes(rs8031916, res, group = rs8031916)) +
  geom_boxplot(aes(fill = rs8031916), outlier.shape = NA) +
  facet_grid(~Species, labeller = label_wrap_gen(multi_line = T, width = 15)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  geom_jitter(size = 1.5, height = 0, width = 0.2, alpha = 0.5) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none", strip.text = element_text(face = "bold.italic", size = 10)) +
  labs(x = "Genotype at rs8031916", y = "Residual Relative Abundance") +
  scale_fill_manual(values = cb)
# tiff(file = "boxplot_taxon_res.tiff", res = 600, width = 3, height = 1.75, units = 'in')
# print(p)
# dev.off()
p5 <- p
p5

## Plot using relative abundance instead of residual relative abundance
# p = ggplot(temp, aes(rs8031916, Abundance, group = rs8031916)) + geom_boxplot(aes(fill = rs8031916), outlier.size = NA) + facet_grid(~Species) +
#   theme(axis.text = element_text(size = 4.5), axis.title = element_text(size = 6)) + geom_jitter(size =1.25, height = 0, width = 0.2, alpha = 0.5) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
#   theme(legend.position = "none", strip.text = element_text(face = "bold.italic", size = 5)) +
#   labs(x = "Genotype at rs8031916", y = "Relative abundance")
# tiff(file = "boxplot_taxon_ra.tiff", res = 600, width = 3, height = 1.75, units = 'in')
# print(p)
# dev.off()
# p5_ra = p

# write.csv(met, file = "all_samples_with_genotypes.csv", row.names = F)

#####
# Beta diversity tests
#####
met <- met[order(as.numeric(as.character(row.names(met)))), ]
otu <- otu[order(as.numeric(as.character(row.names(otu)))), ]

## Ordination/Adonis using genotypes
bc <- vegdist(otu, method = "bray")
adonis(bc ~ met$rs8031916 + met$rs7236481 + met$Age + met$Diabetes)
# adonis(bc ~ met$rs8031916 * met$rs7236481 * met$Age * met$Diabetes) # no sig interactions

# ## quick retest of adonis after dropping znf521 gg individuals # NO DIFFERENCE
# met2 = met[-which(met$rs7236481 == "GG"),]
# bc = as.matrix(bc)
# list = met2$chartno
# bc2 = bc[which(row.names(bc) %in% list), which(row.names(bc) %in% list)]
#
# adonis(bc2 ~ met2$rs8031916 + met2$rs7236481 + met2$Age + met2$Diabetes)
# # adonis(bc2 ~ met2$rs8031916 * met2$rs7236481 * met2$Age * met2$Diabetes) # no sig interactions

d <- capscale(bc ~ 1)
dat <- data.frame(MDS1 = d$CA$u[, 1], MDS2 = d$CA$u[, 2])
dat$chartno <- row.names(dat)
dat <- merge(dat, met, by = "chartno")
ar <- dat[, 2:3]
s <- summary(eigenvals(d))
### species biplot
library(ggrepel)
library(stringr)
c <- t(cor(ar, otu)) # correlation of axes with variables of interest
c <- as.data.frame(na.omit(c))
dists <- sqrt(c[, 1]^2 + c[, 2]^2) # see which genera show maximal combined correlation with axes
keep <- order(dists, decreasing = T)[1:5] # keep the top 5
postn <- rep(0.5, length(keep)) # this is used to offset the labels- if x cor is negative, the hjust = 1, else hjust = 0
postn[c[keep, 1] < 0] <- 1
postn[c[keep, 1] > 0] <- 0
# xend and yend are multiplies by a fraction that is the scaling factor you'd like for plotting
d.arrows <- data.frame(
  x = rep(0, length(keep)), y = rep(0, length(keep)), xend = c[keep, 1] * .5, yend = c[keep, 2] * .5,
  labels = rownames(c)[keep], p = postn
)
d.arrows$labels <- gsub("\\;.*", "", d.arrows$labels)
d.arrows$labels <- gsub("\\.", " ", d.arrows$labels)

## ordination plot
p <- ggplot(dat, aes(MDS1, MDS2)) +
  geom_point(aes(color = rs8031916), size = 2.5, alpha = 0.6) +
  scale_color_manual(values = cb)
p <- p + labs(
  x = paste("MDS1 (", round(s[2, 1] * 100, digits = 2), "%)", sep = ""),
  y = paste("MDS2 (", round(s[2, 2] * 100, digits = 2), "%)", sep = "")
)
p <- p + annotate("segment",
  x = d.arrows$x, xend = 1.1 * d.arrows$xend, y = d.arrows$y,
  yend = 1.1 * d.arrows$yend, size = 1.15, alpha = 0.6, arrow = arrow(length = unit(0.2, "cm"))
)
p <- p + geom_text_repel(
  data = d.arrows, aes(
    x = 1.1 * xend, y = 1.1 * yend,
    label = str_wrap(labels)
  ), size = 3, point.padding = 0.5,
  fontface = "italic"
) +
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme(legend.position = "top", axis.text = element_text(size = 9), axis.title = element_text(size = 10)) +
  theme(legend.title = element_text(size = 10), legend.text = element_text(size = 9))
p12 <- p
p12

## accompanying boxplots per axis
p1 <- ggplot(dat, aes(rs8031916, MDS1, fill = rs8031916)) +
  geom_boxplot(aes(fill = rs8031916)) +
  scale_fill_manual(values = cb) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip() +
  labs(x = "", y = "") +
  theme(legend.position = "none", axis.text = element_blank())
p2 <- ggplot(dat, aes(rs8031916, MDS2, fill = rs8031916)) +
  geom_boxplot(aes(fill = rs8031916)) +
  scale_fill_manual(values = cb) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(x = "", y = "") +
  theme(legend.position = "none", axis.text = element_blank())

# library(grid)
library(gridExtra)
# library(ggpubr)
tiff("bray_betadiv_pcoa.tiff", res = 600, width = 5.5, height = 5.5, units = "in")
print(grid.arrange(p, p1, p2,
  layout_matrix = rbind(
    c(NA, 1, 1, 1, 1, 1, 1),
    c(3, 1, 1, 1, 1, 1, 1),
    c(3, 1, 1, 1, 1, 1, 1),
    c(3, 1, 1, 1, 1, 1, 1),
    c(3, 1, 1, 1, 1, 1, 1),
    c(NA, 2, 2, 2, 2, 2, 2)
  )
))
dev.off()

#####
# data set up for Se/Pa community comparisons
#####
o <- read.csv("clean_species_table_p1.csv", header = TRUE, row.names = 1)
o <- o[complete.cases(o), ]

o <- o / rowSums(o)
o$SampleCode <- row.names(o)
o <- melt(o, id.vars = "SampleCode")

o2 <- read.csv("clean_species_table_p2.csv", header = TRUE, row.names = 1)
o2 <- o2[complete.cases(o2), ]
o2$SampleCode <- row.names(o2)
o2 <- melt(o2, id.vars = "SampleCode")

otu <- rbind(o, o2)
rm(o, o2)
otu <- dcast(data = otu, SampleCode ~ variable, fun.aggregate = sum)
row.names(otu) <- otu[, 1]
otu <- otu[, -1]

#####
# how does pa and staph sp. community diversity compare?
#####

# pull out communities with se and pa
s <- row.names(otu)[which(otu$Staphylococcus.epidermidis > 0)]
pa <- row.names(otu)[which(otu$Pseudomonas.aeruginosa > 0)]

# calculate diversity metrics
one <- exp(vegan::diversity(otu[s, ], index = "shannon"))
four <- exp(vegan::diversity(otu[pa, ], index = "shannon"))
hill <- data.frame(
  rich = c(one, four),
  group = c(
    rep("Se", length(one)),
    rep("Pa", length(four))
  )
)
hill$group <- ordered(hill$group, levels = c("Pa", "Se"))

# p-value annotation from t-test below
p <- ggplot(hill, aes(group, rich, group = group)) +
  geom_boxplot(aes(fill = group), outlier.size = NA) +
  geom_jitter(height = 0, width = 0.2, size = 2.5, alpha = 0.5) +
  ylim(1, 8) +
  labs(x = "Community", y = expression(Hill[1])) +
  theme_bw(base_size = 15) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_fill_manual(values = cb) +
  annotate("text", label = "______________", x = 1.5, y = 7.79, size = 4.25, color = "black") +
  annotate("text", label = "*", x = 1.5, y = 7.74, size = 4.5, color = "black")

# tiff(file = "community_diversity_pa_se.tiff", height = 1.75, width = 2, units = "in", res = 600)
# print(p)
# dev.off()
p6 <- p
p6

shapiro.test(hill$rich)
skewness(hill$rich)

t.test(hill$rich ~ hill$group, alternative = "less")
## welches t-test results
# t = 2.0883, df = 49.3, p = 0.021
# means: pa = 2.5197, se = 3.27309

#####
# merged figures
#####

library(grid)
library(gridExtra)
library(ggpubr)
gl <- lapply(1:3, function(ii) grobTree(rectGrob(), textGrob(ii)))

## Fig 1 with p4, p12 (p1 and p2 accompanying), p5 and p6
tiff("Figure_1.tiff", res = 600, width = 180, height = 185, units = "mm")
print(grid.arrange(arrangeGrob(p4, top = text_grob("a", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
  arrangeGrob(p12, top = text_grob("b", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
  arrangeGrob(p1),
  arrangeGrob(p2),
  arrangeGrob(p5, top = text_grob("c", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), size = 12)),
  arrangeGrob(p6, top = text_grob("d", x = unit(0.025, "npc"), y = unit(1, "npc"), just = c("left", "top"), size = 12)),
  layout_matrix = rbind(
    c(1, 1, 1, 1, NA, 2, 2, 2, 2, 2),
    c(1, 1, 1, 1, 4, 2, 2, 2, 2, 2),
    c(1, 1, 1, 1, 4, 2, 2, 2, 2, 2),
    c(1, 1, 1, 1, 4, 2, 2, 2, 2, 2),
    c(1, 1, 1, 1, 4, 2, 2, 2, 2, 2),
    c(NA, NA, NA, NA, NA, 3, 3, 3, 3, 3),
    c(5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
    c(5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
    c(5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
    c(5, 5, 5, 5, 5, 6, 6, 6, 6, 6),
    c(5, 5, 5, 5, 5, 6, 6, 6, 6, 6)
  )
))
dev.off()

#####
# Reading and setting up data for network correlations
#####
otu <- read.csv("clean_species_table_p1.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu <- otu[complete.cases(otu), ]
otu <- otu / rowSums(otu)

otu2 <- read.csv("clean_species_table_p2.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu2 <- otu2[complete.cases(otu2), ]
## making it so tables can be easily merged. Dropping rare and uncommon species ok for these analyses
otu <- otu[, which(names(otu) %in% names(otu2))]
otu2 <- otu2[, which(names(otu2) %in% names(otu))]
otu <- otu[, order(names(otu))]
otu2 <- otu2[, order(names(otu2))]
otu <- rbind(otu, otu2)
rm(otu2)
## counting occurrences for 5% exclusion
temp <- otu
temp[temp > 0] <- 1
keep <- names(which(colSums(temp) >= (0.05 * nrow(temp))))
rm(temp)
otu <- otu[, keep]
## exclude individuals who lost all species
otu <- otu[rowSums(otu) > 0, ]

otu <- otu[order(as.numeric(row.names(otu))), ]

#####
# geom_net correlation figures
#####
library(geomnet)

keep3 <- data.frame()
for (i in unique(names(otu))) {
  if (length(which(otu[, i] > 0)) > 0.05 * dim(otu)[1]) {
    keep1 <- "yes"
  } else {
    keep1 <- "no"
  }
  occurrences <- length(which(otu[, i] > 0))
  keep2 <- data.frame(i, keep1, occurrences)
  names(keep2) <- c("abs_ids", "keep", "Observed")
  keep3 <- rbind(keep3, keep2)
}
keep <- keep3[which(keep3$keep == "yes"), ]
temp <- otu[, keep$abs_ids]
temp <- temp[rowSums(temp) > 0, ]

cor.matrix <- cor(temp, method = "pearson")
# Convert correlation matrix to binary adjacency matrix
cor.cutoff <- 0.1
cor.adj <- ifelse(abs(cor.matrix) >= cor.cutoff, 1, 0)
features <- row.names(cor.adj)
rownames(cor.adj) <- as.character(1:nrow(cor.adj))
colnames(cor.adj) <- as.character(1:ncol(cor.adj))

flattenSquareMatrix <- function(m) {
  if ((class(m) != "matrix") | (nrow(m) != ncol(m))) stop("Must be a square matrix.")
  if (!identical(rownames(m), colnames(m))) stop("Row and column names must be equal.")
  ut <- upper.tri(m)
  data.frame(
    i = rownames(m)[row(m)[ut]],
    j = rownames(m)[col(m)[ut]],
    cor = t(m)[ut],
    p = m[ut]
  )
}

fa <- flattenSquareMatrix(cor.adj)
fa <- fa[, -4]
bin <- as.data.frame(cor.matrix)
names(bin) <- 1:length(names(bin))
list <- data.frame(names(bin), row.names(bin))
row.names(bin) <- 1:length(row.names(bin))
fc <- flattenSquareMatrix(as.matrix(bin))
net <- cbind(fa, fc[3])
# net = net[,-4]
names(net)[3] <- "plot"

net3 <- net[which(net$plot == 1), ]
names(net3)[1:2] <- c("from_id", "to_id")
net3$dir <- "pos"
net3$dir[which(net3$cor <= 0)] <- "neg"
net3$col <- "goldenrod"
net3$col[which(net3$cor <= 0)] <- "black"
a <- setdiff(unique(net3$to_id), unique(net3$from_id))
temp <- data.frame(a, NA, NA, NA, NA, NA)
names(temp) <- names(net3)
net3 <- rbind(net3, temp)

net3$from_id <- as.character(net3$from_id)
net3$to_id <- as.character(net3$to_id)

net3 <- net3[order(net3$from_id), ]
## numbers were dropped during pearson correlation filtering, need to update number list
net3$from_id2 <- net3$from_id
# net3$to_id2 = net3$to_id
net3$from_id2 <- as.numeric(as.character(net3$from_id2))
# net3$to_id2 = as.numeric(as.character(net3$to_id2))
net3 <- net3[order(net3$from_id2), ]
new <- seq(1, length(unique(net3$from_id)), 1)
drop <- setdiff(new, net3$from_id2) # use these numbers to drop species from list below
old <- unique(net3$from_id2)
new <- data.frame(new, old)
names(new)[2] <- "from_id"
net3 <- merge(net3, new, by = "from_id")
net3$from_id2 <- NULL
names(new)[2] <- "to_id"
names(net3)[length(names(net3))] <- "from_id2"
net3 <- merge(net3, new, by = "to_id", all.x = T)
names(net3)[length(names(net3))] <- "to_id2"

net3$from_id2 <- as.character(net3$from_id2)
net3$to_id2 <- as.character(net3$to_id2)

net3 <- net3[order(net3$from_id2), ]

p <- ggplot(data = net3, aes(from_id = from_id2, to_id = to_id2)) +
  theme_net() +
  geom_net(
    directed = FALSE, labelon = TRUE, repel = FALSE, ecolour = net3$col, labelcolour = "white", size = 10,
    layout.alg = "fruchtermanreingold", linewidth = 1.25, ealpha = 0.6, fontsize = 3.5
  )
p7 <- p

## for making identifier table
list <- list[-drop, ]
list$names.bin. <- seq(1, length(unique(net3$from_id)), 1)
list$Species <- paste(list$names.bin., ") ", list$row.names.bin., sep = "")
names(list)[1:2] <- c("num", "sp")
# list$letter = letters[sort(as.numeric(list$num))]

## connectivity plot
dat.node <- as.data.frame(table(net3$from_id2, net3$col))
names(dat.node) <- c("Id", "cor", "freq")
dat.node$merge <- paste(dat.node$Id, dat.node$cor, sep = "_")
dat.node2 <- as.data.frame(table(net3$to_id2, net3$col))
names(dat.node2) <- c("Id", "cor", "freq2")
dat.node2$merge <- paste(dat.node2$Id, dat.node2$cor, sep = "_")
dat.node <- merge(dat.node, dat.node2, by = "merge", all = T, no.dups = T)
rm(dat.node2)
dat.node$freq2[is.na(dat.node$freq2)] <- 0
dat.node$count <- dat.node$freq + dat.node$freq2
dat.node <- dat.node[, c(1:3, 8)]
names(dat.node) <- c("merge", "Id", "Cor", "Count")
for (i in 1:length(dat.node$merge)) {
  if (dat.node$Cor[i] == "black") {
    dat.node$Count[i] <- dat.node$Count[i] * -1
  }
}
dat.node$Cor <- gsub("goldenrod", "Positive", dat.node$Cor)
dat.node$Cor <- gsub("black", "Negative", dat.node$Cor)
dat.node <- dat.node[order(as.numeric(as.character(dat.node$Id))), ]
dat.node$Id <- ordered(dat.node$Id, levels = sort(as.numeric(as.character(unique(dat.node$Id)))))
p <- ggplot(dat.node, aes(Id, Count, fill = Cor)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = "Correlation Counts", x = "") +
  theme(axis.text = element_text(size = 9), legend.text = element_text(size = 11)) +
  cbbFill2 +
  guides(fill = guide_legend(title = "Correlation")) +
  theme(legend.title = element_text(size = 12), axis.title.x = element_text(size = 12))
p8 <- p

## connectivity plot version 2 -- using sum pearson coefficient
vec <- unique(c(net3$from_id2, net3$to_id2))
vec <- vec[!is.na(vec)]
dat.nod <- data.frame()
for (i in vec) {
  for (j in unique(net3$dir)) {
    a <- net3[which(net3$dir == j), ]
    a1 <- a[which(a$from_id2 == i), ]
    a2 <- a[which(a$to_id2 == i), ]
    s1 <- sum(a1$cor)
    s2 <- sum(a2$cor)
    s <- s1 + s2
    temp <- data.frame(i, j, s)
    dat.nod <- rbind(dat.nod, temp)
  }
}
## nas introduced by net3$dir
dat.nod <- dat.nod[!is.na(dat.nod$j), ]
names(dat.nod) <- c("Id", "cor", "sum")

dat.nod$cor <- gsub("pos", "Positive", dat.nod$cor)
dat.nod$cor <- gsub("neg", "Negative", dat.nod$cor)
# dat.nod = dat.nod[order(as.numeric(as.character(dat.nod$Id))),]
dat.nod$Id <- ordered(dat.nod$Id, levels = sort(as.numeric(as.character(unique(dat.nod$Id)))))
p <- ggplot(dat.nod, aes(Id, sum, fill = cor)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  labs(y = "Correlation Strength\n(Sum Pearson's Coefficient)", x = "") +
  theme(axis.text = element_text(size = 7.5), legend.text = element_text(size = 7.5)) +
  cbbFill2 +
  guides(fill = guide_legend(title = "Correlation")) +
  theme(legend.title = element_text(size = 9), axis.title.x = element_text(size = 9)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p9 <- p
## turn labels vector into dataframe for plotting with grid.table function - apologies for the ugly code
plist <- list$Species
plist <- paste(plist, "   ", sep = "")
a <- ceiling(length(plist) / 3) # number is number of columns desired
b <- 3 * a
bn <- b - length(plist)
plist[(b - bn + 1):b] <- "" # use if bn does not equal zero, comment out otherwise
# tab = data.frame(plist[1:a], plist[(a+1):(2*a)], plist[(2*a+1):(3*a)], plist[(3*a+1):(4*a)], plist[(4*a+1):length(plist)])
# tab = data.frame(plist[1:a], plist[(a+1):(2*a)], plist[(2*a+1):(3*a)], plist[(3*a+1):length(plist)]) # made 4 column table
tab <- data.frame(plist[1:a], plist[(a + 1):(2 * a)], plist[(2 * a + 1):length(plist)]) # made 3 column table

# names(tab) = paste("v", 1:5, sep = "")
# tab$v5 = as.character(tab$v5)
# plist[(b-bn+1):b] = ""
# library(gridExtra)
# t = tableGrob(tab, theme = ttheme_minimal(core = list(fg_params = list(hjust = 0, x = 0.1, fontface = 3)),
#                                           colhead = list(fg_params = list(parse = T))), cols = NULL, rows = NULL)
# t2 = t
# t2$widths = unit(rep(1/ncol(t2), ncol(t2)), "npc")
# # p1 = grid.arrange(p, p2, t, layout_matrix = rbind(c(1,1,2), c(1,1,2), c(3,3, 3)),
# #                   newpage = F)
# png('pearson_network_species_all.png', res = 600, width = 9, height = 7, units = 'in')
# print(grid.arrange(p7, p8, t, layout_matrix = rbind(c(1,1,2), c(1,1,2), c(3,3, 3)),
#                    newpage = F))
# dev.off()

### two versions of figure 2 - with counts and sums in connectivity plot
library(grid)
library(gridExtra)
library(ggpubr)

# tiff('Figure_2_counts.tiff', res = 600, width = 180, height = 115, units = 'mm')
# print(grid.arrange( arrangeGrob(p7, top = text_grob("a", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
#                     arrangeGrob(p8, top = text_grob("b", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
#                     layout_matrix = rbind(c(1, 1 , 1 ,2,2),
#                                           c(1, 1, 1 ,2 ,2),
#                                           c(1, 1, 1 ,2 ,2)),
#                     newpage = F))
# dev.off()

tiff("Figure_2.tiff", res = 600, width = 180, height = 90, units = "mm", compression = "lzw")
print(grid.arrange(arrangeGrob(p7, top = text_grob("a", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
  arrangeGrob(p9, top = text_grob("b", x = unit(0.025, "npc"), y = unit(0.8, "npc"), just = c("left", "top"), size = 12)),
  layout_matrix = rbind(
    c(1, 1, 1, 2, 2),
    c(1, 1, 1, 2, 2),
    c(1, 1, 1, 2, 2)
  ),
  newpage = F
))
dev.off()

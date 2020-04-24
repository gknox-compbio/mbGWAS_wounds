#####
#setting up the data
#####
temp = read.delim("Biobank with Wound Heal Dates 12-18-2018.csv", header = FALSE, sep = ",", stringsAsFactors = F)
temp = temp[,-c(5)]

#there are just a few that have amputation dates (last column), but these dates are usually strange in comparison to wound start, so dropping
drop = which(temp$V8 != "NULL")
temp = temp[-drop,]
temp = temp[,-c(7)]

#strange character from file type, overwriting
temp[1,1] = 2016

temp[,3] = gsub(" .*", "", temp[,3])
temp[,5] = gsub(" .*", "", temp[,5])
temp[,6] = gsub(" .*", "", temp[,6])

temp[,3] = as.Date(temp[,3])
temp[,5] = as.Date(temp[,5])
temp[,6] = as.Date(temp[,6])

names(temp) = c("SampleCode", "chartno", "pcr.date", "location", "wound.start", "wound.end")

#odd sample
temp = temp[-c(grep(2107, temp$SampleCode)),]

#no wound start date, must drop
temp = temp[!is.na(temp$wound.start),]

temp$status = ifelse(is.na(temp$wound.end), "ongoing", "healed")

temp$duration = temp$wound.end - temp$wound.start

#bring in eigen vectors from SNP PCA
tt = read.delim("cov2.txt", sep = " ", stringsAsFactors = FALSE)
tt2 = read.delim("p2_cov2.txt", sep = "\t", stringsAsFactors = FALSE)
tt = tt[,match(names(tt2), names(tt))]
tt = tt[-c(grep("Control", tt$FID)),]
tt$FID = as.integer(tt$FID)
tt = tt[,-2]
tt2 = tt2[,-2]
tt = rbind(tt, tt2)
rm(tt2)
names(tt)[1] = "SampleCode"

#bring in patient ages
a = read.csv("ages_p1_p2.csv", row.names = 1)
tt = merge(tt, a, by = "SampleCode")
temp = merge(temp, tt, by = "SampleCode")
rm(a, tt)
temp = temp[,-c(grep("chartno", names(temp))[2])]
names(temp)[grep("chartno", names(temp))] = "chartno"

#bring in richness
dat = read.csv("met1.csv", stringsAsFactors = FALSE)
dat2 = read.csv("met2.csv", stringsAsFactors = FALSE)
dat = rbind(dat, dat2)
rm(dat2)

dat = dat[dat$chartno %in% temp$chartno,]

temp = temp[temp$SampleCode %in% dat$SampleCode,]
dat = dat[dat$SampleCode %in% temp$SampleCode,]

#sometimes SampleCodes show up more than once in the wound duration data- indicates multiple wounds that can't be directly linked to microbiome data- need to drop
drop = as.numeric(names(which(table(temp$SampleCode) > 1)))

temp = temp[!temp$SampleCode %in% drop,]
dat = dat[!dat$SampleCode %in% drop,]

look = temp[,7:8]
look$duration = as.numeric(look$duration)
look = look[which(look$status == "healed"),]
summary(look$duration)

dat = merge(dat, temp, by = "chartno")

dat$duration = as.integer(dat$duration)

plot(dat$Hill, dat$duration)

row.names(dat) = dat$SampleCode.x
## pull samples from temp where treatment is ongoing.
req = temp[which(temp$status == "ongoing"),]
gen = read.csv(file = "allele_calls_p1.csv", header = T)
gen2 = read.csv(file = "allele_calls_p2.csv", header = T)
gen$SampleCode = as.character(gen$SampleCode)
gen2$SampleCode = as.character(gen2$SampleCode)
gen = rbind(gen, gen2)
rm(gen2)

row.names(gen) = gen$SampleCode
met = merge(dat, gen, by = "row.names")

#####
# heteroscedasticity
#####
m1 = lm(duration ~ Hill, data = dat)
plot(m1, which = c(1), col = 1, add.smooth = FALSE)

het = dat
het = het[complete.cases(het$duration),]
m1 = lm(duration ~ Hill, data = het)
het$res = c(m1$residuals)

## conduct Breusch-Pagan test for heteroscedasticity
bp = lm(res^2 ~ Hill, data = het)
summary(bp)

p = ggplot(het, aes(Hill, res)) + geom_abline(slope = 0, size = 1.5) + geom_point(col = "blue", size = 2.5) + theme_bw() +
  labs(x = expression(Hill["1"]), y = "Residuals") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  annotate("text", label = "BP = 4.91, p = 0.027", x = 6, y = 1750, size = 3, color = "black")
# tiff(file = "variance_assessment_hill_healing_wls.tiff", res = 600, width = 5, height = 5, units = 'in')
# print(p)
# dev.off()

## actual coded test. (same result)
library(lmtest)
bptest(m1)
## studentized Breush-Pagan test
## bp (test statistic) = 4.91, df = 1, p = 0.0267
# adding to above plot from here
# null hypothesis is homoscedastic variance

m1 = lm(duration ~ Hill + EV.4 + Age, data = het)
bptest(m1)
m1 = lm(duration ~ Age, data = het)
bptest(m1)
m1 = lm(duration ~ EV.4, data = het)
bptest(m1)
## repeating residual plot with multiple regression residuals (variables selected through backward selection)
m1 = lm(duration ~ Hill + EV.4 + Age, data = het)
bptest(m1)
het$mres = c(m1$residuals)
# p = ggplot(het, aes(Hill, mres)) + geom_abline(slope = 0, size = 1.5) + geom_point(col = "blue", size = 2.5) + theme_bw() +
#   labs(x = expression(Hill["1"]), y = "Residuals") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
#   annotate("text", label = "BP = 15.4, p = 0.001", x = 6, y = 1750, size = 3, color = "black")
# # tiff(file = "variance_assessment_multivariate_wls.tiff", res = 600, width = 5, height = 5, units = 'in')
# # print(p)
# # dev.off()

## visualize relationship of variables to variance
ggplot(het, aes(Age, mres)) + geom_abline(slope = 0, size = 1.5) + geom_point(col = "blue", size = 2.5) + theme_bw() +
  labs(x = "Age", y = "Residuals") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  annotate("text", label = "BP = 15.4, p = 0.001", x = 35, y = 1750, size = 3, color = "black")
ggplot(het, aes(EV.4, mres)) + geom_abline(slope = 0, size = 1.5) + geom_point(col = "blue", size = 2.5) + theme_bw() +
  labs(x = "Ev.4", y = "Residuals") + theme(axis.title = element_text(size = 14), axis.text = element_text(size = 12)) +
  annotate("text", label = "BP = 15.4, p = 0.001", x = 0.5, y = 1750, size = 3, color = "black")


#####
#weighted regression
#####
m1 = lm(duration ~ Hill, data = het)
w = m1$residuals
w = 1/abs(w)
summary(lm(duration ~ Hill, data = het, weights = w))

#weighted regression with covariates
m1 = lm(duration ~ Hill + EV.4 + Age, data = het)
w = m1$residuals
w  = 1/abs(w)
## starting model before backward selection
# summary(lm(duration ~ Hill + EV.1  + EV.2 + EV.3 + EV.4 + EV.5 + Sex + Diabetes + Age, data = het, weights = w))

summary(lm(duration ~ Hill + EV.4 + Age, data = het, weights = w))
## Final model involving diversity
library(car)
Anova(lm(duration ~ Hill + EV.4 + Age, data = het, weights = w), type = 2)
1-(((1-(6023.7/31559.3))*(58-1))/(58-3-1)) # hill adjusted R2 = 0.146
1-(((1-(4453.8/31559.3))*(58-1))/(58-3-1)) # ev4 adjusted R2 = 0.093
1-(((1-(1705.4/31559.3))*(58-1))/(58-3-1)) # hill adjusted R2 = 0.001

# coplot(duration~Hill|EV.4*Age, panel = function(x,y,...) panel.smooth(x,y,span = 0.8,...), data = het,
#        xlab = c(expression(Hill[1]), "Conditioning variable: EV.4"), ylab = c("Healing Duration (days)", "Conditioning variable: Age"))

p = ggplot(NULL, aes(Hill, duration)) + geom_point(data = het, size = 2.5, alpha = 0.5) + geom_line(data = temp, color = "blue", size = 2) + 
  theme_bw(base_size = 15) + labs(x = expression(Hill[1]), y = "Healing Duration (days)") +
  annotate("text", label = "F = 5.43, p = 0.0025", x = 6, y = 2250, size = 3, color = "black") +
  annotate("text", label = "R2 = 0.189", x = 6, y = 2150, size = 3, color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
tiff(file = "multivariate_wls.tiff", res = 600, width = 4.5, height = 4.5, units = 'in')
print(p)
dev.off()

## for pub version
p = ggplot(NULL, aes(Hill, duration)) + geom_point(data = het, size = 2.5, alpha = 0.5) + 
  labs(x = expression(Hill[1]), y = "Healing Duration (days)") + 
  geom_smooth(data = het, aes(Hill, duration),method = "lm") +
  annotate("text", label = "F = 5.43, p = 0.0025", x = 6, y = 2300, size = 4, color = "black") +
  annotate("text", label = "R2 = 0.189", x = 6, y = 2150, size = 4, color = "black") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11))

tiff(file = "Figure_3.tiff", res = 600, width = 88, height = 90, units = 'mm')
print(p)
dev.off()

# #####
# # Model selection for Healing ~ Hill - QuasiPoisson or Neg. Binomial -- ultimately did not perform better than linear regression
# #####
# library(MASS)
# 
# m1 = glm(duration ~ Hill + EV.1  + EV.2 + EV.3 + EV.4 + EV.5 + Sex + Diabetes + Age, data = het,
#          family = poisson(link = log))
# summary(m1)
# dispersiontest(m1, alternative = "two.sided")
# ## variance is not equidispersed from mean, qp or nb should be used
# 
# # m2 = glm(duration ~ Hill + EV.1  + EV.2 + EV.3 + EV.4 + EV.5 + Sex + Diabetes + Age, data = het,
# #          family = quasipoisson(link = log))
# m2 = glm(duration ~ Hill + EV.1  + EV.3 + EV.4 + Age, data = het,
#          family = quasipoisson(link = log))
# m2 = glm(duration ~ Hill, data = het,
#          family = quasipoisson(link = log))
# summary(m2)
# # m3 = glm.nb(duration ~ Hill + EV.1  + EV.2 + EV.3 + EV.4 + EV.5 + Sex + Diabetes + Age, data = het)
# m3 = glm.nb(duration ~ Hill + EV.3 + EV.4 + EV.5 + Age, data = het)
# summary(m3)
# het$nb.res = m3$residuals
# het$nb.weights = m3$weights
# # m3 = glm.nb(duration ~ Hill + EV.4 + EV.5 + Age, data = het)
# # summary(m3)
# # m3 = glm.nb(duration ~ Hill, data = het)
# # summary(m3)
# histogram(m3$df.residual)
# plot(m3$weights, abs(m3$residuals))
# plot(m2$weights, abs(m2$residuals))
# plot(m2$linear.predictors, m2$weights)
# plot(m3$linear.predictors, m3$weights)
# plot(m2$data$Hill, m2$weights)
# plot(m3$data$Hill, m3$weights)

##### 
# healing by genotype
#####

# drop missing data and rare rs7236481 genotype
met = met[-which(met$rs7236481 == "00"),]
met2 = met[-which(met$rs7236481 == "GG"),]

summary(aov(duration ~ rs8031916, met))
summary(aov(duration ~ rs8031916*rs7236481, met2))

ggplot(met, aes(rs8031916, duration)) + geom_boxplot() + geom_jitter(height = 0)
ggplot(met, aes(rs7236481, duration)) + geom_boxplot() + geom_jitter(height = 0)
ggplot(met2, aes(rs8031916, duration)) + geom_boxplot() + geom_jitter(height = 0) + facet_grid(~rs7236481)

summary(aov(duration ~ rs8031916 + rs7236481 + Age + Hill, met2))
ggplot(met2, aes(rs8031916, duration)) + geom_boxplot() + geom_jitter(height = 0) + facet_grid(~rs7236481)

table(met$rs8031916, met$status)

#####
# exploratory analysis - are any taxa correlated with healing?
#####
otu = read.csv("clean_species_table_p1.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu = otu[complete.cases(otu),]
otu = otu/rowSums(otu)

otu2 = read.csv("clean_species_table_p2.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu2 = otu2[complete.cases(otu2),]
## making it so tables can be easily merged. Dropping rare and uncommon species ok for these analyses
otu = otu[,which(names(otu) %in% names(otu2))]
otu2 = otu2[,which(names(otu2) %in% names(otu))]
otu = otu[,order(names(otu))] 
otu2 = otu2[,order(names(otu2))] 
otu = rbind(otu, otu2)
rm(otu2)
## counting occurrences for 10% exclusion
temp = otu
temp[temp>0] = 1
keep = names(which(colSums(temp) >= (0.1 * nrow(temp))))
rm(temp)
otu = otu[,keep]
## exclude individuals who do lost all species
otu = otu[rowSums(otu) > 0,]

# otu = otu[which(row.names(otu) %in% lmet$chartno),]
otu = otu[order(as.numeric(row.names(otu))),]

### merge with preexisting 'het' file
row.names(het) = het$chartno
temp = otu
# use 10% correlation cut off to identify taxa to include
temp = merge(het, temp, by = "row.names")

## starting model before backward selection
# m = lm(duration ~ Pseudomonas.aeruginosa + Escherichia.coli + Anaerococcus.lactolyticus + 
#          Anaerococcus.vaginalis + Corynebacterium.striatum + Corynebacterium.tuberculostearicum + Finegoldia.magna +
#          Staphylococcus.aureus + Staphylococcus.epidermidis + Staphylococcus.lugdunensis + Streptococcus.agalactiae + Age +
#          EV.1 + EV.2 + EV.3 + EV.4 + EV.5 + Diabetes, data = temp)

m = lm(duration ~ Pseudomonas.aeruginosa + 
             Anaerococcus.vaginalis +
         Age + EV.1 + EV.5 + Diabetes, data = temp)
summary(m)
anova(m)
Anova(m, type = 2)
1-(((1-(5235280/13481430))*(48-1))/(48-6-1)) # pa adjusted R2 = 0.299
1-(((1-(1561056/13481430))*(48-1))/(48-6-1)) # av adjusted R2 = -0.013
1-(((1-(3564712/13481430))*(48-1))/(48-6-1)) # age adjusted R2 = 0.157
1-(((1-(1961518/13481430))*(48-1))/(48-6-1)) # ev1 adjusted R2 = 0.020
1-(((1-(1170204/13481430))*(48-1))/(48-6-1)) # ev5 adjusted R2 = -0.047 
1-(((1-(1367134/13481430))*(48-1))/(48-6-1)) # ev5 adjusted R2 = -0.030

### plotting
# Code below is for making a single melted dataframe where each variable included in either one of two final models
#   will be a facet layer in a plot showing the relationship of each variable to healing duration.

ppa = temp[,c("duration", "Pseudomonas.aeruginosa")]
ppa = ppa[ppa$Pseudomonas.aeruginosa > 0,]
pav = temp[,c("duration", "Anaerococcus.vaginalis")]
pav = pav[pav$Anaerococcus.vaginalis > 0,]
age = temp[,c("duration", "Age")]
dia = temp[,c("duration", "Diabetes")]
pe1 = temp[,c("duration", "EV.1")]
pe4 = temp[,c("duration", "EV.4")]
alp = temp[,c("duration", "Hill")]
pe5 = temp[,c("duration", "EV.5")]

ppa$Species = names(ppa)[2]
pav$Species = names(pav)[2]
age$Species = names(age)[2]
dia$Species = names(dia)[2]
pe1$Species = names(pe1)[2]
pe4$Species = names(pe4)[2]
alp$Species = names(alp)[2]
pe5$Species = names(pe5)[2]
names(ppa)[2] = "ra"
names(pav)[2] = "ra"
names(age)[2] = "ra"
names(dia)[2] = "ra"
names(pe1)[2] = "ra"
names(pe4)[2] = "ra"
names(alp)[2] = "ra"
names(pe5)[2] = "ra"
both = rbind(ppa, pav, age, dia, pe1, pe4, alp, pe5)
both$Species = gsub("\\.", " ", both$Species)
both$Species = ordered(both$Species, levels = c("Hill", "Age", "Pseudomonas aeruginosa", "Anaerococcus vaginalis",
                                                "Diabetes", "EV 1", "EV 4", "EV 5"))

p = ggplot(both, aes(ra, duration)) + geom_point(size = 2.5, alpha = 0.5) + 
  facet_wrap(~Species, labeller = label_wrap_gen(multi_line = T, width = 15), nrow = 4, scales = "free_x") + 
  scale_y_continuous(expand = c(0,0), limits = c(-2000, 3500)) + 
  labs(x = "", y = "Healing Duration (days)") + 
  geom_smooth(aes(ra, duration),method = "lm") + coord_cartesian(ylim = c(0, 3000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(strip.text = element_text(face = "bold.italic", size = 8)) +
  theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10),
        axis.text.x = element_text(hjust = 1, vjust = 1, angle = 45))
tiff(file = "Figure_3.tiff", res = 600, width = 88, height = 180, units = 'mm', compression = 'lzw')
print(p2)
dev.off()



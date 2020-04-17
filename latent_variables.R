#####
#data development
#####
snp = read.delim("extracted_snp_merge.ped", header = F, row.names = 1, sep = " ", stringsAsFactors = FALSE)
snp = snp[,-c(2:5)]
n = read.csv("sig_snps_all.csv", stringsAsFactors = F)
names(snp) = c("SampleCode", n$SNP)
snp[snp == "00"] = NA

otu = read.csv("clean_species_table_p1.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu2 = read.csv("clean_species_table_p2.csv", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
otu = otu[complete.cases(otu),]
otu2 = otu2[complete.cases(otu2),]

p1 = read.csv("allele_calls_p1.csv")
p1 = names(p1)[3:ncol(p1)]
p1 = sub("GSA.", "", p1)

p2 = names(snp)[2:16]
p2 = p2[!p2 %in% p1]

group = data.frame(chartno = c(row.names(otu), row.names(otu2)), plate = c(rep("one", nrow(otu)), rep("two", nrow(otu2))))

ad_to_otu2 = names(otu)[which(!names(otu) %in% names(otu2))]
ad2 = otu2[,1:length(ad_to_otu2)]
names(ad2) = ad_to_otu2
ad2[ad2 > 0 ] = 0
otu2 = cbind(otu2, ad2)

ad_to_otu = names(otu2)[which(!names(otu2) %in% names(otu))]
ad = otu[,1:length(ad_to_otu)]
names(ad) = ad_to_otu
ad[ad > 0 ] = 0
otu = cbind(otu, ad)
rm(ad_to_otu, ad_to_otu2, ad2, ad)

otu2 = otu2[,c(match(names(otu), names(otu2)))]
all(names(otu) == names(otu2))

otu = rbind(otu,otu2)
rm(otu2)

met = read.csv("met1.csv")
met2 = read.csv("met2.csv")
met = rbind(met, met2)
rm(met2)

drop = row.names(otu)[which(!row.names(otu) %in% met$chartno)]
snp = snp[!snp$SampleCode %in% drop,]
otu = otu[!row.names(otu) %in% drop,]

row.names(otu) = met$SampleCode[match(row.names(otu), met$chartno)]
group$SampleCode = met$SampleCode[match(group$chartno, met$chartno)]

snp = snp[snp$SampleCode %in% row.names(otu),]
otu = otu[row.names(otu) %in% snp$SampleCode,]
met = met[met$SampleCode %in% row.names(otu),]

snp = snp[match(row.names(otu), snp$SampleCode),]
met = met[match(row.names(otu), met$SampleCode),]

all(snp$SampleCode == row.names(otu))
all(met$SampleCode == row.names(otu))

temp = otu
temp[temp > 0] = 1
rowSums(temp)
snp$hill1 = exp(diversity(otu))
snp$rich = rowSums(temp)
rm(temp)

snp = snp[complete.cases(snp),]

#####
#score genotypes as integers
#####
clone = snp
for(j in 2:16){
  states = sort(unique(clone[,j]))
  i = 1
  one = ifelse(clone[,j] == states[i], 0, 1) #dominant second allele
  two = ifelse(clone[,j] == states[(i+1)], 1, 0) #heterozygote superiority
  thr = ifelse(clone[,j] == states[(i+2)], 0, 1) #dominant first allele
  
  ad1 = clone[,j]
  ad1 = ifelse(ad1 == states[i], 0, ad1)
  ad1 = ifelse(ad1 == states[(i+1)], 1, ad1)
  ad1 = ifelse(ad1 == states[(i+2)], 2, ad1)
  ad1 = as.numeric(ad1) #forward direction additive
  
  ad2 = clone[,j]
  ad2 = ifelse(ad2 == states[i], 2, ad2)
  ad2 = ifelse(ad2 == states[(i+1)], 1, ad2)
  ad2 = ifelse(ad2 == states[(i+2)], 0, ad2)
  ad2 = as.numeric(ad2) #reverse direction additive
  
  hold = vector("list", 5)
  hold[[1]] = one
  hold[[2]] = two
  hold[[3]] = thr
  hold[[4]] = ad1
  hold[[5]] = ad2
  
  comp = data.frame(one,two,thr,ad1,ad2,clone$hill1)
  comp = cor(comp)
  comp = data.frame(comp[-nrow(comp),])
  keep = which(abs(comp$clone.hill1) == max(abs(comp$clone.hill1)))
  keep = max(comp$clone.hill1[keep])
  keep = which(comp$clone.hill1 == keep)
  if(comp$clone.hill1[keep] < 0) {
    hold[[keep]] = ifelse(hold[[keep]] == 1, 0, 1)
  }
  clone[,j] = hold[[keep]]
}


#####
#regression
#####
comp = cor(clone[,2:17])
comp = data.frame(comp[-nrow(comp),])
comp = comp[order(comp$hill1, decreasing = TRUE),]

summary(lm(hill1 ~ rs11125573, data = clone))

summary(lm(hill1 ~ rs201276730, data = clone))

summary(lm(hill1 ~ rs522554, data = clone))

summary(lm(hill1 ~ rs9307010, data = clone))

summary(lm(hill1 ~ rs3846499, data = clone))

summary(lm(hill1 ~ rs2864680, data = clone))

summary(lm(hill1 ~ clone$'GSA-rs11984782', data = clone))

summary(lm(hill1 ~ clone$'GSA-rs10815368', data = clone))

summary(lm(hill1 ~ clone$'GSA-rs10780221', data = clone))

summary(lm(hill1 ~ clone$'GSA-rs1340141', data = clone))

summary(lm(hill1 ~ clone$'GSA-rs4758411', data = clone))

summary(lm(hill1 ~ rs12307988, data = clone))

summary(lm(hill1 ~ rs1436708, data = clone))

summary(lm(hill1 ~ rs8031916, data = clone))

summary(lm(hill1 ~ rs7236481, data = clone))

#model for those that were sig in individual models
summary(lm(hill1 ~ rs522554  + rs9307010  + 
             rs3846499  + clone$'GSA-rs11984782'  + clone$'GSA-rs10815368'  + 
             clone$'GSA-rs1340141'  + clone$'GSA-rs4758411'  + rs12307988  + 
             rs1436708  + rs8031916  + rs7236481, data = clone))

#sequentially drop those that are not significant
#(-clone$"GSA-rs10815368")
summary(lm(hill1 ~ rs522554  + rs9307010  + 
             rs3846499  + clone$'GSA-rs11984782'  +  
             clone$'GSA-rs1340141'  + clone$'GSA-rs4758411'  + rs12307988  + 
             rs1436708  + rs8031916  + rs7236481, data = clone))

#sequentially drop those that are not significant
#(-rs9307010)
summary(lm(hill1 ~ rs522554  + 
             rs3846499  + clone$'GSA-rs11984782'  +  
             clone$'GSA-rs1340141'  + clone$'GSA-rs4758411'  + rs12307988  + 
             rs1436708  + rs8031916  + rs7236481, data = clone))

#sequentially drop those that are not significant
#(-clone$"GSA-rs4758411")
summary(lm(hill1 ~ rs522554  + 
             rs3846499  + clone$'GSA-rs11984782'  +  
             clone$'GSA-rs1340141'  + rs12307988  + 
             rs1436708  + rs8031916  + rs7236481, data = clone))



#####
#identify adjacent loci of interest
#####
lv  = read.csv("latent_var_requested_snps_cohort1.csv", stringsAsFactors = FALSE)
lv2  = read.csv("latent_var_requested_snps_cohort2.csv", stringsAsFactors = FALSE)
lv = rbind(lv, lv2)
rm(lv2)

for(i in unique(lv$hit)){
  print(paste(i, length(grep(i, lv$hit)), sep = "   "))
}
n = unique(lv$hit)
loc = vector("list", length = length(n))
for(i in 1:length(loc)){
  names(loc)[i] = n[i]
  loc[[i]] = lv[grep(n[i], lv$hit),]
  loc[[i]] = loc[[i]][order(loc[[i]]$P, decreasing = FALSE),]
}

#find elbows
look = c()
for(j in 1:length(loc)){
  temp = c()
  for(i in 1:20){
    h = abs(loc[[j]]$P[i] - loc[[j]]$P[(i+1)])
    temp = c(temp, h)
  }
  stp = which(temp == max(temp))
  look = c(look, stp)
}

#visually assess with plots
al = c()
for(i in 1:length(loc)){
  temp = loc[[i]][1:20,]
  al = rbind(al, temp)
}
al$place = 1:nrow(al)
library(ggplot)
ggplot(al, aes(place, P)) + geom_point() + facet_wrap(~hit, scales = "free")
look

look[c(3,5,8,9,10)] = c(9,5,6,12,14)

for(j in 1:length(loc)){
  loc[[j]] = loc[[j]][1:look[j],]
}

#####
#get the genotype calls for loci of interest
#####

geno = read.delim("extracted_latent_var_snps_cohort1.ped", sep = " ", header = FALSE, row.names = 1, stringsAsFactors = FALSE)
geno = geno[,-c(2:5)]
names(geno)[1] = "SampleCode"

geno2 = read.delim("extracted_latent_var_snps_cohort2.ped", sep = " ", header = FALSE, row.names = 1, stringsAsFactors = FALSE)
geno2 = geno2[,-c(2:5)]
names(geno2)[1] = "SampleCode"

all(geno$SampleCode == geno2$SampleCode)

geno = geno[,-1]
geno2 = geno2[,-1]
geno = cbind(geno, geno2)
rm(geno2)

geno[geno == "00"] = NA

nms = read.delim("latent_var_list_cohort1.txt", header = FALSE, stringsAsFactors = FALSE)
nms2 = read.delim("latent_var_list_cohort2.txt", header = FALSE, stringsAsFactors = FALSE)
nms = c(c(nms$V1, nms2$V1))
rm(nms2)

names(geno) = nms

geno = geno[row.names(geno) %in% clone$SampleCode,]
geno = geno[match(clone$SampleCode, row.names(geno)),]
all(row.names(geno) == clone$SampleCode)

genes = vector("list", length(loc))
names(genes) = names(loc)

for(i in 1:length(genes)){
  x = loc[[i]]$SNP
  x = match(x, names(geno))
  genes[[i]] = geno[,x]
}

#####
#integer code the loci of interest
#####
numgenes = genes
for(z in 1:length(numgenes)){
  for(j in 1:ncol(numgenes[[z]])){
    states = sort(unique(numgenes[[z]][,j]))
    i = 1
    one = ifelse(numgenes[[z]][,j] == states[i], 0, 1) #dominant second allele
    two = ifelse(numgenes[[z]][,j] == states[(i+1)], 1, 0) #heterozygote superiority
    if(length(states) == 3) {
      thr = ifelse(numgenes[[z]][,j] == states[(i+2)], 0, 1) #dominant first allele
      
      ad1 = numgenes[[z]][,j]
      ad1 = ifelse(ad1 == states[i], 0, ad1)
      ad1 = ifelse(ad1 == states[(i+1)], 1, ad1)
      ad1 = ifelse(ad1 == states[(i+2)], 2, ad1)
      ad1 = as.numeric(ad1) #forward direction additive
      
      ad2 = numgenes[[z]][,j]
      ad2 = ifelse(ad2 == states[i], 2, ad2)
      ad2 = ifelse(ad2 == states[(i+1)], 1, ad2)
      ad2 = ifelse(ad2 == states[(i+2)], 0, ad2)
      ad2 = as.numeric(ad2) #reverse direction additive
      
      hold = vector("list", 5)
      hold[[1]] = one
      hold[[2]] = two
      hold[[3]] = thr
      hold[[4]] = ad1
      hold[[5]] = ad2
      
      comp = data.frame(one,two,thr,ad1,ad2,clone$hill1)
      comp = cor(comp, use = "complete.obs")
      comp = data.frame(comp[-nrow(comp),])
      keep = which(abs(comp$clone.hill1) == max(abs(comp$clone.hill1)))
      keep = max(comp$clone.hill1[keep])
      keep = which(comp$clone.hill1 == keep)[1]
      print(paste(keep, z, j, "3"))
      if(comp$clone.hill1[keep] < 0) {
        hold[[keep]] = ifelse(hold[[keep]] == 1, 0, 1)
      }
      numgenes[[z]][,j] = hold[[keep]]
    }
    
    if(length(states) == 2) {
      hold = vector("list", 2)
      hold[[1]] = one
      hold[[2]] = two

      comp = data.frame(one,two,clone$hill1)
      comp = cor(comp, use = "complete.obs")
      comp = data.frame(comp[-nrow(comp),])
      keep = which(abs(comp$clone.hill1) == max(abs(comp$clone.hill1)))
      keep = max(comp$clone.hill1[keep])
      keep = which(comp$clone.hill1 == keep)[1]
      print(paste(keep, z, j, "2"))
      if(comp$clone.hill1[keep] < 0) {
        hold[[keep]] = ifelse(hold[[keep]] == 1, 0, 1)
      }
      numgenes[[z]][,j] = hold[[keep]]
    }
  }
}

#####
#look at correlations amonst adjacent loci of interest
#####

cors = vector("list", length(numgenes))
names(cors) = names(numgenes)
for(i in 1:length(cors)){
  cors[[i]] = cor(numgenes[[i]], use = "pairwise.complete.obs")
}

#####
#organize and sort the correlations
#####

for(j in 1:length(cors)){
  temp = cors[[j]]
  temp = melt(temp)
  temp$Var1 = as.character(temp$Var1)
  temp$Var2 = as.character(temp$Var2)
  temp = temp[temp$Var1 != temp$Var2,]
  
    #keeping only unique correlation entries
  for(i in 1:nrow(temp)){
    one = temp$Var1[i]
    two = temp$Var2[i]
    drop = intersect(which(temp$Var2 == one), which(temp$Var1 == two))
    if(length(drop) != 0){temp = temp[-drop,]}
    print(i)
  }
  
  temp = temp[order(temp$value, decreasing = TRUE),]
  temp$hit = rep(names(cors[j]))
  
  cors[[j]] = temp
}
  
#####
#subset the SNPs for those that may be indicators
#####

cors = numgenes

cors[[1]] = cors[[1]][,c("rs11125573", "GSA-rs10194425"), drop = FALSE]

cors[[2]] = cors[[2]][,c("rs201276730", "rs10469593", "rs10496839", "GSA-rs11894060"), drop = FALSE]
cors[[2]]$rs11894060_rs201276730  = (cors[[2]]$'GSA-rs11894060' + cors[[2]]$rs201276730)/2

cors[[3]] = cors[[3]][,c("GSA-rs10815368", "GSA-rs17432514", "GSA-rs59521969", "GSA-rs56296348"), drop = FALSE]
cors[[3]]$rs59521969_rs56296348  = (cors[[3]]$'GSA-rs59521969' + cors[[3]]$'GSA-rs56296348')/2

cors[[4]] = cors[[4]][,c("GSA-rs10780221", "GSA-rs7865979", "GSA-rs12519", "GSA-rs1331180"), drop = FALSE]
cors[[4]]$rs12519_rs1331180  = (cors[[4]]$'GSA-rs12519' + cors[[4]]$'GSA-rs1331180')/2

cors[[5]] = cors[[5]][,c("GSA-rs4758411", "rs2682111"), drop = FALSE]

cors[[6]] = cors[[6]][,c("rs1436708", "GSA-rs1925053", "rs9597069", "rs3125761"), drop = FALSE]
cors[[6]]$rs1925053_rs9597069 = (cors[[6]]$'GSA-rs1925053'   +  cors[[6]]$rs9597069)/2

cors[[7]] = cors[[7]][,c("rs8031916", "rs4775535", "rs938985", "rs73435666"), drop = FALSE]

cors[[8]] = cors[[8]][,c("rs7236481", "GSA-rs45590038", "rs45456199"), drop = FALSE]
cors[[8]]$rs45590038_rs45456199 = (cors[[8]]$'GSA-rs45590038'   +  cors[[8]]$rs45456199)/2

cors[[9]] = cors[[9]][,c("rs522554", "rs6770339"), drop = FALSE]

cors [[10]] = cors[[10]][,c("rs9307010", "rs9998733"), drop = FALSE]

cors [[11]] = cors[[11]][,c("rs3846499", "GSA-rs32483", "rs6881771", "GSA-rs250244"), drop = FALSE]
cors[[11]]$rs6881771_rs250244 = (cors[[11]]$rs6881771   +  cors[[11]]$'GSA-rs250244')/2

cors [[12]] = cors[[12]][,c("rs2864680"), drop = FALSE]

cors[[13]] = cors[[13]][,c("GSA-rs11984782", "rs11989865"), drop = FALSE]

cors [[14]]  = cors[[14]][,c("GSA-rs1340141", "rs7092247", "GSA-rs17136428", "GSA-rs74124948"), drop = FALSE]

cors[[15]] = cors[[15]][,c("rs12307988", "rs7137605", "exm2267480", "rs35989183", "rs11246976"), drop = FALSE]
cors[[15]]$rs7137605_exm2267480 = (cors[[15]]$rs7137605   *  cors[[15]]$'exm2267480')/2
cors[[15]]$rs11246976_rs35989183 = (cors[[15]]$rs35989183   *  cors[[15]]$'exm2267480')/2

for(i in 1:length(cors)){
  cors[[i]]$average = rowMeans(cors[[i]])
}


for(i in 1:length(cors)){
  cors[[i]]$hill1 = clone$hill1
}

#####
#correlation among potential indicators
#####

for(i in 1:length(cors)){
  cors[[i]] = cor(cors[[i]], use = "pairwise.complete.obs")
}

#####
#refined subset of indicators
#####

cors = numgenes

cors[[1]] = cors[[1]][,c("rs11125573", "GSA-rs10194425"), drop = FALSE]

cors[[2]]$rs11894060_rs201276730  = (cors[[2]]$'GSA-rs11894060' + cors[[2]]$rs201276730)/2
cors[[2]] = cors[[2]][,c("rs10469593", "rs10496839", "rs11894060_rs201276730"), drop = FALSE]

cors[[3]]$rs59521969_rs56296348  = (cors[[3]]$'GSA-rs59521969' + cors[[3]]$'GSA-rs56296348')/2
cors[[3]] = cors[[3]][,c("GSA-rs10815368", "GSA-rs17432514", "rs59521969_rs56296348"), drop = FALSE]

cors[[4]] = cors[[4]][,c("GSA-rs10780221", "GSA-rs7865979", "GSA-rs1331180"), drop = FALSE]

cors[[5]] = cors[[5]][,c("GSA-rs4758411", "rs2682111"), drop = FALSE]

cors[[6]] = cors[[6]][,c("rs1436708"), drop = FALSE]

cors[[7]] = cors[[7]][,c("rs8031916", "rs4775535", "rs938985"), drop = FALSE]

cors[[8]] = cors[[8]][,c("rs7236481", "rs45456199"), drop = FALSE]

cors[[9]] = cors[[9]][,c("rs522554", "rs6770339"), drop = FALSE]

cors [[10]] = cors[[10]][,c("rs9307010"), drop = FALSE]

cors [[11]] = cors[[11]][,c("rs3846499", "GSA-rs32483", "rs6881771"), drop = FALSE]

cors [[12]] = cors[[12]][,c("rs2864680"), drop = FALSE]

cors[[13]] = cors[[13]][,c("GSA-rs11984782", "rs11989865"), drop = FALSE]

cors [[14]]  = cors[[14]][,c("GSA-rs1340141", "GSA-rs17136428", "GSA-rs74124948"), drop = FALSE]

cors[[15]]$rs7137605_exm2267480 = (cors[[15]]$rs7137605   *  cors[[15]]$'exm2267480')/2
cors[[15]] = cors[[15]][,c("rs12307988", "rs7137605_exm2267480"), drop = FALSE]

#if you want to look at the correlations again
# for(i in 1:length(cors)){
#   cors[[i]]$average = rowMeans(cors[[i]])
# }
# 
# 
# for(i in 1:length(cors)){
#   cors[[i]]$hill1 = clone$hill1
# }
# 
# for(i in 1:length(cors)){
#   cors[[i]] = cor(cors[[i]], use = "pairwise.complete.obs")
# }

#####
#lavaan
#####

dat = cors[[1]]
for(i in 2:length(cors)){
  dat = cbind(dat, cors[[i]])
}
dat$hill1 = clone$hill1

names(dat) = sub("GSA-", "", names(dat))

library(lavaan)

mod = '
a =~ rs11125573 + rs10194425 + equal("a=~rs11125573")*rs10194425
hill1 ~ a
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
hill1 ~ b
'
test = cfa(mod, data = dat, std.lv=T) #reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
c =~ rs10815368 + rs17432514 + rs59521969_rs56296348
hill1 ~ c
'
test = cfa(mod, data = dat, std.lv=T)
summary(test, standardized=T, fit.measures=T)

mod = '
d =~ rs10780221 + rs7865979 + rs1331180
hill1 ~ d
'
test = cfa(mod, data = dat, std.lv=T) #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
hill1 ~ e
'
test = cfa(mod, data = dat, std.lv=T) # reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
f =~ rs1436708
hill1 ~ f
'
test = cfa(mod, data = dat, std.lv=T) #reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
g =~ rs8031916 + rs4775535 + rs938985
hill1 ~g
'
test = cfa(mod, data = dat, std.lv=T)  #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
h =~ rs7236481 + rs45456199 + equal("e=~rs7236481")*rs45456199
hill1 ~ h
'
test = cfa(mod, data = dat, std.lv=T)  #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
i =~ rs522554 + rs6770339 + equal("e=~rs522554")*rs6770339
hill1 ~ i
'
test = cfa(mod, data = dat, std.lv=T)   #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
j =~ rs9307010
hill1 ~ j
'
test = cfa(mod, data = dat, std.lv=T)  #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
k =~ rs3846499 + rs32483 + rs6881771
hill1 ~ k
'
test = cfa(mod, data = dat, std.lv=T)  #reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
l =~ rs2864680
hill1 ~ l
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)

mod = '
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
hill1 ~ m
'
test = cfa(mod, data = dat, std.lv=T)   #reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

mod = '
n =~ rs1340141 + rs17136428 + rs74124948
hill1 ~ n
'
test = cfa(mod, data = dat, std.lv=T) #some variances are negative #reg is sig
summary(test, standardized=T, fit.measures=T)

mod = '
o =~ rs12307988 + rs7137605_exm2267480
hill1 ~ o
'
test = cfa(mod, data = dat, std.lv=T) #reg is sig
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

#reverting to original SNPs for those that didn't show improvement for latent as compared to focal SNP lm
#a,c,g,l,n
#a,c,l are actually removed becuase they were n.s. in the single models


mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
d =~ rs10780221 + rs7865979 + rs1331180
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
g =~ rs8031916
h =~ rs7236481 + rs45456199 + equal("e=~rs7236481")*rs45456199
i =~ rs522554 + rs6770339 + equal("e=~rs522554")*rs6770339
j =~ rs9307010
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + d + e + f + g + h + i + j + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)


#h gets dropped because it breaks the covariance matrix
mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
d =~ rs10780221 + rs7865979 + rs1331180
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
g =~ rs8031916
i =~ rs522554 + rs6770339 + equal("e=~rs522554")*rs6770339
j =~ rs9307010
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + d + e + f + g + i + j + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

#drop g

mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
d =~ rs10780221 + rs7865979 + rs1331180
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
i =~ rs522554 + rs6770339 + equal("e=~rs522554")*rs6770339
j =~ rs9307010
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + d + e + f  + i + j + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

#drop i

mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
d =~ rs10780221 + rs7865979 + rs1331180
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
j =~ rs9307010
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + d + e + f + j + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

#drop j

mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
d =~ rs10780221 + rs7865979 + rs1331180
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + d + e + f + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')

#drop d

mod = '
b =~ rs10469593 + rs10496839 + rs11894060_rs201276730
e =~ rs4758411 + rs2682111 + equal("e=~rs4758411")*rs2682111
f =~ rs1436708
k =~ rs3846499 + rs32483 + rs6881771
m =~ rs11984782 + rs11989865 + equal("e=~rs11984782")*rs11989865
o =~ rs12307988 + rs7137605_exm2267480

hill1 ~  b + e + f + k + m + o
'
test = cfa(mod, data = dat, std.lv=T) 
summary(test, standardized=T, fit.measures=T)
inspect(test, 'r2')









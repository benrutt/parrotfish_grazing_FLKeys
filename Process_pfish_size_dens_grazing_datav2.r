##RVC Parrotfish processing for grazing metrics
## B. Ruttenberg May 2016

response <- 'num'
#initial header
header.stnr <- c('year','spcode', 'strat', 'psu','station_nr', response)
header.stnr.1 <- header.stnr[1:(length(header.stnr)-1)]
header.psu <- c(header.stnr[1:(length(header.stnr)-2)],paste("avg.", header.stnr[length(header.stnr)], sep = ""))
header.psu.1 <- c(header.psu[1:(length(header.psu)-1)])
header.strat <- c(header.psu[1:(length(header.psu)-2)], header.psu[length(header.psu)])
header.strat.1 <- c(header.strat[1:(length(header.strat)-1)])
header.domain <- c(header.strat[1:(length(header.strat)-2)], header.strat[length(header.strat)])
header.domain.1 <- c(header.domain[1:(length(header.domain)-1)])
######next header needs to be checked manually each time 'header.stnr' changes
header.ntot <- c('prot', 'strat')

####Basic FL data####

#set wd
setwd('C:/Users/bruttenb/Google Drive/Project Files/Parrotfish/Analyses/Grazing metrics')

library(reshape);library(reshape2)

#standard pfsih dataset
pfish <- read.table('pfish.csv', header = T, sep = ',')

pfish$len <- ifelse(pfish$len < 0, -9, pfish$len)
pfish$len <- round(pfish$len, 0)

denssize <- aggregate(cbind(pfish[response]), pfish[c(header.stnr.1, 'len')], FUN = sum, na.rm = T)

#reshape to wide dataset
denssize.w <- reshape(denssize, timevar = 'len', idvar = header.stnr.1, direction = 'wide')

#split denssize to categories and lengths, and set len = NA to 0
denssize.w.cat <- denssize.w[,1:5]
denssize.w.lens <- denssize.w[,6:ncol(denssize.w)]
denssize.w.lens[is.na(denssize.w.lens)] <- 0
#reform the denssize.wide dataset
denssize.wide <- cbind(denssize.w.cat,denssize.w.lens)


#removes 'num.' from variable names
names(denssize.wide) <- gsub('num.','', names(denssize.wide))


#reshapes to 'long' form again
a <- proc.time()
denssize.l <- reshape(denssize.wide, varying = names(denssize.wide)[6:ncol(denssize.wide)], v.names = 'num', timevar = 'len', times = names(denssize.wide)[6:ncol(denssize.wide)], direction = 'long')
b <- proc.time(); b-a

row.names(denssize.l) <- 1:nrow(denssize.l)
denssize.l$id <- NULL

write.table(denssize.l, 'denssize.l.csv', sep = ',', row.names = F, col.names = T)

denssize.l <- read.table('denssize.l.csv', header = T, sep = ',')


#mean of densitySize over PSU
denssize.psu <- aggregate(cbind(denssize.l[response]), denssize.l[c(header.psu.1, 'len')], FUN = mean, na.rm = T)
denssize.psu$len <- as.numeric(denssize.psu$len)

#sum of all fish of a given species by psu, to create weighting factor for size 
dens.psu <- with(denssize.psu[denssize.psu$len > 0,],
                 aggregate(cbind(denssize.psu[response]), denssize.psu[header.psu.1], FUN = sum, na.rm = T))
names(dens.psu)[names(dens.psu) =='num'] <- 'sumdens'

#merge dens.psu with new sum of num per species/psu 
denssize.sum.psu <- merge(denssize.psu, dens.psu, by = header.psu.1, all.x = T, all.y = F)

#calc weighting factor for size
denssize.sum.psu$wf <- denssize.sum.psu$num/denssize.sum.psu$sumdens
#calc weighted size
denssize.sum.psu$wsize <- with(denssize.sum.psu, ifelse(len <= 0, NA, len*wf)) 
#calc weighted size per species and PSU
meanlen.psu <- aggregate(denssize.sum.psu['wsize'], denssize.sum.psu[header.psu.1], FUN = sum, na.rm = T)
names(meanlen.psu)[names(meanlen.psu) =='wsize'] <- 'meanlen'
#merge denssize.sum.psu with meanlen.psu
denssize.sum.psu.1 <- merge(denssize.sum.psu, meanlen.psu, by = header.psu.1, all.x = T, all.y = F)  
#replaces len of <= 0 with mean len (ROUNDED) for the Species/PSU
denssize.sum.psu.1$len <- with(denssize.sum.psu.1, ifelse(len <= 0, round(meanlen,0), len)) 

denssize.sum.psu.2 <- aggregate(denssize.sum.psu.1['num'], denssize.sum.psu.1[c(header.psu.1, 'len')], FUN = sum, na.rm = T)

#mean of denssize over strat
denssize.strat <- aggregate(cbind(denssize.sum.psu.2[response]), denssize.sum.psu.2[c(header.strat.1, 'len')], FUN = mean, na.rm = T)

write.table(denssize.strat, 'density_size_strat2.csv', row.names = F, col.names = T, sep = ',')

##Combine data

#reads in data frame created above with density and size by strat and year
denssize.strat <- read.table('density_size_strat2.csv', header = T, sep = ',')

#reads in species metrics data
spmetrics <- read.table('grazing_metrics_single_table.csv', header = T, sep = ',')
#species metrics table assuming all species in the same genera are the same
spmetrics <- read.table('grazing_metrics_single_table_by_genera_only.csv', header = T, sep = ',')

#combines species metrics and species by strat data
grazing <- merge(denssize.strat, spmetrics, by = 'spcode', all.x = T, all.y = F)

cyl <- pi*(7.5^2)

#Calc grazing metrics
grazing$indiv_w <- with(grazing, lw_a*(len^lw_b)) #units: g/cyl
grazing$total_w <- with(grazing, indiv_w*num)/cyl #units: g/m2
grazing$algcons <- with(grazing, (0.0342*(indiv_w^0.816))*num*365)/cyl #units: g/m2; allometric relationship from Bruggeman (note: * 365 converts from day to year)
grazing$allmacro <- with(grazing, algcons*propmacro) #units: g/m2/yr 
grazing$brownmacro <- with(grazing, algcons*propbrown) #units: g/m2
grazing$scarsize <- with(grazing, bitescar_m*(len^2))#units: cm2/bite; cut: *num*720)/cyl #720 min/per 12 hr day per m2
grazing$biterate <- with(grazing, (biterate_m*len+biterate_b)*num*600) #units: bites/day/line; 600 min/per 10 hr day
grazing$totalscar <- with(grazing, (biterate*scarsize))/cyl #units: cm2/day/m2
grazing$scarpercent <- with(grazing, totalscar*365/100) #units: percent of reef grazed per year
grazing$bioeros <- with(grazing, ((bioeros_v*(len^3))*biterate*1.7*365)/(cyl*1000)) #using 1.7 g/cm3 to convert from volume grazed to g of bioerosion (Mallela, J. & Perry, C.T. Coral Reefs (2007) 26: 129. doi:10.1007/s00338-006-0169-7), converts to kg/yr/m2




grazing1 <- grazing[grazing$spcode %in% unique(spmetrics$spcode) & grazing$len>10,]
#grazing1$propmacroweight <- grazing1$allmacro/grazing1$algcons

###GRAZING METRICS BY LINE for FL####

#by species sum across strat; for algcons, allmacro, brownmacro, total_w; exclude scar bc of NA values
byyear.sp1 <- with(grazing1, aggregate(cbind(algcons, allmacro, brownmacro, total_w)~spcode+year+strat, FUN = 'sum', na.rm = T))
#by species sum across strat; for totalscar only
byyear.sp2 <- with(grazing1, aggregate(cbind(totalscar, bioeros)~spcode+year+strat, FUN = 'sum', na.rm = T))
#merge byyear.sp datasets
byyear.sp3 <- merge(byyear.sp1, byyear.sp2, by = c('spcode', 'year', 'strat'), all.x = T, all.y = F)
#replaces totalscar with NA values
byyear.sp3$totalscar <- ifelse(is.na(byyear.sp3$totalscar), 0, byyear.sp3$totalscar)
byyear.sp3$bioeros <- ifelse(is.na(byyear.sp3$bioeros), 0, byyear.sp3$bioeros)

byyear.sp <- merge(byyear.sp3, spmetrics[c('spcode', 'speciesabbrev')], by = 'spcode', all.x = T, all.y = F)

strata <- as.data.frame(sort(unique(grazing1$strat)))
names(strata) <- 'strat'
strata$stratname <- c('Forereef Deep', 'Forereef Mid', 'Forereef Shallow', 'Hi-Rel Reef','Inshore Patch','Mid Patch', 'Offshore Patch')
strata$plotorder <- c(7,6,5,4,1,2,3)
strata <- strata[order(strata['plotorder']),]


######Length Freq plots for FL#####

#creates len freq dataset for 2012 only and num > 0
grazing.len.freq <- grazing1[grazing1$year == '2012' & grazing1$num >0, 1:5]
#divides by smallest value of n for each species
for(s in unique(grazing.len.freq$spcode)){
grazing.len.freq$adjnum[grazing.len.freq$spcode == s] <- grazing.len.freq$num[grazing.len.freq$spcode == s]/(min(grazing.len.freq$num[grazing.len.freq$num > 0 & grazing.len.freq$spcode == s]))
}

png(paste('FLKeys_size_freq_per_species.png', sep = ''), width = 5, height = 5, units = 'in', res = 300, bg = "white")

par(mfrow = c(3,3), mar = c(2.5,2,2,1.5), oma = c(0,2,2,0))

big3 <- c('SCA COEL', 'SCA COER', 'SCA GUAC')

for(s in c(unique(as.character(grazing.len.freq$spcode[!grazing.len.freq$spcode %in% big3])), big3)){
tempmat <- grazing.len.freq[grazing.len.freq == s,]
temphist <- rep(tempmat$len, round(tempmat$adjnum))
hist(temphist, main = spmetrics$speciesabbrev[spmetrics$spcode == s], cex.main = 0.95, breaks = seq(0, 75, by = 5),
     freq = F, xlab = 'Length', ylab = 'Rel Freq.', cex.lab = 1, cex.axis = 0.75, col = 'dark red')
mtext(paste('\u03bc','=', round(mean(temphist),1), sep=''), side = 3, adj = .95, line = -1, cex = 0.7)
}
dev.off()    



####Multivariate community plots for FL####

library(vegan)
library(MASS)


multisp <- reshape(lastyrs[,c('spcode', 'stratname', 'total_w')], timevar = 'spcode', idvar = c('stratname'), direction = 'wide')
names(multisp) <- gsub('total_w.','', names(multisp))

y <- 2011

png(paste('FL_NMDS_by_hab_w_fish.png', sep = ''), width = 6, height = 6, units = 'in', res = 300, bg = "white")

par(mfrow = c(1,1), mar = c(2.5,3,2,2), oma = c(3,3,3,0))

for(y in 2003:2012){
ordmat <- multisp[multisp$year == y,-1]
row.names(ordmat) <- ordmat$stratname
ordmat <- ordmat[,-1]

ord <- metaMDS(ordmat)

#creates vectors from ord results
xx <- as.vector(ord$points[,1])
yy <- as.vector(ord$points[,2])
sx <- as.vector(ord$species[,1])
sy <- as.vector(ord$species[,2])

#plots NMDS plots, 
plot(xx,yy, xlim = c(min(xx,sx)*1.3, max(xx,sx)*1.3), ylim = c(min(yy,sy)*1.1, max(yy,sy)*1.1), pch = 21, col = 'dark orange', bg = 'orange', cex = 1.4, xlab = 'MDSA1', ylab = 'MDSA2')
text(xx,yy,labels = strata$stratname, pos = 3, offset = 0.5, cex = 1.1)
text(sx, sy, labels = sort(unique(spmetrics$speciesabbrev)), col = 'blue') #omit this just to show 
mtext(paste('stress=', round(ord[[22]],4),sep = ''),side = 3, line = 0, adj = 0.95, cex = 0.9)
dev.off()

plot(ord, type = 't', main = y)
#points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
text(ord, display = "spec", cex=0.7, col="blue")
mtext(paste('stress=', round(ord[[22]],4),sep = ''),side = 3, line = 0, adj = 0.95, cex = 0.9)
}

dev.off()

##Principal components
par(mfrow = c(4,3), mar = c(2.5,3,2,2), oma = c(0,2,2,0))

for(y in 2001:2012){
  ordmat <- multisp[multisp$year == y,-1]
  row.names(ordmat) <- ordmat$strat
 # ordmat <- ordmat[,-1]
}  
  fl.pca <- princomp(ordmat[,1:6])
plot(fl.pca[[6]][,1], fl.pca[[6]][,2])



par(mfrow = c(2,2), mar = c(2.5,3,2,2), oma = c(0,2,2,0))

for(y in 2009:2012){
  ordmat <- multisp[multisp$year == y,-1]
  row.names(ordmat) <- ordmat$strat
  ordmat <- ordmat[,-1]
  
  d <- dist(ordmat)
  dmat <- isoMDS(d, k = 2)
  
  plot(dmat$points[,1], dmat$points[,2], main = y)
  #points(ord, display = "sites", cex = 0.8, pch=21, col="red", bg="yellow")
  text(dmat$points[,1], dmat$points[,2], labels = row.names(dmat), cex=0.7, col="blue")
  mtext(paste('stress=', round(ord[[22]],4),sep = ''),side = 3, line = 0, adj = 0.95, cex = 0.9)
}


#####Stacked Biomass plots for FL####
resp <- as.data.frame(cbind(c('algcons', 'totalscar', 'allmacro', 'brownmacro', 'total_w', 'scarpercent', 'bioeros'),c('g of C' ~ m^-2 ~ year^-1 , cm^2 ~ day^-1 ~ m^-2,'g of C' ~ m^-2 ~ year^-1, 'g of C' ~ m^-2 ~ year^-1, 'g'~m^-2,'% area grazed'~year^-1, 'kg'~m^-2~year^-1)))
names(resp) <- c('resp', 'units')
resp$places <- c(2.5, 2.5, 2.5, 2.5, 2, 2.25, 2)

#old response unit dataframe
#resp <- as.data.frame(cbind(c('algcons', 'totalscar', 'allmacro', 'brownmacro', 'total_w', 'scarpercent', 'bioeros'),c('g of C/day/m', 'cm2/day/m2', 'g of C/day/m2', 'g of C/day/m2', 'g/m2', '% area grazed/year', 'kg/m2/year')))
#names(resp) <- c('resp', 'units')
#resp$places <- c(2.5, 2.5, 2.5, 2.5, 2, 2.25, 2)

byyear.sp.barplot <- merge(byyear.sp, strata, by = 'strat', all.x = T, all.y = F)

write.table(byyear.sp.barplot, 'byyear_sp_barplot.csv', sep = ',', row.names = T, col.names = T)

#plots group of plots by set of data (by hab or by year)
set <- 'by_hab' #'by_year' or 'by_hab'
  
#for HRSG over 30 years
if(set == 'by_year'){
lastyrs <- byyear.sp.barplot[byyear.sp.barplot$strat == 'HRRF',]
lastyrs <- lastyrs[order(lastyrs$year, lastyrs$spcode),]
}

#for mean of last years
if(set == 'by_hab'){
yearuse <- 2003
lastyrs <- with(byyear.sp.barplot[byyear.sp.barplot$year >= yearuse,], aggregate(cbind(algcons, allmacro, total_w, totalscar, bioeros)~stratname+spcode+speciesabbrev, FUN = 'mean'))
#SE of last years
lastyrs.hab <- with(byyear.sp.barplot[byyear.sp.barplot$year >= yearuse,], aggregate(cbind(algcons, allmacro, total_w, totalscar, bioeros)~stratname+year, FUN = 'sum'))
lastyrsSE <- with(lastyrs.hab, aggregate(cbind(algcons, allmacro, total_w, totalscar, bioeros)~stratname, FUN = 'sd'))
lastyrsSE[,-1] <- lastyrsSE[,-1]/sqrt(length(unique(byyear.sp.barplot$year[byyear.sp.barplot$year >= yearuse])))
lastyrsSE <- lastyrsSE[c(5,6,7,4,3,2,1),]
lastyrsSE$scarpercent <- lastyrsSE$totalscar*(365/100) #converts to percent per year
}
lastyrs$scarpercent <- lastyrs$totalscar*(365/100) #converts to percent per year

for(response in c('total_w', 'allmacro', 'scarpercent', 'bioeros')){
response <- 'scarpercent'

#for 30 year in HRRF
if(set == 'by_year'){
  plotmat <- reshape(lastyrs[, c('speciesabbrev', 'year', response)],timevar = 'year', idvar = 'speciesabbrev', direction = 'wide')
}
if(set == 'by_hab'){
plotmat <- reshape(lastyrs[, c('speciesabbrev', 'stratname', response)],timevar = 'stratname', idvar = 'speciesabbrev', direction = 'wide')#mean of last 3 years
}
  names(plotmat) <- gsub(paste(as.character(response), '.', sep =''),'', names(plotmat))
row.names(plotmat) <- plotmat$speciesabbrev
plotmat$speciesabbrev <- NULL
plotmat <- as.matrix(plotmat)
#reorders habitats
if(set == 'by_hab'){
plotmat <- plotmat[,c(5,6,7,4,3,2,1)] 
}

#stores sums across habitats using species or generic level parameters, and writes the output. Note: need to change object and file names each time through
#macro.genera.params <- colSums(plotmat)
#write.csv(as.data.frame(bioeros.species.params/bioeros.genera.params), "bioerosspovergenera.csv")


#quick code to calc % contribution of each species for each response
View(round(prop.table(plotmat, 2)*100,2))
byhab <- colSums(round(prop.table(plotmat, 2)*100,2))
#diversity calcs
colSums(prop.table(plotmat, 2)*log(prop.table(plotmat, 2)))

#stats on differences between habitats in biomass and 3 responses
habanova <- aov(bioeros~stratname, data = lastyrs.hab)
summary(habanova)
tukey.habaov <- as.data.frame(TukeyHSD(habanova)$stratname)
names(tukey.habaov)[4] <- 'padj'
tukey.habaov[tukey.habaov$padj<0.05,]


#colSums(prop.table(plotmat, 2)[c('Sp. aurofrenatum', 'Sp. rubripinne'),]) #for macro SPAU and SPRU only

png(paste('FL', set, response, 'stacked_w_leg.png',sep ='_'), width = 10, height = 7, bg = 'white', res = 300, units = 'in')

par(mfrow = c(1,1), mar = c(12,8,5,11), oma = c(1,1,1,1), xpd = T)
mycolors <- c(as.character(spmetrics$color[spmetrics$speciesabbrev %in% sort(unique(byyear.sp.barplot$speciesabbrev))]))

#x <- c('dark red', 'red', 'orange',' yellow', 'light green', 'dark green', 'light blue', 'blue', 'dark blue', 'purple')
barcenters <- barplot(plotmat, col = mycolors, las = 2, cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5, main = response, xpd = F, ylim = c(0, 1.1*(max(colSums(plotmat))+max(lastyrsSE[,response]))))#, ylab = resp$units[resp$resp == response]) #legend = rownames(plotmat), 
arrows(barcenters, colSums(plotmat), barcenters, colSums(plotmat)+lastyrsSE[,response], angle = 90, length = 0.1)

#legend(5,15, rownames(plotmat), fill = mycolors, cex = 1.2)
legend("topright", inset = c(-0.35, 0), rownames(plotmat), fill = mycolors, cex = 1.2)
title(ylab = resp$units[resp$resp == response], cex.lab = 2, line = resp$places[resp$resp == response]*1.75)

dev.off()

}


####Panel Scatterplots of metrics##### 

#creates function to insert r2 and p value
panel.cor.p <- function(x, y, digits = 2, cex.cor)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient, r2
  r2 <- (cor(x, y))^2
  txt1 <- format(c(r2, 0.123456789), digits = digits)[1]
  txt2 <- bquote(r^{2}~'='~.(txt1))
  if(r2<0.01) txt2 <- bquote(r^{2}~'< 0.01')
  
  text(0.5, 0.6, txt2)
  
  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt3 <- format(c(p, 0.123456789), digits = digits)[1]
  txt3 <- paste("p = ", txt3, sep = "")
  if(p<0.001) txt3 <- "p < 0.001"
  text(0.5, 0.4, txt3)
}


png('FL_by_hab_scatter_10yravg_no_HRSG.png', width = 6, height = 6, bg = 'white', res = 300, units = 'in')
sum.metrics <- with(lastyrs, aggregate(cbind(scarpercent, allmacro, bioeros, total_w)~stratname, FUN = sum))
sum.metrics <- sum.metrics[!sum.metrics$stratname == 'Hi-Rel Reef',]
pairs(~scarpercent + allmacro + bioeros + total_w, data = sum.metrics, lower.panel = panel.cor.p)
dev.off()

png('FL_by_hab_scatter_10yrs_pts_no_HRSG.png', width = 6, height = 6, bg = 'white', res = 300, units = 'in')
sum.metrics.yr <- with(byyear.sp.barplot[byyear.sp.barplot$year>=2003,], aggregate(cbind(totalscar, allmacro, bioeros, total_w)~stratname+year, FUN = sum))
sum.metrics.yr$scarpercent <- sum.metrics.yr$totalscar*(365/100) 
sum.metrics.yr <- sum.metrics.yr[!sum.metrics.yr$stratname == 'Hi-Rel Reef',]
pairs(~scarpercent + allmacro + bioeros + total_w, data = sum.metrics.yr, lower.panel = panel.cor.p)
dev.off()

#####Multivariate plots using CAP analysis for FL####
library(plyr)
library(vegan)

biomass <- byyear.sp.barplot[byyear.sp.barplot$year > 2002,c('year', 'spcode', 'strat','total_w')]
biomass$stratyr <- with(biomass, paste(strat, year, sep = ' '))

biomass.wide <- reshape(biomass[,c('spcode', 'stratyr', 'total_w')],timevar = 'spcode', idvar = 'stratyr', direction = 'wide')
names(biomass.wide) <- gsub('total_w.', '', names(biomass.wide))


for(j in 1:nrow(biomass.wide)){
biomass.wide$strat[j] <- strsplit(biomass.wide$stratyr, ' ')[[j]][1]
biomass.wide$year[j] <- strsplit(biomass.wide$stratyr, ' ')[[j]][2]
}
biomass.wide$stratyr <- NULL
biomass.wide <- biomass.wide[order(biomass.wide$year, biomass.wide$strat),]

biomass.siteid <- biomass.wide[,c('strat', 'year')]
biomass.sp <- log(biomass.wide[,1:10]+1)

#CAP analysis (Constrained analysis of principal coordinates)
cap <-capscale(biomass.sp~strat, biomass.siteid, dist = 'bray')
#plotting specifications for pch and color
pchs <- c(21,22,23:25,21,22)
#capcolors <- 1:7
capcolors <- rainbow(7)

#creates data frame to plot sites
cap.sites <- as.data.frame(summary(cap)[5][[1]])
cap.sites$strat <- row.names(cap.sites)
cap.sites$strat <- gsub('strat','',cap.sites$strat)
cap.sites <- merge(cap.sites, strata, by = 'strat', all.x = T, all.y = F)
cap.sites$pchs <- pchs
cap.sites$capcolors <- capcolors

#creates data frame to plot arrows
species.arrows <- as.data.frame(summary(cap)[1][[1]])
species.arrows$spcode <- rownames(species.arrows)
species.arrows <- merge(species.arrows, spmetrics[,c('spcode', 'color')], by = 'spcode', all.x = T, all.y = F)
originpts <- rep(0,nrow(species.arrows))
#pchs <- as.data.frame(unique(biomass.wide$strat)); names(pchs) <- 'strat'
#pchs$sym <- 1:7
#summary(cap)

png('FL_RVC_CAP_2panel.png', width = 6, height = 10, bg = 'white', res = 300, units = 'in')

par(mfrow = c(2,1), mar = c(2,2,2,2), oma = c(3,3,1,1))

plot(cap, display = 'sites', type = 'n', xlim = c(-1.25,1.75))
with(cap.sites, points(CAP1, CAP2, pch = pchs, col = 'black', bg = capcolors, cex = 1.5))
#points(cap, display = 'cn', pch = pchs, col = 'black', bg = capcolors, cex = 1.5)
ordiellipse(cap, biomass.siteid$strat, kind = 'se', conf = 0.95, col = capcolors)
#legend('topleft', legend = levels(as.factor(biomass.siteid$strat)), pch = pchs, col = 'black', pt.bg = capcolors, bty = 'n')
legend('topleft', legend = as.factor(cap.sites$stratname[order(cap.sites$plotorder)]), pch = cap.sites$pchs[order(cap.sites$plotorder)], col = 'black', pt.bg = cap.sites$capcolors[order(cap.sites$plotorder)], bty = 'n')

plot(cap, display = 'sites', type = 'n', xlim = c(-1.25,1.75))
arrows(x0 = originpts, x1 = species.arrows$CAP1, y0 = originpts, y1 = species.arrows$CAP2, length = 0.1, lwd = 2, col = as.character(species.arrows$color))
legend("bottomright",  levels(spmetrics$speciesabbrev), fill = as.character(spmetrics$color), cex = 0.8, bty = 'n')
#text(cap, display = 'species', cex = 0.7)
mtext('CAP 1', 1, line = 1, outer = T, cex = 1.25)
mtext('CAP 2', 2, line = 1, outer = T, cex = 1.25)

dev.off()

####PerMANOVA test FL RVC#####

adonis(biomass.sp~strat, biomass.siteid, method = 'bray', sqrt.dist = T)


####NMDS Plot for FL####

ordmat <- biomass.wide[order(biomass.wide$year, biomass.wide$strat),]
ordmat.siteid <- ordmat$strat
row.names(ordmat) <- paste(ordmat$strat,ordmat$year, sep ='')
ordmat <- subset(ordmat, select = -c(strat, year))

ord <- metaMDS(ordmat)

plot(ord)
plot(ord, display = 'sites', main  = 'NMDS', type = 'n')#, xlim = c(-1.25,1.75))
points(ord, display = 'sites', pch = pchs, col = 'black', bg = capcolors, cex = 1.5)
ordiellipse(ord, ordmat.siteid, kind = 'se', conf = 0.95, col = capcolors)
legend('bottomleft', legend = levels(as.factor(ordmat.siteid)), pch = pchs, col = 'black', pt.bg = capcolors, bty = 'n')


#creates vectors from ord results
xx <- as.vector(ord$points[,1])
yy <- as.vector(ord$points[,2])
sx <- as.vector(ord$species[,1])
sy <- as.vector(ord$species[,2])


######Run Rmisc to summarize data below####

sp.bystrat <- with(byyear.sp[byyear.sp$year >= 2010,], aggregate(cbind(algcons, allmacro, brownmacro, totalscar, total_w)~spcode+strat, FUN = 'mean', na.rm = T))

sp.bystrat.sd <- with(byyear.sp[byyear.sp$year >= 2010,], aggregate(cbind(algcons, allmacro, brownmacro, totalscar, total_w)~spcode+strat, FUN = 'sd', na.rm = T))

#inserts sd values for sp.bystrat.sd and converts to SE
for(i in 1:5){
sp.bystrat[,i+7] <- sp.bystrat.sd[,i+2]/sqrt(3)  
  }
names(sp.bystrat)[8:12] <- paste(names(sp.bystrat[3:7]), '_SE', sep = '')

#plots responses by species.
for(j in resp$resp){
  png(paste('sp_w_', j, '_per_strat.png', sep = ''), width = 5, height = 7, units = 'in', res = 300, bg = "white")
  
  par(mfrow = c(4,2), mar = c(2.5,3,2,2), oma = c(0,2,2,0))
  
  for(i in unique(sp.bystrat$strat)) {
    barplot(sp.bystrat[sp.bystrat$strat == i,j], names.arg = sp.bystrat[sp.bystrat$strat == i,c('spcode')] , col = 'blue', ylim = c(0,max(sp.bystrat[j])*1.05), xlab = "", las = 2, cex.names = 0.8)
mtext(i, side = 3, line = -1.5)
  }
  title(main=toupper(j),outer=T)
  mtext(resp$units[resp$resp == j], side = 2, outer = T)
  
  dev.off()
}



#all species sum across strat
byyear <- with(grazing1, aggregate(cbind(algcons, totalscar, allmacro, brownmacro, total_w)~year+strat, FUN = 'sum', na.rm = T))

#byyear <- byyear[byyear$year > 1995,]
byyear$scarpercent <- byyear$totalscar/100
#removes outlier biomass values
#byyear <- byyear[byyear$total_w < 650,]

resp <- as.data.frame(cbind(c('algcons', 'totalscar', 'allmacro', 'brownmacro', 'total_w', 'scarpercent', 'bioeros'),c('g of C/day/m2', 'cm2/day/m2', 'g of C/day/m2', 'g of C/day/m2', 'g/m2', '% area grazed/year', 'kg/m2/year')))
names(resp) <- c('resp', 'units')


#Time series of each response by strat
for(j in resp$resp){
png(paste('year_w_', j, '_per_strat.png', sep = ''), width = 5, height = 7, units = 'in', res = 300, bg = "white")

par(mfrow = c(4,2), mar = c(2.5,2.5,2,2), oma = c(0,2,2,0))

for(i in unique(byyear$strat)) {
  plot(byyear[byyear$strat == i,c('year')], byyear[byyear$strat == i,j], pch = 21, col = 'black', bg = 'blue', main = i, xlim = c(min(byyear$year),max(byyear$year)),ylim = c(0,max(byyear[j])*1.05), xlab = "")
}
title(main=toupper(j),outer=T)
mtext(resp$units[resp$resp == j], side = 2, outer = T)

dev.off()
}

####Total parrotfish biomass by response by habitat for FL####

for(j in resp$resp[resp$resp !='total_w']){
png(paste('biomass_w_', j, '_per_strat.png', sep = ''), width = 5, height = 7, units = 'in', res = 300, bg = "white")

par(mfrow = c(4,2), mar = c(2.5,2.5,2,2), oma = c(2,2,2,0))

for(i in unique(byyear$strat)) {
plot(byyear[byyear$strat == i,c('total_w')], byyear[byyear$strat == i,j], pch = 21, col = 'black', bg = 'red', main = strata$stratname[strata$strat == i], xlim = c(0,max(byyear$total_w)),ylim = c(0,max(byyear[j])*1.05), xlab = "")
mtext(paste('r2=',round(summary(lm(byyear[byyear$strat == i,j]~byyear[byyear$strat == i,c('total_w')]))[[8]],3),sep = ''), cex = .7, line = -1.25, adj = 0.05)

}
mtext('Total parrotfish biomass (g/m2)', side = 1, outer = T)
mtext(resp$units[resp$resp == j], side = 2, outer = T)
title(main=toupper(j),outer=T)

dev.off()
}


#reponse by biomass for all habitats
#holds initial 'byyear' dataframe 
byyearbystrat <- with(byyear.sp, aggregate(cbind(algcons,allmacro, total_w, totalscar)~year+strat, FUN = 'sum'))
byyearbystrat$scarpercent <- byyearbystrat$totalscar/100

byyearhold <- byyearbystrat[byyearbystrat$year >= 1996,]

#replaces byyear with byyearhold after modifying by year hold
byyear <- byyearhold

byyear <- byyear[byyear$year>=1996,]

png('FL_all_habitats_grazing_response_by_biomass.png', width = 7, height = 3.5, units = 'in', res = 300, bg = "white")

par(mfrow = c(1,3), mar = c(2.5,4.5,2,2), oma = c(4,3,2,0), cex.lab = 1.4, cex.axis = 1.4)

for(j in c('algcons', 'allmacro', 'scarpercent')){
plot(byyearhold$total_w, byyearhold[,j], pch = 21, col = 'black', bg = 'green', main = j, xlim = c(0,max(byyearhold$total_w)),ylim = c(0,max(byyearhold[j])*1.05), xlab = '', ylab = resp$units[resp$resp == j])
mtext(paste('r2=',round(summary(lm(byyearhold[,j]~byyearhold$total_w))[[8]],3),sep = ''), cex = .7, line = -1.25, adj = 0.05) 
#mtext(resp$units[resp$resp == j], side = 2, line = 2.25,cex = 1.2, outer = F)
}
mtext('Total parrotfish biomass (g/m2)', side = 1, outer = T, cex = 1.1)
title(main='Grazing response by total parrotfish biomass',outer=T)

dev.off()


####Process STX 2012 data####
library(reshape)
stx <- read.table('STX_2012_fish_sizeclass.csv', header = T, sep = ',')
stx$survey_type <- NULL
names(stx)[names(stx) == 'fish_count'] <- 'num'
names(stx)[names(stx) == 'fish_size'] <- 'len'


stx$len <- with(stx, ifelse(len == '0--5', 2.5,
                      ifelse(len == '5--10', 7.5,
                        ifelse(len == '10--15', 12.5,
                          ifelse(len == '15--20', 17.5,
                            ifelse(len == '20--25', 22.5,
                              ifelse(len == '25--30', 27.5,
                                ifelse(len == '30--35', 32.5,
                                  ifelse(len == 'NULL', 'nulllen', as.character(stx$len))))))))))

spcodes <- read.table('Biogeo RVC Species Codes.csv', header = T, sep = ',')

stx1 <- merge(stx, spcodes, by = 'species_code', all.x = T, all.y = F)
#read spmetrics
#spmetrics <- read.table('grazing_metrics_single_table.csv', header = T, sep = ',') #FL Keys metrics
spmetrics <- read.table('grazing_metrics_single_table_STX.csv', header = T, sep = ',') #STX metrics

stx2 <- stx1[stx1$spcode %in% unique(spmetrics$spcode),]

#removes observation with >15 indiv >15cm; one obs from SDR
stx2 <- stx2[!(stx2$num > 15 & as.numeric(as.character(stx2$len)) > 15),]


#list of stations with lat lon
stx.latlon <- with(stx2, aggregate(cbind(latitude, longitude)~station_code+strata+zone_2+depth_2+habitat_2, FUN = mean))
stx.latlon$habtype <- with(stx.latlon, paste(zone_2, habitat_2, depth_2, sep = ' '))
allstrat <- unique(stx.latlon$habtype)

#with(stx.latlon, plot(longitude, latitude))
#with(stx.latlon, text(longitude, latitude, zone_2, cex = 0.5))

#ggplot(stx.latlon, aes(longitude, latitude, label = stx.latlon$zone_2)) + geom_text(aes(color = factor(stx.latlon$zone_2)))

#remove some unnecessary variables
removals <- c('region', 'survey_year', 'batch_code', 'survey_index', 'management', 'family', 'genus', 'scientific_name', 'common_name', 'trophic')
stx2[removals] <- NULL

vars <- names(stx2); vars <- vars[vars != c('num')]; vars <- vars[vars != c('len')]

#reshape to add zeroes to sizes with no data
stx.wide <- reshape(stx2, timevar = 'len', idvar = vars, direction = 'wide')
#splits dataset to add '0' values to NAs in length categories




lens.vars <- names(stx.wide)[(grep('num.', names(stx.wide)))]
lens.vars.nums <- gsub('num.','')

stx.wide.cat <- stx.wide[,!names(stx.wide) %in% lens.vars]
stx.wide.lens <- stx.wide[,lens.vars]
stx.wide.lens[is.na(stx.wide.lens)] <- 0

#reform the stx.wide dataset
stx.wide <- cbind(stx.wide.cat,stx.wide.lens)
#removes 'num.' from variable names
names(stx.wide) <- gsub('num.','', names(stx.wide))

#reshapes to 'long' form again
stx.l <- reshape(stx.wide, varying = names(stx.wide)[names(stx.wide) %in% lens.vars.nums], v.names = 'num', timevar = 'len', times = names(stx.wide)[names(stx.wide) %in% lens.vars.nums], direction = 'long')
row.names(stx.l) <- 1:nrow(stx.l)
stx.l$id <- NULL

stx.sp.stncode <- as.data.frame(matrix(nrow = length(unique(stx.l$spcode))*length(unique(stx.l$len)), ncol = length(unique(stx$station_code))+2))

names(stx.sp.stncode) <- c('spcode', 'len', unique(as.character(stx$station_code)))

stx.sp.stncode[,] <- 0
stx.sp.stncode$spcode <- unique(stx.l$spcode)
stx.sp.stncode$len <- rep(sort(unique(stx.l$len)),each = length(unique(stx.l$spcode))) 

stx.sp.stncode.l <- reshape(stx.sp.stncode, varying = names(stx.sp.stncode)[c(-1,-2)], v.names = 'num1',timevar = 'station_code', times = names(stx.sp.stncode)[c(-1,-2)], direction = 'long')
row.names(stx.sp.stncode.l) <- 1:nrow(stx.sp.stncode.l)
stx.sp.stncode.l$id <- NULL

stx.l.all <- merge(stx.sp.stncode.l, stx.l, by = c('spcode','station_code', 'len'), all.x = T, all.y = T)

stx.hab.vars <- as.data.frame(with(stx, table(station_code, strata, zone_2, depth_strata, depth_2, habitat_strata, habitat_2)))
stx.hab.vars <- stx.hab.vars[stx.hab.vars$Freq >0,]
stx.hab.vars$Freq <- NULL

stx.l.all[names(stx.hab.vars)[-1]] <- NULL

stx.l.all1 <- merge(stx.l.all, stx.hab.vars, by = 'station_code', all = T)

#writes stx long table
write.table(stx.l.all1, 'stx_long.csv', row.names = F, col.names = T, sep = ',')

#reads STX LONG table
stx.l <- read.table('stx_long.csv', header = T, sep = ',')

#reads RVC vs Biogeo sp codes
spcodes <- read.table('Biogeo RVC Species Codes.csv', header = T, sep = ',')
#adds RVC spcodes
#stx.l <- merge(stx.1, spcodes, by = 'species_code', all.x = T, all.y = F); rm(stx.1)
#subsets to parrotfish and surgeonfish only
#spmetrics <- read.table('grazing_metrics_single_table.csv', header = T, sep = ',') #FL Keys metrics
spmetrics <- read.table('grazing_metrics_single_table_STX.csv', header = T, sep = ',') #STX metrics

stx.pfish <- stx.l[stx.l$spcode %in% unique(spmetrics$spcode) & !stx.l$len == 'nulllen',]

#exludes all fish <= 10 cm
stx.pfish <- stx.pfish[as.numeric(as.character(stx.pfish$len)) > 10,]


#stx.pfish$len <- as.numeric(as.character(stx.pfish$len))

#removes all data where len > max len of all indiv seen
#stx.pfish <- stx.pfish[as.numeric(as.character(stx.pfish$len)) < (max(as.numeric(as.character(stx.pfish$len[stx.pfish$num >0])))+ 10),]


library(reshape)
#creates data frame with unique habitat code combinations; may not be needed, but good for reference
habitatcodes <- sort(unique(with(stx.pfish, paste(strata, zone, zone_2, depth_strata, depth_2, habitat_strata, habitat_2, sep = ' '))))
habitatcodes <- as.data.frame(colsplit(habitatcodes, split = ' ', c('strata', 'zone', 'zone_2', 'depth_strata', 'depth_2', 'habitat_strata', 'habitat_2')))

stx.pfish <- stx.pfish[!stx.pfish$habitat_2 == 'Bedrock',]
stx.pfish$zone_3 <- with(stx.pfish, ifelse(zone_2 == 'SARI', 'NE',
                                      ifelse(zone_2 == 'West', 'SW',
                                        ifelse(zone_2 == 'South', 'SW',
                                          ifelse(zone_2 == 'East', 'NE',
                                            ifelse(zone_2 == 'North', 'NE',
                                              ifelse(zone_2 == 'EEMP','NE', as.character(zone_2))))))))

#stx.pfish$zone_3 <- with(stx.pfish, ifelse(zone_2 == 'SARI', 'North',
#                                      ifelse(zone_2 == 'EEMP','East', as.character(zone_2))))


stx.pfish$newhab <- with(stx.pfish, ifelse(habitat_2 == 'SCR', 'Patchreef', as.character(habitat_2)))
stx.pfish$newhab <- with(stx.pfish, ifelse(newhab == 'Reef', 'Hi-rel_reef', as.character(newhab)))

#headers for stx analysis. Add 'zone_3' to include region of island
#stx.stn <- c('spcode', 'zone_3', 'newhab', 'station_code')
stx.stn <- c('spcode', 'zone_3', 'newhab', 'depth_2', 'station_code')
stx.hab <- stx.stn[1:length(stx.stn)-1]

stx.by.station <- aggregate(cbind(stx.pfish['num']), stx.pfish[c(stx.stn, 'len')], FUN = sum, na.rm = T) 
stx.by.station$strat <- with(stx.by.station, paste(zone_3, newhab, depth_2, sep = ' '))

#to include depth
#stx.by.station$newstrat <- with(stx.by.station, paste(strat, depth_2, sep = ' '))

#use strat to drop depth
#stx.by.station$strat <- with(stx.by.station, ifelse(zone_3 == 'West', 'West',
#                                              ifelse(strat == 'North Patchreef', 'North Pavement', as.character(strat))))
#use newstrat to use depth
#stx.by.station$newstrat <- with(stx.by.station, ifelse(strat == 'Buck Reef Shallow', 'Buck Reef',
#                                              ifelse(strat == 'Buck Reef Deep', 'Buck Reef',
#                                                ifelse(strat == 'East Reef Deep', 'East Reef',
#                                                  ifelse(strat == 'East Reef Shallow', 'East Reef',
#                                                    ifelse(strat == 'North Patchreef Deep', 'North Pavement',
#                                                      ifelse(strat == 'North Patchreef Shallow', 'North Pavement',
#                                                        ifelse(strat == 'North Pavement Deep', 'North Pavement',
#                                                          ifelse(strat == 'North Pavement Shallow', 'North Pavement',
#                                                            ifelse(strat == 'South Patchreef Shallow', 'South Pavement',
#                                                              ifelse(strat == 'South Patchreef Deep', 'South Pavement',
#                                                                ifelse(strat == 'South Reef Deep', 'South Reef',
#                                                                  ifelse(strat == 'South Reef Shallow', 'South Reef', 
#                                                                    ifelse(strat == 'West', 'West', 
#                                                                      ifelse(strat == 'South Pavement Deep', 'South Pavement',
#                                                                        ifelse(strat == 'South Pavement Shallow', 'South Pavement',as.character(strat)))))))))))))))))

#stx.by.station$strat <- stx.by.station$newstrat

stx.by.hab <- aggregate(cbind(stx.by.station['num']), stx.by.station[c('spcode', 'strat', 'len')], FUN = mean, na.rm = T)
#calculates n of samples per strat
stx.by.hab$n <- aggregate(cbind(stx.by.station['num']), stx.by.station[c('spcode', 'strat', 'len')], FUN = 'length')[,4]
#checks number of habclasses
table(with(stx.by.hab[stx.by.hab$num >0,], paste(zone_3, newhab, depth_2, sep = ' ')))

#stx.by.hab$strat <- with(stx.by.hab, paste(zone_3, newhab, depth_2, sep = ' ')) #NOT NEEDED

habs <- with(stx.by.hab, aggregate(n~strat, FUN = 'mean'))

#remove Buck Reef Shallow
stx.by.hab <- stx.by.hab[!stx.by.hab$strat == 'Buck Hi-rel_reef Shallow',]

#to use STX parameters
#spmetrics <- read.table('grazing_metrics_single_table.csv', header = T, sep = ',') #FL Keys metrics
spmetrics <- read.table('grazing_metrics_single_table_STX.csv', header = T, sep = ',') #STX metrics

stx.graz <- merge(stx.by.hab, spmetrics, by = 'spcode', all.x = T, all.y = F)
#stx.graz <- stx.graz[!stx.graz$len == 'nulllen',]
stx.graz$len <- as.numeric(as.character(stx.graz$len))

stx.graz$indiv_w <- with(stx.graz, lw_a*(len^lw_b)) #units: g/cyl
stx.graz$total_w <- with(stx.graz, indiv_w*num)/100 #units: g/m2
stx.graz$algcons <- with(stx.graz, (0.0342*(indiv_w^0.816))*num)/100 #units: g/m2; allometric relationship from Bruggeman
stx.graz$allmacro <- with(stx.graz, algcons*propmacro) #units: g/m2
stx.graz$brownmacro <- with(stx.graz, algcons*propbrown) #units: g/m2
stx.graz$scarsize <- with(stx.graz, bitescar_m*(len^2)) #units: cm2/bite; cut: *num*720)/100 so it is area per bite, no time, no density
stx.graz$biterate <- with(stx.graz, (biterate_m*len+biterate_b)*600*num) #units: bites/day/line; 600 min/per 10 hr day
stx.graz$totalscar <- with(stx.graz, (biterate*scarsize))/100 #units: cm2/day/m2
stx.graz$bioeros <- with(stx.graz, ((bioeros_v*(len^3))*biterate*2*365)/(100*1000)) ##using 2 g/cm3 to convert from volume grazed to g of bioerosion, converts to kg/yr/m2; 365 is days/yr; 100 is m2 per transect, 1000 is g per kg

stx.graz <- stx.graz[!stx.graz$spcode == 'SCA GUAC',]



####creates len freq dataset for num > 0
stx.graz.len.freq <- stx.graz[stx.graz$num >0, 1:5]
#divides by smallest value of n for each species
for(s in unique(stx.graz.len.freq$spcode)){
  stx.graz.len.freq$adjnum[stx.graz.len.freq$spcode == s] <- stx.graz.len.freq$num[stx.graz.len.freq$spcode == s]/(min(stx.graz.len.freq$num[stx.graz.len.freq$num > 0 & stx.graz.len.freq$spcode == s]))
}

png(paste('STX_size_freq_per_species.png', sep = ''), width = 5, height = 5, units = 'in', res = 300, bg = "white")

par(mfrow = c(3,3), mar = c(2.5,2,2,1.5), oma = c(0,2,2,0))

big3 <- c('SCA COEL', 'SCA COER', 'SCA GUAC')

for(s in c(unique(as.character(stx.graz.len.freq$spcode[!stx.graz.len.freq$spcode %in% big3])), big3)){
  tempmat <- stx.graz.len.freq[stx.graz.len.freq == s,]
  temphist <- rep(tempmat$len, round(tempmat$adjnum))
  hist(temphist, main = spmetrics$speciesabbrev[spmetrics$spcode == s], breaks = seq(0, 75, by = 5),
       freq = F, xlab = 'Length', ylab = 'Rel Freq.', cex.lab = 1, cex.axis = 0.75, col = 'gold3')
  mtext(paste('\u03bc','=', round(mean(temphist),1), sep=''), side = 3, adj = .95, line = -1, cex = 0.7)
}
dev.off()   


#summarizes by strata
stx.sp.bystrat1 <- with(stx.graz, aggregate(cbind(algcons, allmacro, brownmacro, total_w)~spcode+strat, FUN = 'sum', na.rm = T))

stx.sp.bystrat2 <- with(stx.graz, aggregate(cbind(totalscar, bioeros)~spcode+strat, FUN = 'sum', na.rm = T))
stx.sp.bystrat3 <- merge(stx.sp.bystrat1, stx.sp.bystrat2, by = c('spcode', 'strat'), all.x = T, all.y = F)
stx.sp.bystrat3$totalscar <- ifelse(is.na(stx.sp.bystrat3$totalscar), 0, stx.sp.bystrat3$totalscar)
stx.sp.bystrat3$bioeros <- ifelse(is.na(stx.sp.bystrat3$bioeros), 0, stx.sp.bystrat3$bioeros)
stx.sp.bystrat <- merge(stx.sp.bystrat3, spmetrics[c('spcode', 'speciesabbrev')], by = 'spcode')
stx.sp.bystrat$scarpercent <- stx.sp.bystrat$totalscar/100

stx.resp <- as.data.frame(cbind(c('algcons', 'totalscar', 'allmacro', 'brownmacro', 'total_w'),c('g of C/day/m2', 'cm2/day/m2', 'g of C/day/m2', 'g of C/day/m2', 'g/m2')))
names(stx.resp) <- c('resp', 'units')



#sp.bystrat.sd <- with(byyear.sp[byyear.sp$year >= 2010,], aggregate(cbind(algcons, allmacro, brownmacro, totalscar, total_w)~spcode+strat, FUN = 'sd', na.rm = T))

#inserts sd values for sp.bystrat.sd and converts to SE
#for(i in 1:5){
#  sp.bystrat[,i+7] <- sp.bystrat.sd[,i+2]/sqrt(3)  
#}
#names(sp.bystrat)[8:12] <- paste(names(sp.bystrat[3:7]), '_SE', sep = '')

#splits 'strat' variable back into region, habitat, and depth categories
for(i in 1:nrow(stx.sp.bystrat)){
stx.sp.bystrat$reg[i] <- strsplit(stx.sp.bystrat$strat[i], ' ')[[1]][1]
stx.sp.bystrat$habtype[i] <- strsplit(stx.sp.bystrat$strat[i], ' ')[[1]][2]
stx.sp.bystrat$depth[i] <- strsplit(stx.sp.bystrat$strat[i], ' ')[[1]][3]
}

#ggplot(stx.sp.bystrat, aes(habtype, total_w, fill = spcode)) + geom_bar(stat = 'identity') + facet_grid(depth~reg) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

####STX plots biomass stacked barplot####
##responses, units, and plotting places dataframe
resp <- as.data.frame(cbind(c('algcons', 'totalscar', 'allmacro', 'brownmacro', 'total_w', 'scarpercent', 'bioeros'),c('g of C' ~ day^-1 ~ m^-2, cm^2 ~ day^-1 ~ m^-2,'g of C' ~ day^-1 ~ m^-2, 'g of C' ~ day^-1 ~ m^-2, 'g'~m^-2,'% area grazed'~year^-1, 'kg'~m^-2~year^-1)))
names(resp) <- c('resp', 'units')
resp$places <- c(2.5, 2.5, 2.5, 2.5, 2, 2.25, 2)

#loop for all 4 responses
for(response in c('total_w', 'allmacro', 'scarpercent', 'bioeros')){
#response <- 'total_w'

plotmat <- reshape(stx.sp.bystrat[,c('speciesabbrev', 'strat', response)],timevar = 'strat', idvar = 'speciesabbrev', direction = 'wide')
names(plotmat) <- gsub(paste(as.character(response), '.', sep =''),'', names(plotmat))
row.names(plotmat) <- plotmat$speciesabbrev
plotmat$spcode <- NULL; plotmat$speciesabbrev <- NULL
plotmat <- as.matrix(plotmat)

png(paste('STX', response, 'stacked_by_hab_leg.png', sep = '_'), width = 10, height = 7, bg = 'white', res = 300, units = 'in')

par(mfrow = c(1,1), mar = c(12,8,5,11), oma = c(1,1,1,1), xpd = T, cex.lab = .8)
mycolors <- c(as.character(spmetrics$color[spmetrics$speciesabbrev %in% sort(unique(stx.sp.bystrat$speciesabbrev))]))

barcenters <- barplot(plotmat, col = mycolors, las = 2, cex.axis = 1.5, cex.names = 1.5, cex.lab = 1.5, main = response, xpd = F)

legend("topright", inset = c(-0.35, 0), rownames(plotmat), fill = mycolors, cex = 1.2)
title(ylab = resp$units[resp$resp == response], cex.lab = 2, line = resp$places[resp$resp == response]*1.75)

#barplot(plotmat[,-1], col = mycolors, las = 2, cex.axis =1, cex.names = 1, cex.lab = 1, main = response)#, legend = rownames(plotmat))
#legend(5,15, rownames(plotmat), fill = mycolors, cex = 1.2)

dev.off()
}

par(mar = c(12,6,4,2))
mycolors <- c('dark red', 'red', 'orange',' yellow', 'light green', 'dark green', 'light blue', 'blue', 'dark blue')
par(mfrow = c(1,1))
barplot(plotmat[,-1], col = mycolors, ylab = expression(paste('g m'^'-2')), las = 2, cex.axis =1.5, cex.names = 1.5, cex.lab = 1.5)#legend = rownames(plotmat),
#legend(5,15, rownames(plotmat), fill = mycolors, cex = 1.2)

dev.off()





##GGplot!
library(ggplot2)

png('STX_stacked_spbiomass_by_hab.png', width = 8, height = 6, bg = 'white', res = 300, units = 'in')

ggplot(stx.sp.bystrat[,c('speciesabbrev', 'strat','total_w')], aes(x = strat, y = total_w, fill = speciesabbrev), las = 2) +
  geom_bar(stat = 'identity') +
  xlab("") +
  ylab(expression(paste('g m'^'-2'))) +
  theme(legend.title=element_blank(),
        panel.background=element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15))

dev.off()

#plots responses by species.
for(j in stx.resp$resp){
  png(paste('Params.stx.sp_w_', j, '_per_strat.png', sep = ''), width = 8, height = 12, units = 'in', res = 300, bg = "white")
  
  par(mfrow = c(4,3), mar = c(2.5,3,2,2), oma = c(2,2,2,0))
  
  for(i in unique(stx.sp.bystrat$strat)) {
    barplot(stx.sp.bystrat[stx.sp.bystrat$strat == i,j], names.arg = stx.sp.bystrat[stx.sp.bystrat$strat == i,c('spcode')] , col = 'blue', ylim = c(0,max(stx.sp.bystrat[j])*1.05), xlab = "", las = 2, cex.names = 0.8)
    mtext(i, side = 3, line = -1.5)
  }
  title(main=toupper(j),outer=T)
  mtext(stx.resp$units[stx.resp$resp == j], side = 2, outer = T)
  
  dev.off()
}


####STX NMDS plots####
library(vegan)
library(MASS)

ordmat <- reshape(stx.sp.bystrat[,c('speciesabbrev', 'strat','total_w')],timevar = 'speciesabbrev', idvar = 'strat', direction = 'wide')
names(ordmat) <- gsub('total_w.','', names(ordmat))
row.names(ordmat) <- ordmat$strat
ordmat$strat <- NULL
ordmat <- as.matrix(ordmat)


png(paste('STX_NMDS_by_hab.png', sep = ''), width = 6, height = 6, units = 'in', res = 300, bg = "white")

par(mfrow = c(1,1), mar = c(2.5,3,2,2), oma = c(3,3,3,0))


  ord <- metaMDS(ordmat)
  
  #creates vectors from ord results
  xx <- as.vector(ord$points[,1])
  yy <- as.vector(ord$points[,2])
  sx <- as.vector(ord$species[,1])
  sy <- as.vector(ord$species[,2])
  
  #plots NMDS plots, 
  plot(xx,yy, xlim = c(min(xx,sx)*1.3, max(xx,sx)*1.5), ylim = c(min(yy,sy)*1.2, max(yy,sy)*1.2), pch = 21, col = 'dark orange', bg = 'orange', cex = 1.4, xlab = 'MDSA1', ylab = 'MDSA2')
  text(xx,yy,labels = row.names(ordmat), pos = 3, offset = 0.5, cex = 1.1)
#  text(sx, sy, labels = sort(unique(spmetrics$speciesabbrev)), col = 'blue') #omit this just to show 
  mtext(paste('stress=', round(ord[[22]],4),sep = ''),side = 3, line = 0, adj = 0.95, cex = 0.9)
  dev.off()
  








barplot(stx.sp.bystrat[stx.sp.bystrat$strat == i,j], names.arg = stx.sp.bystrat[stx.sp.bystrat$strat == i,c('spcode')] , col = 'blue', ylim = c(0,max(stx.sp.bystrat[j])*1.05), xlab = "", las = 2, cex.names = 0.8)

barplot(stx.sp.bystrat)


stx.mean1 <- with(stx.graz, aggregate(cbind(algcons, allmacro, brownmacro, total_w)~newstrat, FUN = 'sum', na.rm = T))

stx.mean2 <- with(stx.graz, aggregate(cbind(totalscar)~newstrat, FUN = 'sum', na.rm = T))

stx.mean <- merge(stx.mean1, stx.mean2, by = c('newstrat'), all.x = T, all.y = F)

write.table(stx.mean, 'stx_grazing_data1.csv', sep = ',', row.names = F, col.names = T)



#new calc of stx.mean
stx.mean1 <- with(stx.sp.bystrat, aggregate(cbind(algcons, allmacro, brownmacro, total_w)~strat, FUN = 'sum'))
stx.mean2 <- with(stx.sp.bystrat, aggregate(cbind(totalscar)~strat, FUN = 'sum'))
stx.mean <- merge(stx.mean1, stx.mean2, by = 'strat', all.x = T, all.y = F)
stx.mean$scarpercent <- stx.mean$totalscar/100
stx.mean <- stx.mean[stx.mean$scarpercent < 12,]

png('STX_all_habitats_grazing_response_by_biomass.png', width = 7, height = 3.5, units = 'in', res = 300, bg = "white")

par(mfrow = c(1,3), mar = c(2.5,4.5,2,2), oma = c(4,3,2,0), cex.lab = 1.4, cex.axis = 1.4)

for(k in c('algcons', 'allmacro', 'scarpercent')){

  plot(stx.mean$total_w, stx.mean[,k], pch = 21, col = 'black', bg = 'dark orange', cex =1.2, main = k, xlim = c(0,33),ylim = c(0,max(stx.mean[k])*1.05), xlab = '', ylab = resp$units[resp$resp == k])
  mtext(paste('r2=',round(summary(lm(stx.mean[,k]~stx.mean$total_w))[[8]],3),sep = ''), cex = .7, line = -1.25, adj = 0.05)
#mtext(resp$units[resp$resp == k], side = 2, line = 2.25,cex = 0.75, outer = F)  
  }
mtext('Total parrotfish biomass (g/m2)', side = 1, outer = T)
title(main='Grazing response by total parrotfish biomass',outer=T)

dev.off()  
  
  
####plot both FL and STX

png('STX_all_habitats_grazing_response_by_biomass.png', width = 7, height = 3.5, units = 'in', res = 300, bg = "white")

par(mfrow = c(1,3), mar = c(2.5,4.5,2,2), oma = c(4,3,2,0), cex.lab = 1.4, cex.axis = 1.4)

for(k in c('algcons', 'allmacro', 'scarpercent')){
  
  plot(byyearhold$total_w, byyearhold[,k], pch = 21, col = 'white', bg = 'white', cex =1.1, main = k, xlim = c(0,39),ylim = c(0,max(c(stx.mean[[k]],byyearhold[[k]]))*1.05), xlab = '', ylab = resp$units[resp$resp == k])
  #mtext(paste('r2=',round(summary(lm(stx.mean[,k]~stx.mean$total_w))[[8]],3),sep = ''), cex = .7, line = -1.25, adj = 0.05)

  points(stx.mean$total_w, stx.mean[,k], pch = 21, col = 'black', bg = 'dark orange', main = '', cex = 1.3)
   # mtext(paste('r2=',round(summary(lm(byyear[,j]~byyear$total_w))[[8]],3),sep = ''), cex = .7, line = -1.25, adj = 0.05) 
    
  
  
    #mtext(resp$units[resp$resp == k], side = 2, line = 2.25,cex = 0.75, outer = F)  
}
mtext('Total parrotfish biomass (g/m2)', side = 1, outer = T)
title(main='Grazing response by total parrotfish biomass',outer=T)

dev.off()  


####Grazing metrics by species####

stxmetrics <- read.table('grazing_metrics_single_table_STX.csv', header = T, sep = ',') #STX metrics
flmetrics <- read.table('grazing_metrics_single_table.csv', header = T, sep = ',')

allmetrics1 <- merge(flmetrics,stxmetrics[c('spcode', 'propmacro','propbrown','bitescar_m', 'biterate_m', 'biterate_b')], by = 'spcode')

fl <- with(grazing1[grazing1$num > 0 & grazing1$year >= 2010,], aggregate(len~spcode, FUN = 'max', na.rm =T))
stc <- with(stx.graz[stx.graz$num > 0,], aggregate(len~spcode, FUN = 'max', na.rm =T))
maxsize <- merge(fl,stc, by = 'spcode', all.x = T, all.y = F)

allmetrics <- merge(allmetrics1, maxsize, by = 'spcode', all.x = y, all.y = F)

names(allmetrics) <- gsub('.x','.fl', names(allmetrics))
names(allmetrics) <- gsub('.y','.stx', names(allmetrics))


lenmat <- as.data.frame(sort(c(rep(5:70,times = length(unique(allmetrics$spcode))))))
names(lenmat) <- 'len'
lenmat$spcode <- unique(sort(allmetrics$spcode))

metricsmat <- merge(lenmat, allmetrics, by = 'spcode', all.x = T, all.y = F)

metricsmat <- metricsmat[order(metricsmat$len),]

#need to generate size response plots
metricsmat$indiv_w <- with(metricsmat,lw_a*(len^lw_b))
metricsmat$indiv_w.fl <- with(metricsmat, ifelse(len > len.fl, 'NA',lw_a*(len^lw_b)))
metricsmat$indiv_w.stx <- with(metricsmat, ifelse(len > len.stx, 'NA',lw_a*(len^lw_b)))
metricsmat$max_w.fl <- with(metricsmat, lw_a*(len.fl^lw_b))
metricsmat$max_w.fl <- with(metricsmat, lw_a*(len.stx^lw_b))
metricsmat$algcons <- with(metricsmat, (0.0342*(indiv_w^0.816)))
metricsmat$allmacro.fl <- with(metricsmat, algcons*propmacro.fl)
metricsmat$allmacro.stx <- with(metricsmat, algcons*propmacro.stx)
#metricsmat$brownmacro <- with(metricsmat, algcons*propbrown)
metricsmat$scarsize.fl <- with(metricsmat, bitescar_m.fl*(len^2)*720)
metricsmat$scarsize.stx <- with(metricsmat, bitescar_m.stx*(len^2)*720)
metricsmat$biterate.fl <- with(metricsmat, (biterate_m.fl*len+biterate_b.fl)*720)
metricsmat$biterate.stx <- with(metricsmat, (biterate_m.stx*len+biterate_b.stx)*720)
metricsmat$totalscar.fl <- with(metricsmat, (biterate.fl*scarsize.fl))
metricsmat$totalscar.stx <- with(metricsmat, (biterate.stx*scarsize.stx))

write.table(metricsmat, 'metricsmatrix.csv', sep = ',', row.names = F, col.names = T)

a <- as.data.frame(with(metricsmat, table(spcode, len.fl)))
a <- a[a$Freq>0,-3]
b <- as.data.frame(with(metricsmat, table(spcode, len.stx)))
b <- b[b$Freq > 0,-3]
maxlen <- merge(a,b, by = 'spcode', all.x = T, all.y. = F)


s <- 'SPA VIRI'

r <- 'allmacro'

png(paste('metricsplot', s, r, '.png', sep = '_'), width = 4, height = 3, units = 'in', res = 300, bg = 'white')

par(mfrow = c(1,1), mar = c(4,4,3,1), cex.lab = 1.2, cex.axis = 1.2)

with(metricsmat[metricsmat$spcode == s & metricsmat$len.fl <= as.numeric(as.character(maxlen$len.fl[maxlen$spcode == s])),], plot(indiv_w, allmacro.fl, type = 'l', lwd = 3, col = 'blue', main = unique(metricsmat$speciesabbrev[metricsmat$spcode == s]), xlab = 'body mass (g)', ylab = stx.resp$units[stx.resp$resp == r]))
with(metricsmat[metricsmat$spcode == s & metricsmat$len.stx <= as.numeric(as.character(maxlen$len.stx[maxlen$spcode == s])),], lines(indiv_w, allmacro.stx, lwd = 3, col = 'red', main = unique(metricsmat$speciesabbrev[metricsmat$spcode == s],), xlab = 'body mass (g)', ylab = stx.resp$units[stx.resp$resp == r]))

dev.off()


###TO HERE
#checks parrotfish species
with(stx.l[stx.l$family == 'Scaridae',], unique(species_code))

x<- stx.l[is.na(stx.l$spcode) & stx.l$num >0,]


### RVC Parrotfish data processing junk
dim(denssize.sum.psu.1[denssize.sum.psu.1$num > 0 & is.na(denssize.sum.psu.1$len),])



#junk testing below
dim(denssize.sum.psu[denssize.sum.psu$len <= 0 | denssize.sum.psu$num == 0,])


sum(denssize.sum.psu$wsize[denssize.sum.psu$year == 1980 & denssize.sum.psu$spcode == 'SPA AURO' & denssize.sum.psu$psu == '044U'])
#


psu <- as.data.frame(table(denssize.psu$num[denssize.psu$len<=0 & !denssize.psu$num == 0]))
names(psu) <- c('density','freq')

#negative length troubleshooting
xx <- with(denssize.psu[denssize.psu$num > 0,], aggregate(num~year+spcode+prot+psu, FUN = length))

yy <- with(denssize.psu[denssize.psu$len <= 0 & denssize.psu$num >0,], aggregate(num~year+spcode+prot+psu, FUN = length))

zz <- merge(xx, yy, by = c('year', 'spcode', 'prot', 'psu'), all.x = T, all.y = F)
zz$num.y[is.na(zz$num.y)] <- 0
zz$prop <- zz$num.y/zz$num.x

length(zz$prop[zz$prop > 0])
length(zz$year[zz$prop == 1  & zz$num.y == 1])



length(zz$prop[zz$prop<0.05])/length(zz$prop)

hist(zz$num.y, breaks = seq(0,max(zz$num.y)+1,.005))
unique(zz$num.y)


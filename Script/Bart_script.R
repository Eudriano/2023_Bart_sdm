#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Eudriano Costa
# Manuscrispt: Using BART to predict Callinectes sapidus distribution 
# SDM - Species Distribution Models  
# R version 4.0.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

#Cuidado
save.image("2023_BART.RData")

setwd("~/efc/R_Stat/2022_NEMA/SDM_Bayes/1_Species")
load(file = "2023_.RData") # save work space


#--------------------------- Labraries/ packages 
packs <-c('rgdal','spdep', 'tmap','tidyverse','dismo','maptools','sf','report', 
          'performance','sp','rgeos','ggplot2','plotly','raster','tidyr',
          'dplyr','sdm', 'rworldmap','usdm', 'mapview', 
          'evaluate','corrplot','GGally','embarcadero')

lapply(packs, FUN = function(X) {
                    do.call("require", list(X)) 
                                  })
(.packages()) # Check loaded packages

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                           Species data                                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/efc/R_Stat/2022_NEMA/SDM_Bayes/1_Species")
spg<-read.csv("Csapidus_filtered.csv",dec=".",sep=",",header=T)
head(spg, 20) 
nrow(spg) # 41.122 occurrences

#--------------------------- Duplicated coordinates data 
dups2 <- duplicated(spg[, c("Latitude","Longitude")]);dups2
sum(dups2) # duplicated coordinates = 13.470

spg <- spg[!dups2, ]
nrow(spg) # final data = 27.652 occurrences

#--------------------------- Covert data frame to spatial data
spg1 <- SpatialPointsDataFrame(spg[, 1:2], data.frame(spg[, 3]))
names(spg1@data) <- 'Presence'
nrow(spg1)
head(spg1)

save(spg1, file='spg1.rds')
load('spg1.rds')

spg1<-as.data.frame(spg1)
Sapidus<-spg1
Sapidus<-Sapidus[, c('Longitude', 'Latitude')]
Sapidus$Presence = 1
nrow(Sapidus) #27652

write.csv(Sapidus,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

writeOGR(EEZ, driver="ESRI Shapefile",
         dsn = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes" , layer = "EEZFinal")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Set color pallette                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cl<-colorRampPalette(c("#3E49BB","#3498DB","yellow",    
                                "orange", "red", "darkred"))(200)
                                
save(cl, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/cl.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/cl.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                      PRESENT- 2000 to 2014: Covariates                       #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
Present_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/Present", 
                         pattern="\\.tif$", full.names=T)
Present_Cov

#--------------------------- Stack rasters and set crs
Present_Cov <- stack(Present_Cov)
Present_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(Present_Cov[[3]],col=cl)

#--------------------------- Crop both dataset manually
#e <- drawExtent() # ATTETION 
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
Present_bioc<- crop(Present_Cov,e)

#---------------------------  Crop by EEZ
Start <- Sys.time()
masked <- mask(x = Present_bioc, mask = EEZ)
plot(masked, col=cl)
Present_bioc <- crop(x = masked, y = extent(EEZ))
plot(Present_bioc, col=cl)
# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 3.6 min

#---------------------------  Rename and CRS
Present_covs <- Present_bioc # covariates
Present_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots
plot(stack(Present_covs),col=cl) #cropped files

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Covs.tiff",
     width=8, height=4, units="in", res=600)

plot(Present_covs, col=cl)

dev.off()

#---------------------------  Rasterizations and point by grid cell
Sapidus <- SpatialPointsDataFrame(Sapidus[, 1:2], data.frame(Presence = Sapidus[, 3]))
tmp <- rasterize(Sapidus, Present_covs[[1]], field = "Presence", fun = "min")
pts.sp1 <- rasterToPoints(tmp,
                          fun = function(x) {
                            x > 0
                          })
nrow(pts.sp1)# 2121

projection(Sapidus) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus) 
st_crs(Sapidus)

write.csv(pts.sp1,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Presence_Rasterization.csv')


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
pres.cov <- raster::extract(Present_covs, pts.sp1[, 1:2])
pres.cov <- na.omit(pres.cov)
head(pres.cov) #602, 668
nrow(pres.cov)
nrow(pts.sp1) #2121
which(is.na(pres.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsences                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence <- randomPoints(Present_covs, nrow(pres.cov))
#absence1 <- randomPoints(Present_covs, p=2121, 2121 )
#plot(absence1)

plot(absence)


nrow(absence)

abs.cov <- raster::extract(Present_covs, absence) #668

write.csv(absence,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')
#write.csv(absence1,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence1_0.csv')

#--------------------------- Code the response
pres.cov <- data.frame(pres.cov)
pres.cov$Sapidus <- 1
abs.cov <- data.frame(abs.cov)
abs.cov$Sapidus <- 0

#--------------------------- And one to bind them
Present_all.cov <- rbind(pres.cov, abs.cov)
Present_all.cov <- Present_all.cov[complete.cases(Present_all.cov), ]
nrow(Present_all.cov)
head(Present_all.cov)
which(is.na(Present_all.cov))

Present_xvars <- names(Present_all.cov)[!(names(Present_all.cov) == 'Sapidus')]
Present_xvars

# Save files
write.csv(Present_all.cov,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Present_all.cov.csv')

save(Present_all.cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_all.cov.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_all.cov.rds') 

save(Present_xvars, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_xvars.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_xvars.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                          Correlation and collinearity                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Present_vif<-vifstep(Present_all.cov[,1:4])
Present_vif

library(GGally)
ggpairs(Present_all.cov[,1:4])

library(ppcor)
pcor(Present_all.cov[,1:4], method = "pearson")

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(Present_all.cov[,-7]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cchf = Sapidus
#--------------------------- Variable selection and returns the model
Present.model <- bart.step(x.data = Present_all.cov[, Present_xvars],
                           y.data = Present_all.cov[, 'Sapidus'],
                           full = TRUE,
                           quiet = TRUE)

save(Present.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_model.rds') 

# Rugosity verification : Rug was not recommended
dbarts::pd2bart(Present.model,plotquants=T,
                xind = c('Pres_Temp', 'Pres_Rugosity'),
                pl = TRUE)

#----------------------------------------TEsting
#https://github.com/cjcarlson/embarcadero/blob/master/vignettes/virtualbart.Rmd
sdm_Present<- bart(y.train=Present_all.cov[, 'Sapidus'],
                      x.train=Present_all.cov[, Present_xvars],
                      keeptrees = TRUE)
summary(sdm_Present)

map <- predict(sdm_Present, Present_covs, quiet=TRUE)
map <- predict(sdm_Present, Present_covs,quantiles=c(0.025, 0.975), quiet=TRUE)

par(mfrow=c(2,2))
par(mar=c(2,1,2,5))
plot(map[[1]], 'Posterior mean', col=cl,
     box=F, axes=F)
plot(map[[2]], 'Lower 95% CI bound',col=cl, 
     box=F, axes=F)
plot(map[[3]], 'Upper 95% CI bound', col=cl,
     box=F, axes=F)
plot(map[[3]]-map[[2]], 'Credible interval width',col=cl, 
     box=F, axes=F)


#-----------------------------------------

#--------------------------- Retune: cross-validation
RET_Present.model<-retune(Present.model, reps = 10) 
summary(RET_Present.model)

save(RET_Present.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_retune.model.rds') 

#--------------------------- Spatial prediction
Present_Pred.map <- predict(
  object = RET_Present.model,
  x.layers = Present_covs,
  quantiles = c(0.025, 0.975),
  splitby = 20,
  quiet = TRUE)

save(Present_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_Pred.map.rds') 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------Present: PLOT probs.: mean, 2.5% , 97.5% and uncertainty
#>>>>>>>>>>>>>>>>>>>>>>>>>>>


#~~~~~~~~~~ Present_Mean predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_MeanPred.tiff",
     width=6, height=4, units="in", res=600)
plot(Present_Pred.map[[1]],col=cl,
     box = F,
     axes = F,
     main = 'Present: mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
points(pts.sp1[!is.na(raster::extract(Present_Pred.map[[1]], 
                                      SpatialPoints(pts.sp1[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ Present_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(Present_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'Present: 2.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ Present_97.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_97.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(Present_Pred.map[[3]],col=cl,
     box = F,
     axes = F,
     main = 'Present: 97.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ Present_Uncertainty predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Uncertainty.tiff",
     width=6, height=4, units="in", res=600)
plot(Present_Pred.map[[3]]-Present_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'Present Uncertainty',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Posterior width', side=2, line=1.3))
dev.off()

#----- Present: save rasters
writeRaster(Present_Pred.map[[1]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.MeanPred.tif"
            ,overwrite=TRUE)
writeRaster(Present_Pred.map[[2]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.2.5Pred.tif"
            ,overwrite=TRUE)
writeRaster(Present_Pred.map[[3]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.97.5Pred.tif"
            ,overwrite=TRUE)
Present_Uncert<-Present_Pred.map[[3]]-Present_Pred.map[[2]]# Uncertainty
writeRaster(Present_Uncert, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.Uncertainty.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present                           ANALYTICS                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~ Load modified functions
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/varimpEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/summaryEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/partialsEU.R")

#~~~~~~~~~~ Varimp diagonal
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_VarimpDiag.tiff",
     width=6, height=4, units="in", res=600)
varimp.diag(Present_all.cov[,Present_xvars], Present_all.cov[,'Sapidus'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
#varimp(RET_Present.model, plots=TRUE)
varimpEU(RET_Present.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_Present.model)
dev.off()

library(tidytreatment)
posterior_fitted <- fitted_draws(RET_Present.model, value = "fit", include_newdata = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present              TWO-DIMENSIONAL NICHES                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Presente Temperature ~ salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Temp_Sal.tiff",
     width=8.5, height=4, units="in", res=600)
dbarts::pd2bart(RET_Present.model,plotquants=T, cexlab= 1.5, axes = F,
                xind = c('Pres_Temp', 'Pres_Sal'), contour.color ='black',
                pl = TRUE)
dev.off()

# Presente Temperature ~ current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Temp_Cur.tiff",
     width=8.5, height=4, units="in", res=600)
dbarts::pd2bart(RET_Present.model,plotquants=T,
                xind = c('Pres_Temp', 'Pres_Cur'),
                pl = TRUE)
dev.off()

# Presente Temperature ~ PrimProd
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Temp_PrimProd.tiff",
     width=8.5, height=4, units="in", res=600)
dbarts::pd2bart(RET_Present.model,plotquants=T,
                xind = c('Pres_Temp', 'Pres_PrimProd'),
                pl = TRUE)
dev.off()

# Presente Temperature ~ pH
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Temp_pH.tiff",
     width=8.5, height=4, units="in", res=600)
dbarts::pd2bart(RET_Present.model,plotquants=T,
                xind = c('Pres_Temp', 'Pres_pH'),
                pl = TRUE)
dev.off()

# Presente Temperature ~ Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_Temp_Rugosity.tiff",
     width=8.5, height=4, units="in", res=600)
dbarts::pd2bart(RET_Present.model,plotquants=T,
                xind = c('Pres_Temp', 'Pres_Rugosity'),
                pl = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present              PARTIAL PLOTS: Binary ~ Predictors                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Presente all variables
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_AllEU.tiff",
     width=13, height=9.5, units="in", res=600)
partialEU(RET_Present.model,
          c('Pres_Temp','Pres_Sal','Pres_Cur'),#,'Pres_Rugosity'
          trace = FALSE,
          ci = TRUE,
          panel = TRUE,
          smooth = 5)
dev.off()


# Present Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_Temp.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_Temp',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

# Presente Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_Sal.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_Sal',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

# Presente Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_Cur.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_Cur',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

# Primary porductivity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_PrimProd.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_PrimProd',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

# pH
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_pH.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_pH',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_Rugosity.tiff",
     width=6, height=6, units="in", res=600)
partialEU(RET_Present.model,
          'Pres_Rugosity',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 5)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Present Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatTemp.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.Temp <- spartial(RET_Present.model, Present_covs, 
                            x.vars = 'Pres_Temp', equal = TRUE)
plot(Present_sp.Temp,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Temp',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.Temp, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatTemp.tif"
            ,overwrite=TRUE)

# Present Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatSal.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.Sal <- spartial(RET_Present.model, Present_covs, 
                           x.vars = 'Pres_Sal', equal = TRUE)
plot(Present_sp.Sal,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Sal',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.Sal, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatSal.tif"
            ,overwrite=TRUE)

# Present Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatCur.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.Cur <- spartial(RET_Present.model, Present_covs, 
                           x.vars = 'Pres_Cur', equal = TRUE)
plot(Present_sp.Cur,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Cur',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.Cur, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatCur.tif"
            ,overwrite=TRUE)

# Primary porductivity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatPrimProd.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.PrimProd <- spartial(RET_Present.model, Present_covs, 
                                x.vars = 'Pres_PrimProd', equal = TRUE)
plot(Present_sp.PrimProd,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: PrimProd',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.PrimProd, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatPrimProd.tif"
            ,overwrite=TRUE)

# pH
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatpH.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.pH <- spartial(RET_Present.model, Present_covs, 
                          x.vars = 'Pres_pH', equal = TRUE)
plot(Present_sp.pH,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: PrimProd',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.pH, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatpH.tif"
            ,overwrite=TRUE)

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/Present_ParPlot_EffSpatRugosity.tiff",
     width=6, height=4, units="in", res=600)
Present_sp.Rugosity <- spartial(RET_Present.model, Present_covs, 
                                x.vars = 'Pres_Rugosity', equal = TRUE)
plot(Present_sp.Rugosity,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Rugosity',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(Present_sp.Rugosity, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/Present.EffSpatRugosity.tif"
            ,overwrite=TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                      RCP4.5 - 2050 : Covariates                              #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP45_2050_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP45_2050", 
                            pattern="\\.tif$", full.names=T)
RCP45_2050_Cov

#--------------------------- Stack rasters and set crs
RCP45_2050_Cov <- stack(RCP45_2050_Cov)
RCP45_2050_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP45_2050_Cov[[4]],col=cl)
plot(Sapidus,add=T, pch=19,cex=0.3, col="black")

#--------------------------- Crop both dataset manually
#e <- drawExtent() # ATTETION 
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP45_2050_bioc<- crop(RCP45_2050_Cov,e)

#---------------------------  Crop by EEZ
Start <- Sys.time()
masked <- mask(x = RCP45_2050_bioc, mask = EEZ)
plot(masked, col=cl)
RCP45_2050_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP45_2050_bioc, col=cl)
# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 15.2min

#---------------------------  Rename and CRS
RCP45_2050_covs <- RCP45_2050_bioc # covariates
RCP45_2050_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots
plot(stack(RCP45_2050_covs),col=cl) #cropped files


save(RCP45_2050_covs, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_Covs.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_Covs.rds') 

#---------------------------  Rasterizations and point by grid cell
Sapidus <-as.data.frame(Sapidus)
nrow(Sapidus)
Sapidus <- SpatialPointsDataFrame(Sapidus[, 2:3], data.frame(RCP45_2050ence = Sapidus[, 1]))
tmp <- rasterize(Sapidus, RCP45_2050_covs[[1]], field = "RCP45_2050ence", fun = "min")
pts.sp1 <- rasterToPoints(tmp,
                          fun = function(x) {
                            x > 0
                          })
nrow(pts.sp1)# 2130

projection(Sapidus) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus) 
st_crs(Sapidus)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract RCP45_2050ence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract RCP45_2050ence values
RCP45_2050.cov <- raster::extract(RCP45_2050_covs, pts.sp1[, 1:2])
RCP45_2050.cov <- na.omit(RCP45_2050.cov)
head(RCP45_2050.cov) #602
nrow(RCP45_2050.cov)
nrow(pts.sp1)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsences                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence <- randomPoints(RCP45_2050_covs, nrow(RCP45_2050.cov))
abs.cov <- raster::extract(RCP45_2050_covs, absence)
plot(absence)
#--------------------------- Code the response
RCP45_2050.cov <- data.frame(RCP45_2050.cov)
RCP45_2050.cov$Sapidus <- 1
abs.cov <- data.frame(abs.cov)
abs.cov$Sapidus <- 0

#--------------------------- And one to bind them
RCP45_2050_all.cov <- rbind(RCP45_2050.cov, abs.cov)
RCP45_2050_all.cov <- RCP45_2050_all.cov[complete.cases(RCP45_2050_all.cov), ]
nrow(RCP45_2050_all.cov)
head(RCP45_2050_all.cov)

RCP45_2050_xvars <- names(RCP45_2050_all.cov)[!(names(RCP45_2050_all.cov) == 'Sapidus')]
RCP45_2050_xvars

# Save files
write.csv(RCP45_2050_all.cov,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP45_2050_all.cov.csv')


save(RCP45_2050_all.cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_all.cov.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_all.cov.rds') 

save(RCP45_2050_xvars, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_xvars.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_xvars.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                          Correlation and collinearity                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RCP45_2050_vif<-vifstep(RCP45_2050_all.cov[,-5])
RCP45_2050_vif

library(GGally)
ggpairs(RCP45_2050_all.cov[,-5])

library(ppcor)
pcor(RCP45_2050_all.cov[,1:3], method = "pearson")

corrplot.mixed(RCP45_2050_all.cov[,1:3], lower.col = 'black', number.cex = .7)

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP45_2050_all.cov[,-5]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cchf = Sapidus
#--------------------------- Variable selection and returns the model
RCP45_2050.model <- bart.step(x.data = RCP45_2050_all.cov[, RCP45_2050_xvars],
                              y.data = RCP45_2050_all.cov[, 'Sapidus'],
                              full = TRUE,
                              quiet = TRUE)



save(RCP45_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP45_2050.model<-retune(RCP45_2050.model, reps = 10) 
summary(RET_RCP45_2050.model)
summary(RCP45_2050.model)

save(RET_RCP45_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_retune.model.rds') 

#--------------------------- Spatial prediction
RCP45_2050_Pred.map <- predict(
  object = RET_RCP45_2050.model,
  x.layers = RCP45_2050_covs,
  quantiles = c(0.025, 0.975),
  splitby = 20,
  quiet = TRUE)

save(RCP45_2050_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_Pred.map.rds') 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------RCP45_2050: PLOT probs.: mean, 2.5% , 97.5% and uncertainty
#>>>>>>>>>>>>>>>>>>>>>>>>>>>

#~~~~~~~~~~ RCP45_2050_Mean predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_MeanPred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2050_Pred.map[[1]],col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2050: mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
points(pts.sp1[!is.na(raster::extract(RCP45_2050_Pred.map[[1]], 
                                      SpatialPoints(pts.sp1[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ RCP45_2050_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2050_Pred.map[[2]], col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2050: 2.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()


#~~~~~~~~~~ RCP45_2050_97.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_97.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2050_Pred.map[[3]], col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2050: 97.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ RCP45_2050_Uncertainty predict

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Uncertainty.tiff",
     width=6, height=4, units="in", res=600)

plot(RCP45_2050_Pred.map[[3]]-RCP45_2050_Pred.map[[2]], col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2050 Uncertainty',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Posterior width', side=2, line=1.3))
dev.off()

#----- RCP45_2050: save rasters
writeRaster(RCP45_2050_Pred.map[[1]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.MeanPred.tif"
            ,overwrite=TRUE)
writeRaster(RCP45_2050_Pred.map[[2]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.2.5Pred.tif"
            ,overwrite=TRUE)
writeRaster(RCP45_2050_Pred.map[[3]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.97.5Pred.tif"
            ,overwrite=TRUE)
RCP45_2050_Uncert<-RCP45_2050_Pred.map[[3]]-RCP45_2050_Pred.map[[2]]# Uncertainty
writeRaster(RCP45_2050_Uncert, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.Uncertainty.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050                           ANALYTICS                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~ Load modified functions
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/varimpEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/summaryEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/partialsEU.R")

#~~~~~~~~~~ Varimp diagonal
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_VarimpDiag.tiff",
     width=6, height=4, units="in", res=600)
varimp.diag(RCP45_2050_all.cov[,RCP45_2050_xvars], RCP45_2050_all.cov[,'Sapidus'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
#varimp(RET_RCP45_2050.model, plots=TRUE)
varimpEU(RET_RCP45_2050.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_RCP45_2050.model)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050              TWO-DIMENSIONAL NICHES                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050e Temperature ~ salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Temp_Sal.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2050.model,plotquants=T,
                xind = c('RCP45_2050_Temp', 'RCP45_2050_Sal'),
                pl = TRUE)
dev.off()

# RCP45_2050e Temperature ~ current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Temp_Cur.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2050.model,plotquants=T,
                xind = c('RCP45_2050_Temp', 'RCP45_2050_Cur'),
                pl = TRUE)
dev.off()

# RCP45_2050e Temperature ~ PrimProd
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Temp_PrimProd.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2050.model,plotquants=T,
                xind = c('RCP45_2050_Temp', 'RCP45_2050_PrimProd'),
                pl = TRUE)
dev.off()

# RCP45_2050e Temperature ~ pH
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Temp_pH.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2050.model,plotquants=T,
                xind = c('RCP45_2050_Temp', 'RCP45_2050_pH'),
                pl = TRUE)
dev.off()

# RCP45_2050e Temperature ~ Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Temp_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2050.model,plotquants=T,
                xind = c('RCP45_2050_Temp', 'RCP45_2050_Rugosity'),
                pl = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050           PARTIAL PLOTS: Binary ~ Predictors                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050e all variables
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_All.tiff",
     width=8, height=5, units="in", res=600)
partial(RET_RCP45_2050.model,
        c('RCP45_2050_Temp','RCP45_2050_Sal','RCP45_2050_Cur','RCP45_2050_Rugosity'),
        trace = FALSE,
        ci = TRUE,
        panel = TRUE,
        smooth = 8)
dev.off()

# RCP45_2050e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_Temp.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2050.model,
          'RCP45_2050_Temp',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP45_2050e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_Sal.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2050.model,
          'RCP45_2050_Sal',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP45_2050e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_Cur.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2050.model,
          'RCP45_2050_Cur',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2050.model,
          'RCP45_2050_Rugosity',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP45_2050e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_EffSpatTemp.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2050_sp.Temp <- spartial(RET_RCP45_2050.model, RCP45_2050_covs, 
                               x.vars = 'RCP45_2050_Temp', equal = TRUE)
plot(RCP45_2050_sp.Temp,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Temp',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2050_sp.Temp, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.EffSpatTemp.tif"
            ,overwrite=TRUE)

# RCP45_2050e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_EffSpatSal.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2050_sp.Sal <- spartial(RET_RCP45_2050.model, RCP45_2050_covs, 
                              x.vars = 'RCP45_2050_Sal', equal = TRUE)
plot(RCP45_2050_sp.Sal,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Sal',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2050_sp.Sal, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.EffSpatSal.tif"
            ,overwrite=TRUE)

# RCP45_2050e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_EffSpatCur.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2050_sp.Cur <- spartial(RET_RCP45_2050.model, RCP45_2050_covs, 
                              x.vars = 'RCP45_2050_Cur', equal = TRUE)
plot(RCP45_2050_sp.Cur,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Cur',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2050_sp.Cur, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.EffSpatCur.tif"
            ,overwrite=TRUE)

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_ParPlot_EffSpatRugosity.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2050_sp.Rugosity <- spartial(RET_RCP45_2050.model, RCP45_2050_covs, 
                                   x.vars = 'RCP45_2050_Rugosity', equal = TRUE)
plot(RCP45_2050_sp.Rugosity,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Rugosity',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2050_sp.Rugosity, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2050.EffSpatRugosity.tif"
            ,overwrite=TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                      RCP4.5 - 2100 : Covariates                              #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP45_2100_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP45_2100", 
                            pattern="\\.tif$", full.names=T)
RCP45_2100_Cov

#--------------------------- Stack rasters and set crs
RCP45_2100_Cov <- stack(RCP45_2100_Cov)
RCP45_2100_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP45_2100_Cov[[4]],col=cl)
plot(Sapidus,add=T, pch=19,cex=0.3, col="black")

#--------------------------- Crop both dataset manually
#e <- drawExtent() # ATTETION 
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP45_2100_bioc<- crop(RCP45_2100_Cov,e)

#---------------------------  Crop by EEZ
Start <- Sys.time()
masked <- mask(x = RCP45_2100_bioc, mask = EEZ)
plot(masked, col=cl)
RCP45_2100_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP45_2100_bioc, col=cl)
# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 15.2min

#---------------------------  Rename and CRS
RCP45_2100_covs <- RCP45_2100_bioc # covariates
RCP45_2100_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots
plot(stack(RCP45_2100_covs),col=cl) #cropped files

save(RCP45_2100_covs, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Covs.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Covs.rds') 

#---------------------------  Rasterizations and point by grid cell
Sapidus<-as.data.frame(Sapidus)
Sapidus <- SpatialPointsDataFrame(Sapidus[, 2:3], data.frame(RCP45_2100ence = Sapidus[, 1]))
tmp <- rasterize(Sapidus, RCP45_2100_covs[[1]], field = "RCP45_2100ence", fun = "min")
pts.sp1 <- rasterToPoints(tmp,
                          fun = function(x) {
                            x > 0
                          })
nrow(pts.sp1)# 2130

projection(Sapidus) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus) 
st_crs(Sapidus)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract RCP45_2100ence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract RCP45_2100ence values
RCP45_2100.cov <- raster::extract(RCP45_2100_covs, pts.sp1[, 1:2])
RCP45_2100.cov <- na.omit(RCP45_2100.cov)
head(RCP45_2100.cov) #602
nrow(RCP45_2100.cov)
nrow(pts.sp1)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsences                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence <- randomPoints(RCP45_2100_covs, nrow(RCP45_2100.cov))
abs.cov <- raster::extract(RCP45_2100_covs, absence)
plot(absence)
#--------------------------- Code the response
RCP45_2100.cov <- data.frame(RCP45_2100.cov)
RCP45_2100.cov$Sapidus <- 1
abs.cov <- data.frame(abs.cov)
abs.cov$Sapidus <- 0

#--------------------------- And one to bind them
RCP45_2100_all.cov <- rbind(RCP45_2100.cov, abs.cov)
RCP45_2100_all.cov <- RCP45_2100_all.cov[complete.cases(RCP45_2100_all.cov), ]
nrow(RCP45_2100_all.cov)
head(RCP45_2100_all.cov)

RCP45_2100_xvars <- names(RCP45_2100_all.cov)[!(names(RCP45_2100_all.cov) == 'Sapidus')]
RCP45_2100_xvars

# Save files
write.csv(RCP45_2100_all.cov,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP45_2100_all.cov.csv')


save(RCP45_2100_all.cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_all.cov.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_all.cov.rds') 

save(RCP45_2100_xvars, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_xvars.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_xvars.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                          Correlation and collinearity                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RCP45_2100_vif<-vifstep(RCP45_2100_all.cov[,-5])
RCP45_2100_vif

library(GGally)
ggpairs(RCP45_2100_all.cov[,-5])

library(ppcor)
pcor(RCP45_2100_all.cov[,1:3], method = "pearson")

corrplot.mixed(RCP45_2100_all.cov[,-5], lower.col = 'black', number.cex = .7)

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP45_2100_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cchf = Sapidus
#--------------------------- Variable selection and returns the model
RCP45_2100.model <- bart.step(x.data = RCP45_2100_all.cov[, RCP45_2100_xvars],
                              y.data = RCP45_2100_all.cov[, 'Sapidus'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP45_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP45_2100.model<-retune(RCP45_2100.model, reps = 10) 
summary(RET_RCP45_2100.model)
summary(RCP45_2100.model)

save(RET_RCP45_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_retune.model.rds') 

#--------------------------- Spatial prediction
RCP45_2100_Pred.map <- predict(
  object = RET_RCP45_2100.model,
  x.layers = RCP45_2100_covs,
  quantiles = c(0.025, 0.975),
  splitby = 20,
  quiet = TRUE)

save(RCP45_2100_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Pred.map.rds') 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------RCP45_2100: PLOT probs.: mean, 2.5% , 97.5% and uncertainty
#>>>>>>>>>>>>>>>>>>>>>>>>>>>

#~~~~~~~~~~ RCP45_2100_Mean predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_MeanPred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2100_Pred.map[[1]],col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2100: mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
points(pts.sp1[!is.na(raster::extract(RCP45_2100_Pred.map[[1]], 
                                      SpatialPoints(pts.sp1[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ RCP45_2100_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2100_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2100: 2.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()


#~~~~~~~~~~ RCP45_2100_97.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_97.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2100_Pred.map[[3]],col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2100: 97.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ RCP45_2100_Uncertainty predict

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Uncertainty.tiff",
     width=6, height=4, units="in", res=600)

plot(RCP45_2100_Pred.map[[3]]-RCP45_2100_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP45_2100 Uncertainty',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Posterior width', side=2, line=1.3))
dev.off()

#----- RCP45_2100: save rasters
writeRaster(RCP45_2100_Pred.map[[1]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.MeanPred.tif"
            ,overwrite=TRUE)
writeRaster(RCP45_2100_Pred.map[[2]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.2.5Pred.tif"
            ,overwrite=TRUE)
writeRaster(RCP45_2100_Pred.map[[3]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.97.5Pred.tif"
            ,overwrite=TRUE)
RCP45_2100_Uncert<-RCP45_2100_Pred.map[[3]]-RCP45_2100_Pred.map[[2]]# Uncertainty
writeRaster(RCP45_2100_Uncert, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.Uncertainty.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100                           ANALYTICS                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~ Load modified functions
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/varimpEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/summaryEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/partialsEU.R")

#~~~~~~~~~~ Varimp diagonal
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_VarimpDiag.tiff",
     width=6, height=4, units="in", res=600)
varimp.diag(RCP45_2100_all.cov[,RCP45_2100_xvars], RCP45_2100_all.cov[,'Sapidus'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
#varimp(RET_RCP45_2100.model, plots=TRUE)
varimpEU(RET_RCP45_2100.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_RCP45_2100.model)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100              TWO-DIMENSIONAL NICHES                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100e Temperature ~ salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Temp_Sal.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2100.model,plotquants=T,
                xind = c('RCP45_2100_Temp', 'RCP45_2100_Sal'),
                pl = TRUE)
dev.off()

# RCP45_2100e Temperature ~ current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Temp_Cur.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2100.model,plotquants=T,
                xind = c('RCP45_2100_Temp', 'RCP45_2100_Cur'),
                pl = TRUE)
dev.off()

# RCP45_2100e Temperature ~ Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Temp_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP45_2100.model,plotquants=T,
                xind = c('RCP45_2100_Temp', 'RCP45_2100_Rugosity'),
                pl = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100              PARTIAL PLOTS: Binary ~ Predictors                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100e all variables
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_All.tiff",
     width=8, height=5, units="in", res=600)
partial(RET_RCP45_2100.model,
        c('RCP45_2100_Temp','RCP45_2100_Sal','RCP45_2100_Cur','RCP45_2100_Rugosity'),
        trace = FALSE,
        ci = TRUE,
        panel = TRUE,
        smooth = 8)
dev.off()

# RCP45_2100e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_Temp.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2100.model,
          'RCP45_2100_Temp',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP45_2100e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_Sal.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2100.model,
          'RCP45_2100_Sal',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP45_2100e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_Cur.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2100.model,
          'RCP45_2100_Cur',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP45_2100.model,
          'RCP45_2100_Rugosity',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP45_2100e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_EffSpatTemp.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2100_sp.Temp <- spartial(RET_RCP45_2100.model, RCP45_2100_covs, 
                               x.vars = 'RCP45_2100_Temp', equal = TRUE)
plot(RCP45_2100_sp.Temp,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Temp',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2100_sp.Temp, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.EffSpatTemp.tif"
            ,overwrite=TRUE)

# RCP45_2100e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_EffSpatSal.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2100_sp.Sal <- spartial(RET_RCP45_2100.model, RCP45_2100_covs, 
                              x.vars = 'RCP45_2100_Sal', equal = TRUE)
plot(RCP45_2100_sp.Sal,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Sal',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2100_sp.Sal, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.EffSpatSal.tif"
            ,overwrite=TRUE)

# RCP45_2100e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_EffSpatCur.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2100_sp.Cur <- spartial(RET_RCP45_2100.model, RCP45_2100_covs, 
                              x.vars = 'RCP45_2100_Cur', equal = TRUE)
plot(RCP45_2100_sp.Cur,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Cur',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2100_sp.Cur, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.EffSpatCur.tif"
            ,overwrite=TRUE)

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_ParPlot_EffSpatRugosity.tiff",
     width=6, height=4, units="in", res=600)
RCP45_2100_sp.Rugosity <- spartial(RET_RCP45_2100.model, RCP45_2100_covs, 
                                   x.vars = 'RCP45_2100_Rugosity', equal = TRUE)
plot(RCP45_2100_sp.Rugosity,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Rugosity',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP45_2100_sp.Rugosity, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP45_2100.EffSpatRugosity.tif"
            ,overwrite=TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                      RCP8.5 - 2050 : Covariates                              #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP85_2050_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP85_2050", 
                            pattern="\\.tif$", full.names=T)
RCP85_2050_Cov

#--------------------------- Stack rasters and set crs
RCP85_2050_Cov <- stack(RCP85_2050_Cov)
RCP85_2050_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP85_2050_Cov[[4]],col=cl)
plot(Sapidus,add=T, pch=19,cex=0.3, col="black")


#--------------------------- Crop both dataset manually
#e <- drawExtent() # ATTETION 
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP85_2050_bioc<- crop(RCP85_2050_Cov,e)

#---------------------------  Crop by EEZ
Start <- Sys.time()
masked <- mask(x = RCP85_2050_bioc, mask = EEZ)
plot(masked, col=cl)
RCP85_2050_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP85_2050_bioc, col=cl)
# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 15.2min

#---------------------------  Rename and CRS
RCP85_2050_covs <- RCP85_2050_bioc # covariates
RCP85_2050_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots
plot(stack(RCP85_2050_covs),col=cl) #cropped files

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Covs.tiff",
     width=6, height=4, units="in", res=600)

plot(RCP85_2050_covs[[4]], col=cl)
plot(Sapidus,add=T, pch=19,cex=0.3, col="black")

dev.off()

save(RCP85_2050_cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_Covs.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_Covs.rds') 

#---------------------------  Rasterizations and point by grid cell
Sapidus<-as.data.frame(Sapidus)
Sapidus <- SpatialPointsDataFrame(Sapidus[, 2:3], data.frame(RCP85_2050ence = Sapidus[, 1]))
tmp <- rasterize(Sapidus, RCP85_2050_covs[[1]], field = "RCP85_2050ence", fun = "min")
pts.sp1 <- rasterToPoints(tmp,
                          fun = function(x) {
                            x > 0
                          })
nrow(pts.sp1)# 2130

projection(Sapidus) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus) 
st_crs(Sapidus)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract RCP85_2050ence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract RCP85_2050ence values
RCP85_2050.cov <- raster::extract(RCP85_2050_covs, pts.sp1[, 1:2])
RCP85_2050.cov <- na.omit(RCP85_2050.cov)
head(RCP85_2050.cov) #602
nrow(RCP85_2050.cov)
nrow(pts.sp1)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsences                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence <- randomPoints(RCP85_2050_covs, nrow(RCP85_2050.cov))
abs.cov <- raster::extract(RCP85_2050_covs, absence)
plot(absence)
#--------------------------- Code the response
RCP85_2050.cov <- data.frame(RCP85_2050.cov)
RCP85_2050.cov$Sapidus <- 1
abs.cov <- data.frame(abs.cov)
abs.cov$Sapidus <- 0

#--------------------------- And one to bind them
RCP85_2050_all.cov <- rbind(RCP85_2050.cov, abs.cov)
RCP85_2050_all.cov <- RCP85_2050_all.cov[complete.cases(RCP85_2050_all.cov), ]
nrow(RCP85_2050_all.cov)
head(RCP85_2050_all.cov)


RCP85_2050_xvars <- names(RCP85_2050_all.cov)[!(names(RCP85_2050_all.cov) == 'Sapidus')]
RCP85_2050_xvars

# Save files
write.csv(RCP85_2050_all.cov,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP85_2050_all.cov.csv')


save(RCP85_2050_all.cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_all.cov.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_all.cov.rds') 

save(RCP85_2050_xvars, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_xvars.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_xvars.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                          Correlation and collinearity                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
RCP85_2050_vif<-vifstep(RCP85_2050_all.cov[,-5])
RCP85_2050_vif

library(GGally)
ggpairs(RCP85_2050_all.cov[,-5])

library(ppcor)
pcor(RCP85_2050_all.cov[,-5], method = "pearson")

corrplot.mixed(RCP85_2050_all.cov[,-5], lower.col = 'black', number.cex = .7)

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP85_2050_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# cchf = Sapidus
#--------------------------- Variable selection and returns the model
RCP85_2050.model <- bart.step(x.data = RCP85_2050_all.cov[, RCP85_2050_xvars],
                              y.data = RCP85_2050_all.cov[, 'Sapidus'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP85_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP85_2050.model<-retune(RCP85_2050.model, reps = 10) 
summary(RET_RCP85_2050.model)
summary(RCP85_2050.model)

save(RET_RCP85_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_retune.model.rds') 

#--------------------------- Spatial prediction
RCP85_2050_Pred.map <- predict(
  object = RET_RCP85_2050.model,
  x.layers = RCP85_2050_covs,
  quantiles = c(0.025, 0.975),
  splitby = 20,
  quiet = TRUE)

save(RCP85_2050_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_Pred.map.rds') 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------RCP85_2050: PLOT probs.: mean, 2.5% , 97.5% and uncertainty
#>>>>>>>>>>>>>>>>>>>>>>>>>>>

#~~~~~~~~~~ RCP85_2050_Mean predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_MeanPred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2050_Pred.map[[1]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2050: mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
points(pts.sp1[!is.na(raster::extract(RCP85_2050_Pred.map[[1]], 
                                      SpatialPoints(pts.sp1[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ RCP85_2050_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2050_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2050: 2.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()


#~~~~~~~~~~ RCP85_2050_97.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_97.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2050_Pred.map[[3]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2050: 97.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ RCP85_2050_Uncertainty predict

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Uncertainty.tiff",
     width=6, height=4, units="in", res=600)

plot(RCP85_2050_Pred.map[[3]]-RCP85_2050_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2050 Uncertainty',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Posterior width', side=2, line=1.3))
dev.off()

#----- RCP85_2050: save rasters
writeRaster(RCP85_2050_Pred.map[[1]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.MeanPred.tif"
            ,overwrite=TRUE)
writeRaster(RCP85_2050_Pred.map[[2]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.2.5Pred.tif"
            ,overwrite=TRUE)
writeRaster(RCP85_2050_Pred.map[[3]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.97.5Pred.tif"
            ,overwrite=TRUE)
RCP85_2050_Uncert<-RCP85_2050_Pred.map[[3]]-RCP85_2050_Pred.map[[2]]# Uncertainty
writeRaster(RCP85_2050_Uncert, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.Uncertainty.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050                           ANALYTICS                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~ Load modified functions
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/varimpEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/summaryEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/partialsEU.R")

#~~~~~~~~~~ Varimp diagonal
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_VarimpDiag.tiff",
     width=6, height=4, units="in", res=600)
varimp.diag(RCP85_2050_all.cov[,RCP85_2050_xvars], RCP85_2050_all.cov[,'Sapidus'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
#varimp(RET_RCP85_2050.model, plots=TRUE)
varimpEU(RET_RCP85_2050.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_RCP85_2050.model)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050              TWO-DIMENSIONAL NICHES                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050e Temperature ~ salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Temp_Sal.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2050.model,plotquants=T,
                xind = c('RCP85_2050_Temp', 'RCP85_2050_Sal'),
                pl = TRUE)
dev.off()

# RCP85_2050e Temperature ~ current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Temp_Cur.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2050.model,plotquants=T,
                xind = c('RCP85_2050_Temp', 'RCP85_2050_Cur'),
                pl = TRUE)
dev.off()

# RCP85_2050e Temperature ~ Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Temp_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2050.model,plotquants=T,
                xind = c('RCP85_2050_Temp', 'RCP85_2050_Rugosity'),
                pl = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050          PARTIAL PLOTS: Binary ~ Predictors                       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050e all variables
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_All.tiff",
     width=8, height=5, units="in", res=600)
partial(RET_RCP85_2050.model,
        c('RCP85_2050_Temp','RCP85_2050_Sal','RCP85_2050_Cur','RCP85_2050_Rugosity'),
        trace = FALSE,
        ci = TRUE,
        panel = TRUE,
        smooth = 8)
dev.off()

# RCP85_2050e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_Temp.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2050.model,
          'RCP85_2050_Temp',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP85_2050e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_Sal.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2050.model,
          'RCP85_2050_Sal',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP85_2050e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_Cur.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2050.model,
          'RCP85_2050_Cur',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2050.model,
          'RCP85_2050_Rugosity',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP85_2050e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_EffSpatTemp.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2050_sp.Temp <- spartial(RET_RCP85_2050.model, RCP85_2050_covs, 
                               x.vars = 'RCP85_2050_Temp', equal = TRUE)
plot(RCP85_2050_sp.Temp,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Temp',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2050_sp.Temp, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.EffSpatTemp.tif"
            ,overwrite=TRUE)

# RCP85_2050e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_EffSpatSal.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2050_sp.Sal <- spartial(RET_RCP85_2050.model, RCP85_2050_covs, 
                              x.vars = 'RCP85_2050_Sal', equal = TRUE)
plot(RCP85_2050_sp.Sal,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Sal',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2050_sp.Sal, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.EffSpatSal.tif"
            ,overwrite=TRUE)

# RCP85_2050e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_EffSpatCur.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2050_sp.Cur <- spartial(RET_RCP85_2050.model, RCP85_2050_covs, 
                              x.vars = 'RCP85_2050_Cur', equal = TRUE)
plot(RCP85_2050_sp.Cur,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Cur',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2050_sp.Cur, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.EffSpatCur.tif"
            ,overwrite=TRUE)


# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_ParPlot_EffSpatRugosity.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2050_sp.Rugosity <- spartial(RET_RCP85_2050.model, RCP85_2050_covs, 
                                   x.vars = 'RCP85_2050_Rugosity', equal = TRUE)
plot(RCP85_2050_sp.Rugosity,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Rugosity',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2050_sp.Rugosity, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2050.EffSpatRugosity.tif"
            ,overwrite=TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                      RCP8.5 - 2100 : Covariates                              #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#RCP85_2100_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP85_2100", 
                            pattern="\\.tif$", full.names=T)
#RCP85_2100_Cov

RCP85_2100.model <- bart.step(x.data = RCP85_2100_all.cov[, RCP85_2100_xvars],
                              y.data = RCP85_2100_all.cov[, 'Sapidus'],
                              full = TRUE,
                              quiet = TRUE)

#--------------------------- Retune: cross-validation
RET_RCP85_2100.model<-retune(RCP85_2100.model, reps = 10) 
summary(RET_RCP85_2100.model)

save(RET_RCP85_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_retune.model.rds') 

#--------------------------- Spatial prediction
RCP85_2100_Pred.map <- predict(
                          object = RCP85_2100.model,
                          x.layers = RCP85_2100_covs,
                          quantiles = c(0.025, 0.975),
                          splitby = 20,
                          quiet = TRUE)

save(RCP85_2100_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_Pred.map.rds') 

#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#-------------RCP85_2100: PLOT probs.: mean, 2.5% , 97.5% and uncertainty
#>>>>>>>>>>>>>>>>>>>>>>>>>>>

#~~~~~~~~~~ RCP85_2100_Mean predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_MeanPred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2100_Pred.map[[1]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2100: mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
points(pts.sp1[!is.na(raster::extract(RCP85_2100_Pred.map[[1]], 
                                      SpatialPoints(pts.sp1[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ RCP85_2100_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2100_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2100: 2.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()


#~~~~~~~~~~ RCP85_2100_97.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_97.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP85_2100_Pred.map[[3]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2100: 97.5% mean predict',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Probability', side=2, line=1.3))
dev.off()

#~~~~~~~~~~ RCP85_2100_Uncertainty predict

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Uncertainty.tiff",
     width=6, height=4, units="in", res=600)

plot(RCP85_2100_Pred.map[[3]]-RCP85_2100_Pred.map[[2]],col=cl,
     box = F,
     axes = F,
     main = 'RCP85_2100 Uncertainty',
     zlim = c(0,1),
     axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Posterior width', side=2, line=1.3))
dev.off()

#----- RCP85_2100: save rasters
writeRaster(RCP85_2100_Pred.map[[1]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.MeanPred.tif"
            ,overwrite=TRUE)
writeRaster(RCP85_2100_Pred.map[[2]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.2.5Pred.tif"
            ,overwrite=TRUE)
writeRaster(RCP85_2100_Pred.map[[3]], 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.97.5Pred.tif"
            ,overwrite=TRUE)
RCP85_2100_Uncert<-RCP85_2100_Pred.map[[3]]-RCP85_2100_Pred.map[[2]]# Uncertainty
writeRaster(RCP85_2100_Uncert, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.Uncertainty.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100                           ANALYTICS                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~ Load modified functions
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/varimpEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/summaryEU.R")
source("~/efc/R_Stat/2022_NEMA/SDM_Bayes/partialsEU.R")

#~~~~~~~~~~ Varimp diagonal
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_VarimpDiag.tiff",
     width=6, height=4, units="in", res=600)
varimp.diag(RCP85_2100_all.cov[,RCP85_2100_xvars], RCP85_2100_all.cov[,'Sapidus'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
varimpEU(RET_RCP85_2100.model, plots=TRUE)
#varimpEU(RET_RCP85_2100.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_RCP85_2100.model)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100              TWO-DIMENSIONAL NICHES                                  #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100e Temperature ~ salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Temp_Sal.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2100.model,plotquants=T,
                xind = c('RCP85_2100_Temp', 'RCP85_2100_Sal'),
                pl = TRUE)
dev.off()

# RCP85_2100e Temperature ~ current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Temp_Cur.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2100.model,plotquants=T,
                xind = c('RCP85_2100_Temp', 'RCP85_2100_Cur'),
                pl = TRUE)
dev.off()

# RCP85_2100e Temperature ~ Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Temp_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
dbarts::pd2bart(RET_RCP85_2100.model,plotquants=T,
                xind = c('RCP85_2100_Temp', 'RCP85_2100_Rugosity'),
                pl = TRUE)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100              PARTIAL PLOTS: Binary ~ PRedictors                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100e all variables
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_All.tiff",
     width=8, height=5, units="in", res=600)
partial(RET_RCP85_2100.model,
        c('RCP85_2100_Temp','RCP85_2100_Sal','RCP85_2100_Cur','RCP85_2100_Rugosity'),
        trace = FALSE,
        ci = TRUE,
        panel = TRUE,
        smooth = 8)
dev.off()

# RCP85_2100e Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_Temp.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2100.model,
          'RCP85_2100_Temp',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP85_2100e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_Sal.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2100.model,
          'RCP85_2100_Sal',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# RCP85_2100e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_Cur.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2100.model,
          'RCP85_2100_Cur',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_Rugosity.tiff",
     width=6, height=4, units="in", res=600)
partialEU(RET_RCP85_2100.model,
          'RCP85_2100_Rugosity',
          trace = FALSE,
          ci = TRUE,
          panel = FALSE,
          smooth = 8)
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP85_2100 Temperature
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_EffSpatTemp.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2100_sp.Temp <- spartial(RET_RCP85_2100.model, RCP85_2100_covs, 
                               x.vars = 'RCP85_2100_Temp', equal = TRUE)
plot(RCP85_2100_sp.Temp,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Temp',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2100_sp.Temp, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.EffSpatTemp.tif"
            ,overwrite=TRUE)

# RCP85_2100e Salinity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_EffSpatSal.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2100_sp.Sal <- spartial(RET_RCP85_2100.model, RCP85_2100_covs, 
                              x.vars = 'RCP85_2100_Sal', equal = TRUE)
plot(RCP85_2100_sp.Sal,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Sal',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2100_sp.Sal, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.EffSpatSal.tif"
            ,overwrite=TRUE)

# RCP85_2100e Current
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_EffSpatCur.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2100_sp.Cur <- spartial(RET_RCP85_2100.model, RCP85_2100_covs, 
                              x.vars = 'RCP85_2100_Cur', equal = TRUE)
plot(RCP85_2100_sp.Cur,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Cur',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2100_sp.Cur, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.EffSpatCur.tif"
            ,overwrite=TRUE)

# Rugosity
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_ParPlot_EffSpatRugosity.tiff",
     width=6, height=4, units="in", res=600)
RCP85_2100_sp.Rugosity <- spartial(RET_RCP85_2100.model, RCP85_2100_covs, 
                                   x.vars = 'RCP85_2100_Rugosity', equal = TRUE)
plot(RCP85_2100_sp.Rugosity,col=cl, 
     box = FALSE,
     axes = FALSE,
     main = 'Spartial plot: Rugosity',
     #zlim = c(0,1),
     #axis.args=list(at=pretty(0:1), labels=pretty(0:1)),
     legend.args=list(text='Partial effect', side=2, line=1.3))
dev.off()
writeRaster(RCP85_2100_sp.Rugosity, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCP85_2100.EffSpatRugosity.tif"
            ,overwrite=TRUE)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                              RCPS_ANALYSIS                                   #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~ RCP_4.5:2050
mfrow=c(2,3)
Current<-Present_Pred.map[[1]]> 0.7
RCP45_2050<-RCP45_2050_Pred.map[[1]]> 0.7
RCP45_2100<-RCP45_2100_Pred.map[[1]]> 0.7
RCP85_2050<-RCP85_2050_Pred.map[[1]]> 0.7
RCP85_2100<-RCP85_2100_Pred.map[[1]]> 0.7

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_plots/RCP_4.5_2050.tiff",
     width=6, height=4, units="in", res=600)

RCP4.5_2050 <- RCP45_2050-Current
plot(RCP4.5_2050)

dev.off()
writeRaster(RCP4.5_2050, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCPs_Differences/RCP4.5_2050_70.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~ RCP_4.5:2100
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_plots/RCP_4.5_2100.tiff",
     width=6, height=4, units="in", res=600)

RCP4.5_2100 <- RCP45_2100-Current
plot(RCP4.5_2100)

dev.off()
writeRaster(RCP4.5_2100, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCPs_Differences/RCP4.5_2100_70.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~ RCP_8.5:2050
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_plots/RCP_8.5_2050.tiff",
     width=6, height=4, units="in", res=600)

RCP8.5_2050 <- RCP85_2050-Current
plot(RCP8.5_2050)

dev.off()
#Test
RCP8.5_2050_Uncert<-RCP8.5_2050-RCP85_2100_Uncert
plot(RCP8.5_2050_Uncert)
writeRaster(RCP8.5_2050, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCPs_Differences/RCP8.5_2050_70.tif"
            ,overwrite=TRUE)

#~~~~~~~~~~~~~~~~ RCP_8.5:2100
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_plots/RCP_4.5_2100.tiff",
     width=6, height=4, units="in", res=600)

RCP8.5_2100 <- RCP85_2100-Current
plot(RCP4.5_2100)

dev.off()
writeRaster(RCP8.5_2100, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCPs_Differences/RCP8.5_2100_70.tif"
            ,overwrite=TRUE)



# create classification matrix
m=c(-0.6, 0, 1,0, .01, 2,0.1, Inf, 3)
mat=matrix(m,ncol=3,byrow = T)
elev=reclassify(RCP8.5_2100, mat)
plot(elev)
writeRaster(elev, 
            filename = "~/efc/R_Stat/2022_NEMA/SDM_Bayes/4_Rasters/RCPs_Differences/elev.tif"
            ,overwrite=TRUE)


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#               Plot: Prediction ~  lat Long Native and nonnative              #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/efc/R_Stat/2022_NEMA/SDM_Bayes/8_PredAnalysis")
MPred<-read.csv("DataMeanPred.csv",dec=".",sep=";",header=T)
str(MPred)


ggplot(MPred, aes(x=Scenary, y=Mean_Pred, fill=Habitat)) +
  geom_boxplot(outlier.shape = NA) +
  coord_cartesian(ylim = quantile(MPred$Mean_Pred, c(0, 1)))+ 
  facet_wrap(~Habitat)




library(mgcv)
require(gam)
mod_gam1 = gam(Mean_Pred ~ s(Lat, bs="cs")+Habitat, data = MPred)
summary(mod_gam1)
plot(mod_gam1,se = TRUE)


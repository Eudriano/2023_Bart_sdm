#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Author: Eudriano Costa
# Manuscrispt: Aquatic species shows asymmetric distribution range shifts in native and
# non-native areas 
# SDM - Species Distribution Models  
# R version 4.0.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())

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
#
# Scenario model:    Current/present (2000 - 2014)                             #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables

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

spg1<-as.data.frame(spg1)
Sapidus_Pres<-spg1
Sapidus_Pres<-Sapidus_Pres[, c('Longitude', 'Latitude')]
Sapidus_Pres$Presence = 1
nrow(Sapidus_Pres) #27652

write.csv(Sapidus_Pres,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_Pres_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

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
#                      PRESENT- 1830 to 2022: Covariates                       #
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
Time.masked # 0.1073102 mins

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
Sapidus_Pres <- SpatialPointsDataFrame(Sapidus_Pres[, 1:2], data.frame(Presence = Sapidus_Pres[, 3]))
tmp_Pres <- rasterize(Sapidus_Pres, Present_covs[[1]], field = "Presence", fun = "min")
pts.sp1_Pres <- rasterToPoints(tmp_Pres,
                               fun = function(x) {
                                 x > 0
                               })
nrow(pts.sp1_Pres)# 2121

projection(Sapidus_Pres) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus_Pres) 
st_crs(Sapidus_Pres)

write.csv(pts.sp1_Pres,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Presence_Rasterization.csv')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
pres.cov <- raster::extract(Present_covs, pts.sp1_Pres[, 1:2])
pres.cov <- na.omit(pres.cov)
head(pres.cov) #602, 668
nrow(pres.cov)
nrow(pts.sp1_Pres) #2121
which(is.na(pres.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsence_Press                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence_Pres <- randomPoints(Present_covs, nrow(pres.cov))

plot(absence_Pres)

nrow(absence_Pres)

abs.cov_Pres <- raster::extract(Present_covs, absence_Pres) #668

write.csv(absence_Pres,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')

#--------------------------- Code the response
pres.cov <- data.frame(pres.cov)
pres.cov$Sapidus_Pres <- 1
abs.cov_Pres <- data.frame(abs.cov_Pres)
abs.cov_Pres$Sapidus_Pres <- 0

#--------------------------- And one to bind them
Present_all.cov <- rbind(pres.cov, abs.cov_Pres)
Present_all.cov <- Present_all.cov[complete.cases(Present_all.cov), ]
nrow(Present_all.cov)
head(Present_all.cov)
which(is.na(Present_all.cov))

Present_xvars <- names(Present_all.cov)[!(names(Present_all.cov) == 'Sapidus_Pres')]
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
# Sapidus_Pres
#--------------------------- Variable selection and returns the model
Present.model <- bart.step(x.data = Present_all.cov[, Present_xvars],
                           y.data = Present_all.cov[, 'Sapidus_Pres'],
                           full = TRUE,
                           quiet = TRUE)

save(Present.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/Present_model.rds') 

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
points(pts.sp1_Pres[!is.na(raster::extract(Present_Pred.map[[1]], 
                                           SpatialPoints(pts.sp1_Pres[,1:2]))),],
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
varimp.diag(Present_all.cov[,Present_xvars], Present_all.cov[,'Sapidus_Pres'], 
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Present              TWO-DIMENSIONAL "NICHE"                                  #
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

#
#----------------------------------- / / --------------------------------------#
#
rm(list = ls())
Start <- Sys.time()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#  Scenario model:         RCP45_2050                                          #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
names(spg1@data) <- 'RCP45_2050'
nrow(spg1)
head(spg1)

spg1<-as.data.frame(spg1)
Sapidus_RCP45_2050<-spg1
Sapidus_RCP45_2050<-Sapidus_RCP45_2050[, c('Longitude', 'Latitude')]
Sapidus_RCP45_2050$RCP45_2050 = 1
nrow(Sapidus_RCP45_2050) #27652

write.csv(Sapidus_RCP45_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_RCP45_2050_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

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
#                       RCP45_2050: Covariates                                 #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP45_2050_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP45_2050", 
                            pattern="\\.tif$", full.names=T)
RCP45_2050_Cov

#--------------------------- Stack rasters and set crs
RCP45_2050_Cov <- stack(RCP45_2050_Cov)
RCP45_2050_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP45_2050_Cov[[3]],col=cl)

#--------------------------- Crop both dataset manually
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP45_2050_bioc<- crop(RCP45_2050_Cov,e)

#---------------------------  Crop by EEZ
masked <- mask(x = RCP45_2050_bioc, mask = EEZ)
plot(masked, col=cl)
RCP45_2050_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP45_2050_bioc, col=cl)

#---------------------------  Rename and CRS
RCP45_2050_covs <- RCP45_2050_bioc # covariates
RCP45_2050_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_Covs.tiff",
     width=8, height=4, units="in", res=600)

plot(RCP45_2050_covs, col=cl)

dev.off()

#---------------------------  Rasterizations and point by grid cell
Sapidus_RCP45_2050 <- SpatialPointsDataFrame(Sapidus_RCP45_2050[, 1:2], data.frame(RCP45_2050 = Sapidus_RCP45_2050[, 3]))
tmp_RCP45_2050 <- rasterize(Sapidus_RCP45_2050, RCP45_2050_covs[[1]], field = "RCP45_2050", fun = "min")
pts.sp1_RCP45_2050 <- rasterToPoints(tmp_RCP45_2050,
                                     fun = function(x) {
                                       x > 0
                                     })
nrow(pts.sp1_RCP45_2050)# 2121

projection(Sapidus_RCP45_2050) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus_RCP45_2050) 
st_crs(Sapidus_RCP45_2050)

write.csv(pts.sp1_RCP45_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP45_2050_Rasterization.csv')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
RCP45_2050.cov <- raster::extract(RCP45_2050_covs, pts.sp1_RCP45_2050[, 1:2])
RCP45_2050.cov <- na.omit(RCP45_2050.cov)
head(RCP45_2050.cov) #602, 668
nrow(RCP45_2050.cov)
nrow(pts.sp1_RCP45_2050) #2121
which(is.na(RCP45_2050.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsence_RCP45_2050s                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence_RCP45_2050 <- randomPoints(RCP45_2050_covs, nrow(RCP45_2050.cov))

plot(absence_RCP45_2050)

nrow(absence_RCP45_2050)

abs.cov_RCP45_2050 <- raster::extract(RCP45_2050_covs, absence_RCP45_2050) #668

write.csv(absence_RCP45_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')

#--------------------------- Code the response
RCP45_2050.cov <- data.frame(RCP45_2050.cov)
RCP45_2050.cov$Sapidus_RCP45_2050 <- 1
abs.cov_RCP45_2050 <- data.frame(abs.cov_RCP45_2050)
abs.cov_RCP45_2050$Sapidus_RCP45_2050 <- 0

#--------------------------- And one to bind them
RCP45_2050_all.cov <- rbind(RCP45_2050.cov, abs.cov_RCP45_2050)
RCP45_2050_all.cov <- RCP45_2050_all.cov[complete.cases(RCP45_2050_all.cov), ]
nrow(RCP45_2050_all.cov)
head(RCP45_2050_all.cov)
which(is.na(RCP45_2050_all.cov))

RCP45_2050_xvars <- names(RCP45_2050_all.cov)[!(names(RCP45_2050_all.cov) == 'Sapidus_RCP45_2050')]
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

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP45_2050_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2050                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#--------------------------- Variable selection and returns the model
RCP45_2050.model <- bart.step(x.data = RCP45_2050_all.cov[, RCP45_2050_xvars],
                              y.data = RCP45_2050_all.cov[, 'Sapidus_RCP45_2050'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP45_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2050_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP45_2050.model<-retune(RCP45_2050.model, reps = 10) 
summary(RET_RCP45_2050.model)

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
points(pts.sp1_RCP45_2050[!is.na(raster::extract(RCP45_2050_Pred.map[[1]], 
                                                 SpatialPoints(pts.sp1_RCP45_2050[,1:2]))),],
       col = 'black',
       pch = 16,
       cex = 0.3)
dev.off()

#~~~~~~~~~~ RCP45_2050_2.5% predict
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2050_2.5_Pred.tiff",
     width=6, height=4, units="in", res=600)
plot(RCP45_2050_Pred.map[[2]],col=cl,
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
plot(RCP45_2050_Pred.map[[3]],col=cl,
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
plot(RCP45_2050_Pred.map[[3]]-RCP45_2050_Pred.map[[2]],col=cl,
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
varimp.diag(RCP45_2050_all.cov[,RCP45_2050_xvars], RCP45_2050_all.cov[,'Sapidus_RCP45_2050'], 
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
# RCP45_2050       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP45_2050 Temperature
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

# RCP45_2050 Salinity
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

# RCP45_2050 Current
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

# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 0.1073102 mins

#
#----------------------------------- / / --------------------------------------#
#
rm(list = ls())
Start <- Sys.time()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#  Scenario model:         RCP45_2100                                          #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
names(spg1@data) <- 'RCP45_2100'
nrow(spg1)
head(spg1)

spg1<-as.data.frame(spg1)
Sapidus_RCP45_2100<-spg1
Sapidus_RCP45_2100<-Sapidus_RCP45_2100[, c('Longitude', 'Latitude')]
Sapidus_RCP45_2100$RCP45_2100 = 1
nrow(Sapidus_RCP45_2100) #27652

write.csv(Sapidus_RCP45_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_RCP45_2100_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

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
#                       RCP45_2100: Covariates                                 #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP45_2100_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP45_2100", 
                            pattern="\\.tif$", full.names=T)
RCP45_2100_Cov

#--------------------------- Stack rasters and set crs
RCP45_2100_Cov <- stack(RCP45_2100_Cov)
RCP45_2100_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP45_2100_Cov[[3]],col=cl)

#--------------------------- Crop both dataset manually
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP45_2100_bioc<- crop(RCP45_2100_Cov,e)

#---------------------------  Crop by EEZ
masked <- mask(x = RCP45_2100_bioc, mask = EEZ)
plot(masked, col=cl)
RCP45_2100_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP45_2100_bioc, col=cl)

#---------------------------  Rename and CRS
RCP45_2100_covs <- RCP45_2100_bioc # covariates
RCP45_2100_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_Covs.tiff",
     width=8, height=4, units="in", res=600)

plot(RCP45_2100_covs, col=cl)

dev.off()

#---------------------------  Rasterizations and point by grid cell
Sapidus_RCP45_2100 <- SpatialPointsDataFrame(Sapidus_RCP45_2100[, 1:2], data.frame(RCP45_2100 = Sapidus_RCP45_2100[, 3]))
tmp_RCP45_2100 <- rasterize(Sapidus_RCP45_2100, RCP45_2100_covs[[1]], field = "RCP45_2100", fun = "min")
pts.sp1_RCP45_2100 <- rasterToPoints(tmp_RCP45_2100,
                                     fun = function(x) {
                                       x > 0
                                     })
nrow(pts.sp1_RCP45_2100)# 2121

projection(Sapidus_RCP45_2100) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus_RCP45_2100) 
st_crs(Sapidus_RCP45_2100)

write.csv(pts.sp1_RCP45_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP45_2100_Rasterization.csv')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
RCP45_2100.cov <- raster::extract(RCP45_2100_covs, pts.sp1_RCP45_2100[, 1:2])
RCP45_2100.cov <- na.omit(RCP45_2100.cov)
head(RCP45_2100.cov) #602, 668
nrow(RCP45_2100.cov)
nrow(pts.sp1_RCP45_2100) #2121
which(is.na(RCP45_2100.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsence_RCP45_2100s                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence_RCP45_2100 <- randomPoints(RCP45_2100_covs, nrow(RCP45_2100.cov))

plot(absence_RCP45_2100)

nrow(absence_RCP45_2100)

abs.cov_RCP45_2100 <- raster::extract(RCP45_2100_covs, absence_RCP45_2100) #668

write.csv(absence_RCP45_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')

#--------------------------- Code the response
RCP45_2100.cov <- data.frame(RCP45_2100.cov)
RCP45_2100.cov$Sapidus_RCP45_2100 <- 1
abs.cov_RCP45_2100 <- data.frame(abs.cov_RCP45_2100)
abs.cov_RCP45_2100$Sapidus_RCP45_2100 <- 0

#--------------------------- And one to bind them
RCP45_2100_all.cov <- rbind(RCP45_2100.cov, abs.cov_RCP45_2100)
RCP45_2100_all.cov <- RCP45_2100_all.cov[complete.cases(RCP45_2100_all.cov), ]
nrow(RCP45_2100_all.cov)
head(RCP45_2100_all.cov)
which(is.na(RCP45_2100_all.cov))

RCP45_2100_xvars <- names(RCP45_2100_all.cov)[!(names(RCP45_2100_all.cov) == 'Sapidus_RCP45_2100')]
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

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP45_2100_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP45_2100_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP45_2100                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#--------------------------- Variable selection and returns the model
RCP45_2100.model <- bart.step(x.data = RCP45_2100_all.cov[, RCP45_2100_xvars],
                              y.data = RCP45_2100_all.cov[, 'Sapidus_RCP45_2100'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP45_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP45_2100.model<-retune(RCP45_2100.model, reps = 10) 
summary(RET_RCP45_2100.model)

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

plot(RCP45_2100_Pred.map)

save(RCP45_2100_Pred.map,
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Pred.map.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP45_2100_Pred.map.rds') 

l#>>>>>>>>>>>>>>>>>>>>>>>>>>>
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
points(pts.sp1_RCP45_2100[!is.na(raster::extract(RCP45_2100_Pred.map[[1]], 
                                                 SpatialPoints(pts.sp1_RCP45_2100[,1:2]))),],
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
varimp.diag(RCP45_2100_all.cov[,RCP45_2100_xvars], RCP45_2100_all.cov[,'Sapidus_RCP45_2100'], 
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
# RCP45_2100       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP45_2100 Temperature
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

# RCP45_2100 Salinity
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

# RCP45_2100 Current
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

# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 0.1073102 mins

#
#----------------------------------- / / --------------------------------------#
#
rm(list = ls())
Start <- Sys.time()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#  Scenario model:         RCP85_2050                                          #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
names(spg1@data) <- 'RCP85_2050'
nrow(spg1)
head(spg1)

spg1<-as.data.frame(spg1)
Sapidus_RCP85_2050<-spg1
Sapidus_RCP85_2050<-Sapidus_RCP85_2050[, c('Longitude', 'Latitude')]
Sapidus_RCP85_2050$RCP85_2050 = 1
nrow(Sapidus_RCP85_2050) #27652

write.csv(Sapidus_RCP85_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_RCP85_2050_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

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
#                       RCP85_2050: Covariates                                 #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP85_2050_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP85_2050", 
                            pattern="\\.tif$", full.names=T)
RCP85_2050_Cov

#--------------------------- Stack rasters and set crs
RCP85_2050_Cov <- stack(RCP85_2050_Cov)
RCP85_2050_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP85_2050_Cov[[3]],col=cl)

#--------------------------- Crop both dataset manually
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP85_2050_bioc<- crop(RCP85_2050_Cov,e)

#---------------------------  Crop by EEZ
masked <- mask(x = RCP85_2050_bioc, mask = EEZ)
plot(masked, col=cl)
RCP85_2050_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP85_2050_bioc, col=cl)

#---------------------------  Rename and CRS
RCP85_2050_covs <- RCP85_2050_bioc # covariates
RCP85_2050_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_Covs.tiff",
     width=8, height=4, units="in", res=600)

plot(RCP85_2050_covs, col=cl)

dev.off()

#---------------------------  Rasterizations and point by grid cell
Sapidus_RCP85_2050 <- SpatialPointsDataFrame(Sapidus_RCP85_2050[, 1:2], data.frame(RCP85_2050 = Sapidus_RCP85_2050[, 3]))
tmp_RCP85_2050 <- rasterize(Sapidus_RCP85_2050, RCP85_2050_covs[[1]], field = "RCP85_2050", fun = "min")
pts.sp1_RCP85_2050 <- rasterToPoints(tmp_RCP85_2050,
                                     fun = function(x) {
                                       x > 0
                                     })
nrow(pts.sp1_RCP85_2050)# 2121

projection(Sapidus_RCP85_2050) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus_RCP85_2050) 
st_crs(Sapidus_RCP85_2050)

write.csv(pts.sp1_RCP85_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP85_2050_Rasterization.csv')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
RCP85_2050.cov <- raster::extract(RCP85_2050_covs, pts.sp1_RCP85_2050[, 1:2])
RCP85_2050.cov <- na.omit(RCP85_2050.cov)
head(RCP85_2050.cov) #602, 668
nrow(RCP85_2050.cov)
nrow(pts.sp1_RCP85_2050) #2121
which(is.na(RCP85_2050.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsence_RCP85_2050s                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence_RCP85_2050 <- randomPoints(RCP85_2050_covs, nrow(RCP85_2050.cov))

plot(absence_RCP85_2050)

nrow(absence_RCP85_2050)

abs.cov_RCP85_2050 <- raster::extract(RCP85_2050_covs, absence_RCP85_2050) #668

write.csv(absence_RCP85_2050,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')

#--------------------------- Code the response
RCP85_2050.cov <- data.frame(RCP85_2050.cov)
RCP85_2050.cov$Sapidus_RCP85_2050 <- 1
abs.cov_RCP85_2050 <- data.frame(abs.cov_RCP85_2050)
abs.cov_RCP85_2050$Sapidus_RCP85_2050 <- 0

#--------------------------- And one to bind them
RCP85_2050_all.cov <- rbind(RCP85_2050.cov, abs.cov_RCP85_2050)
RCP85_2050_all.cov <- RCP85_2050_all.cov[complete.cases(RCP85_2050_all.cov), ]
nrow(RCP85_2050_all.cov)
head(RCP85_2050_all.cov)
which(is.na(RCP85_2050_all.cov))

RCP85_2050_xvars <- names(RCP85_2050_all.cov)[!(names(RCP85_2050_all.cov) == 'Sapidus_RCP85_2050')]
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

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2050_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP85_2050_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2050                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#--------------------------- Variable selection and returns the model
RCP85_2050.model <- bart.step(x.data = RCP85_2050_all.cov[, RCP85_2050_xvars],
                              y.data = RCP85_2050_all.cov[, 'Sapidus_RCP85_2050'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP85_2050.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2050_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP85_2050.model<-retune(RCP85_2050.model, reps = 10) 
summary(RET_RCP85_2050.model)

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
points(pts.sp1_RCP85_2050[!is.na(raster::extract(RCP85_2050_Pred.map[[1]], 
                                                 SpatialPoints(pts.sp1_RCP85_2050[,1:2]))),],
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
varimp.diag(RCP85_2050_all.cov[,RCP85_2050_xvars], RCP85_2050_all.cov[,'Sapidus_RCP85_2050'], 
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
# RCP85_2050       SPARTIAL EFFECT - PARTIAL MAPS                     #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# RCP85_2050 Temperature
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

# RCP85_2050 Salinity
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

# RCP85_2050 Current
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

# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 0.1073102 mins

#
#----------------------------------- / / --------------------------------------#
#
rm(list = ls())
Start <- Sys.time()
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#  Scenario model:         RCP85_2100                                          #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
names(spg1@data) <- 'RCP85_2100'
nrow(spg1)
head(spg1)

spg1<-as.data.frame(spg1)
Sapidus_RCP85_2100<-spg1
Sapidus_RCP85_2100<-Sapidus_RCP85_2100[, c('Longitude', 'Latitude')]
Sapidus_RCP85_2100$RCP85_2100 = 1
nrow(Sapidus_RCP85_2100) #27652

write.csv(Sapidus_RCP85_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Sapidus_RCP85_2100_1.csv')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               EEZ to crop data                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
EEZ <-read_sf(
  '~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/EEZ_shapes/EEZ_0_200.shp'#EEZ_NoISla
)
plot(EEZ,col='aliceblue')

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
#                       RCP85_2100: Covariates                                 #
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
RCP85_2100_Cov <-list.files(path="~/efc/R_Stat/2022_NEMA/SDM_Bayes/2_Predictors/RCP85_2100", 
                            pattern="\\.tif$", full.names=T)
RCP85_2100_Cov

#--------------------------- Stack rasters and set crs
RCP85_2100_Cov <- stack(RCP85_2100_Cov)
RCP85_2100_Cov@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

plot(RCP85_2100_Cov[[3]],col=cl)

#--------------------------- Crop both dataset manually
# --- crop by LAT and LONG
e <- extent(c(-99.89664,45.00737, -47.29764, 62.02496))
RCP85_2100_bioc<- crop(RCP85_2100_Cov,e)

#---------------------------  Crop by EEZ
masked <- mask(x = RCP85_2100_bioc, mask = EEZ)
plot(masked, col=cl)
RCP85_2100_bioc <- crop(x = masked, y = extent(EEZ))
plot(RCP85_2100_bioc, col=cl)

#---------------------------  Rename and CRS
RCP85_2100_covs <- RCP85_2100_bioc # covariates
RCP85_2100_covs@crs <- crs('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0') 

#---------------------------  Covariate Plots

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Covs.tiff",
     width=8, height=4, units="in", res=600)

plot(RCP85_2100_covs, col=cl)

dev.off()

#---------------------------  Rasterizations and point by grid cell
Sapidus_RCP85_2100 <- SpatialPointsDataFrame(Sapidus_RCP85_2100[, 1:2], data.frame(RCP85_2100 = Sapidus_RCP85_2100[, 3]))
tmp_RCP85_2100 <- rasterize(Sapidus_RCP85_2100, RCP85_2100_covs[[1]], field = "RCP85_2100", fun = "min")
pts.sp1_RCP85_2100 <- rasterToPoints(tmp_RCP85_2100,
                                     fun = function(x) {
                                       x > 0
                                     })
nrow(pts.sp1_RCP85_2100)# 2121

projection(Sapidus_RCP85_2100) <- CRS('+proj=longlat +ellps=WGS84 
	               +datum=WGS84 +no_defs +towgs84=0,0,0')
class(Sapidus_RCP85_2100) 
st_crs(Sapidus_RCP85_2100)

write.csv(pts.sp1_RCP85_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP85_2100_Rasterization.csv')

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Extract presence values                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract presence values
RCP85_2100.cov <- raster::extract(RCP85_2100_covs, pts.sp1_RCP85_2100[, 1:2])
RCP85_2100.cov <- na.omit(RCP85_2100.cov)
head(RCP85_2100.cov) #602, 668
nrow(RCP85_2100.cov)
nrow(pts.sp1_RCP85_2100) #2121
which(is.na(RCP85_2100.cov))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                               Generate pseudoabsence_RCP85_2100s                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Loading env variables
absence_RCP85_2100 <- randomPoints(RCP85_2100_covs, nrow(RCP85_2100.cov))

plot(absence_RCP85_2100)

nrow(absence_RCP85_2100)

abs.cov_RCP85_2100 <- raster::extract(RCP85_2100_covs, absence_RCP85_2100) #668

write.csv(absence_RCP85_2100,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/Absence_0.csv')

#--------------------------- Code the response
RCP85_2100.cov <- data.frame(RCP85_2100.cov)
RCP85_2100.cov$Sapidus_RCP85_2100 <- 1
abs.cov_RCP85_2100 <- data.frame(abs.cov_RCP85_2100)
abs.cov_RCP85_2100$Sapidus_RCP85_2100 <- 0

#--------------------------- And one to bind them
RCP85_2100_all.cov <- rbind(RCP85_2100.cov, abs.cov_RCP85_2100)
RCP85_2100_all.cov <- RCP85_2100_all.cov[complete.cases(RCP85_2100_all.cov), ]
nrow(RCP85_2100_all.cov)
head(RCP85_2100_all.cov)
which(is.na(RCP85_2100_all.cov))

RCP85_2100_xvars <- names(RCP85_2100_all.cov)[!(names(RCP85_2100_all.cov) == 'Sapidus_RCP85_2100')]
RCP85_2100_xvars

# Save files
write.csv(RCP85_2100_all.cov,'~/efc/R_Stat/2022_NEMA/SDM_Bayes/6_csv_files/RCP85_2100_all.cov.csv')

save(RCP85_2100_all.cov, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_all.cov.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_all.cov.rds') 

save(RCP85_2100_xvars, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_xvars.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_xvars.rds') 

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#                          Correlation and collinearity                        #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_all.cov.tiff",
     width=6, height=4, units="in", res=600)

corrplot(cor(RCP85_2100_all.cov[,1:3]),addCoef.col = 1,number.cex = 0.5,
         order = 'AOE',type = 'upper') 
dev.off()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# RCP85_2100                         BART MODELLING                               #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#--------------------------- Variable selection and returns the model
RCP85_2100.model <- bart.step(x.data = RCP85_2100_all.cov[, RCP85_2100_xvars],
                              y.data = RCP85_2100_all.cov[, 'Sapidus_RCP85_2100'],
                              full = TRUE,
                              quiet = TRUE)

save(RCP85_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_model.rds') 

#--------------------------- Retune: cross-validation
RET_RCP85_2100.model<-retune(RCP85_2100.model, reps = 10) 
summary(RET_RCP85_2100.model)

save(RET_RCP85_2100.model, 
     file='~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_retune.model.rds')
load('~/efc/R_Stat/2022_NEMA/SDM_Bayes/5_R_files/RCP85_2100_retune.model.rds') 

#--------------------------- Spatial prediction
RCP85_2100_Pred.map <- predict(
  object = RET_RCP85_2100.model,
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
points(pts.sp1_RCP85_2100[!is.na(raster::extract(RCP85_2100_Pred.map[[1]], 
                                                 SpatialPoints(pts.sp1_RCP85_2100[,1:2]))),],
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
varimp.diag(RCP85_2100_all.cov[,RCP85_2100_xvars], RCP85_2100_all.cov[,'Sapidus_RCP85_2100'], 
            ri.data = NULL, iter = 50, quiet = FALSE)
dev.off()

#~~~~~~~~~~ Varimp
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_varimpEU.tiff",
     width=6, height=4, units="in", res=600)
#varimp(RET_RCP85_2100.model, plots=TRUE)
varimpEU(RET_RCP85_2100.model, plots=TRUE)
dev.off()

#~~~~~~~~~~ Summary
tiff(file="~/efc/R_Stat/2022_NEMA/SDM_Bayes/3_Plots/RCP85_2100_Summary.tiff",
     width=8, height=6, units="in", res=600)
summary.bartEU(RET_RCP85_2100.model)
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

# RCP85_2100 Salinity
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

# RCP85_2100 Current
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

# Computation time
Time.masked <- difftime(Sys.time(),Start,units="min") # Time difference
Time.masked # 0.1073102 mins

#
#----------------------------------- / / --------------------------------------#
#

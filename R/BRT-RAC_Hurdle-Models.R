
## "Species distribution of anchovy, sardine & sardinella with delta log-normal BRT-RAC models" 
## Laura Julià Melis (lauraj@icm.csic.es)


#############
# 1. SET UP #
#############

# 1.1. Set working directory.
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# 1.2. Load necessary libraries.
library(dismo)
library(raster)
library(ggplot2)
library(gridExtra)
library(gbm)
library(fields)
library(gstat)
library(sf)
library(spdep)
library(maptools)
library(ggridges)
library(colorspace)




##########################################
# 2. READ SPECIES AND ENVIRONMENTAL DATA #
##########################################

# If, MEDITS Data:
# ----------------
load("../data/Anchovy_environment_MEDSEA.RData")
load("../data/Sardina_environment_MEDSEA.RData")
sardine_env <- sar_env
rm(list=setdiff(ls(), c("anc_env", "sardine_env", "depth", 
                        "npp_year", "sal_year","tem_year")))
load("../data/Sardinella_environment_MEDSEA.RData")
sardinella_env <- sar_env
rm(list=setdiff(ls(), c("anc_env", "sardine_env", "sardinella_env", "depth",
                        "npp_year", "sal_year", "tem_year")))


# 2.1. In the data frames with the sample data, we leave only the necessary columns
sample_anchovy <- anc_env[, c(2,10,11,5:9)]
sample_anchovy$Biomasa <- sample_anchovy$Biomasa/0.011 # Catchability coefficient q= 0.011

sample_sardine <- sardine_env[, c(2,10,11,5:9)]
sample_sardine$Biomasa <- sample_sardine$Biomasa/0.016 # Catchability coefficient q= 0.016

sample_sardinella <- sardinella_env[, c(2,10,11,5:9)]
sample_sardinella$Biomasa <- sample_sardinella$Biomasa/0.0135 # Catchability coefficient q=(0.011+0.016)/2
remove(list=c("anc_env", "sardine_env", "sardinella_env"))


# IF, VMS Data:
# -------------
load("../data/CORYs_environment_MEDSEA.RData")
df_env_ <- df_env_[complete.cases(df_env_[6:9]),]
rm(list=setdiff(ls(), c("df_env_", "depth","npp_year", "sal_year", "tem_year")))

# 2.1. We separate data frame by species and leave only necessary variables.
sample_anchovy <- df_env_[, c(2,10,11,3,6:9)]
sample_sardine <- df_env_[, c(2,10,11,4,6:9)]
sample_sardinella <- df_env_[, c(2,10,11,5,6:9)]
colnames(sample_anchovy)[2:4] <- c("lon", "lat", "Biomasa")
colnames(sample_sardine)[2:4] <- c("lon", "lat", "Biomasa")
colnames(sample_sardinella)[2:4] <- c("lon", "lat", "Biomasa")
rm("df_env_")


# FROM HERE, IT IS THE SAME FOR BOTH MEDITS AND VMS: 
# --------------------------------------------------
# 2.2. We remove points that are outside the GSA06.
gsa <- readShapePoly("../Shapefiles GSA-ZEPA/GSAs_simplified.shp",
                     proj4string=CRS("+init=epsg:4326"))
gsa06 <- subset(gsa, gsa@data$SECT_COD %in% c("GSA06"))

# ------------------------- anchovy ------------------------- #
sample_anchovy_gsa_sf <- as(st_as_sf(x = sample_anchovy, coords = c("lon", "lat"), 
                                     crs="+proj=longlat +datum=WGS84"),"Spatial")
sample_anchovy_gsa_sf <- sample_anchovy_gsa_sf[gsa06,]
sample_anchovy <- as.data.frame(sample_anchovy_gsa_sf)
sample_anchovy <- sample_anchovy[,c(1,7,8,2:6)]
colnames(sample_anchovy)[2:3] <- c("lon", "lat")

# ------------------------- sardine ------------------------- #
sample_sardine_gsa_sf <- as(st_as_sf(x = sample_sardine, coords = c("lon", "lat"), 
                                     crs="+proj=longlat +datum=WGS84"),"Spatial")
sample_sardine_gsa_sf <- sample_sardine_gsa_sf[gsa06,]
sample_sardine <- as.data.frame(sample_sardine_gsa_sf)
sample_sardine <- sample_sardine[,c(1,7,8,2:6)]
colnames(sample_sardine)[2:3] <- c("lon", "lat")

# ------------------------- sardinella ------------------------- #
sample_sardinella_gsa_sf <- as(st_as_sf(x = sample_sardinella,
                                        coords = c("lon", "lat"),
                                        crs="+proj=longlat +datum=WGS84"),
                               "Spatial")
sample_sardinella_gsa_sf <- sample_sardinella_gsa_sf[gsa06,]
sample_sardinella <- as.data.frame(sample_sardinella_gsa_sf)
sample_sardinella <- sample_sardinella[,c(1,7,8,2:6)]
colnames(sample_sardinella)[2:3] <- c("lon", "lat")

remove(list=c("gsa", "sample_anchovy_gsa_sf", "sample_sardine_gsa_sf",
              "sample_sardinella_gsa_sf"))


# 2.3. We create a presence/absence column (1/0).
sample_anchovy$presence <- ifelse(sample_anchovy$Biomasa == 0, 0, 1)
sample_sardine$presence <- ifelse(sample_sardine$Biomasa == 0, 0, 1)
sample_sardinella$presence <- ifelse(sample_sardinella$Biomasa == 0, 0, 1)


# 2.4. We create dataframes with only positive values to model the biomass.
sample_anchovy_bio <- subset(sample_anchovy, presence ==1, select= -presence)
sample_sardine_bio <- subset(sample_sardine, presence ==1, select= -presence)
sample_sardinella_bio <- subset(sample_sardinella, presence ==1, select= -presence)


# 2.5. We can save the data in a new workspace.
save.image(file = "data/data.RData")
load("data/data.RData") # run to load the workspace




###########################
# 3. EXPLORATORY ANALYSIS # 
###########################

# 3.1. Map of the study area.
## Load the SpatialPolygonsDataframe of Spain and France
spain <- raster::getData('GADM',country="ESP",level=0)
france <- getData('GADM',country="FRA",level=0)

## Loas ZEPA shapefile
zepa <-  readShapePoly("../Shapefiles GSA-ZEPA/ZEPA_Ebre_WGS84.shp",
                       proj4string=CRS("+init=epsg:4326"))

## Create plot
plot(spain, xlim=c(-1,4), ylim=c(37,42), axes=TRUE,col='lightgray')
plot(gsa06, add=TRUE,col=alpha('olivedrab3',0.4))
points(sample_anchovy$lon, sample_anchovy$lat, pch=19, cex=0.5, col ="red")
contour(depth, add=T, col=sequential_hcl(8))
legend(2.9, 37.8, legend=c("ZEPA ES0000512", "Profundidad"), 
       col=c(alpha('olivedrab3',0.4), sequential_hcl(8)), lty=c(NA,1),pch=c(16,NA), cex=0.8, bty="n")
text(x=2.3, y=42, labels= "Cabo de Creus", cex= 0.8, font=4)
text(x=-0.15, y=40.73, labels= "Delta del río Ebro", cex= 0.8, font=4)
text(x=-0.8, y=38.8, labels= "Cabo de la Nao", cex= 0.8, font=4)
text(x=-1.7, y=37.7, labels= "Cabo de Palos", cex= 0.8, font=4)
box()


# 3.2. Response variable distribution.
par(mfrow=c(3,3))
## Biomass density by species.
plot(density(sample_anchovy$Biomasa),
     main="Anchoa", xlab = "Biomasa (kg)", ylab= "Densidad")
plot(density(sample_sardine$Biomasa), 
     main="Sardina", xlab = "Biomasa (kg)", ylab= "Densidad")
plot(density(sample_sardinella$Biomasa), 
     main="Alacha", xlab = "Biomasa (kg)", ylab= "Densidad")

## Biomass density removing zeros.
plot(density(sample_anchovy_bio$Biomasa), 
     main="Anchoa", xlab = "Biomasa positiva (kg)", ylab= "Densidad")
plot(density(sample_sardine_bio$Biomasa), 
     main="Sardina", xlab = "Biomasa positiva  (kg)", ylab= "Densidad")
plot(density(sample_sardinella_bio$Biomasa),
     main="Alacha", xlab = "Biomasa positiva (kg)", ylab= "Densidad")

## Biomass density log-transformed.
plot(density(log(sample_anchovy_bio$Biomasa)), 
     main="Anchoa", xlab = "log(Biomasa)", ylab= "Densidad")
plot(density(log(sample_sardine_bio$Biomasa)),
     main="Sardina", xlab = "log(Biomasa)", ylab= "Densidad")
plot(density(log(sample_sardinella_bio$Biomasa)),
     main="Alacha", xlab = "log(Biomasa)", ylab= "Densidad")
par(mfrow=c(1,1))


# 3.3. Temporal trends of log-transformed biomass.

## First the average biomass per year is calculated
df1 <-aggregate(sample_anchovy$Biomasa, list(sample_anchovy$year), FUN=mean) 
df2 <-aggregate(sample_sardine$Biomasa, list(sample_sardine$year), FUN=mean) 
df3 <-aggregate(sample_sardinella$Biomasa, list(sample_sardinella$year), FUN=mean) 
colnames(df1) <- colnames(df2) <- colnames(df3) <- c("Year", "Biomass")
df1$Year <- df2$Year <- df3$Year <- 1998:2017

## Then we create a plot per species.
p1 <- ggplot(df1, aes(Year, Biomass)) + 
  geom_point() + geom_smooth() + theme_bw() +
  labs(subtitle= "Anchoa", y= "Biomasa (kg)", x = "Año", tag="A") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_x_continuous(breaks=c(1998,2008,2017))

p2 <- ggplot(df2, aes(Year, Biomass)) + 
  geom_point() + geom_smooth() + theme_bw()+
  labs(subtitle= "Sardina", y= "Biomasa (kg)", x = "Año", tag="B") + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+ 
  scale_x_continuous(breaks=c(1998,2008,2017))

p3 <- ggplot(df3, aes(Year, Biomass)) +
  geom_point() + geom_smooth() + theme_bw()+
  labs(subtitle= "Alacha", y= "Biomasa (kg)", x = "Año", tag="C")+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))+ 
  scale_x_continuous(breaks=c(1998,2008,2017))

grid.arrange(p1, p2, p3, ncol=3)
remove(list=c("p1", "p2", "p3", "df1", "df2", "df3"))

# 3.4. Plots of the environmental variables.

## First, cut rasters by GSA06.
for (i in 1:length(sal_year)){
  if(i == 1){
    sal_year.gsa <- mask(sal_year[[i]], gsa06)
    tem_year.gsa <- mask(tem_year[[i]], gsa06)
    npp_year.gsa <- mask(npp_year[[i]], gsa06)
    depth.gsa <- mask(depth, gsa06)
    
  } else{
    sal_year.gsa <- stack(sal_year.gsa, mask(sal_year[[i]], gsa06))
    tem_year.gsa <- stack(tem_year.gsa, mask(tem_year[[i]], gsa06))
    npp_year.gsa <- stack(npp_year.gsa, mask(npp_year[[i]], gsa06))
  }
}

## Then, plot.
# ------------------------ depth ------------------------ #
plot(depth.gsa, col=tim.colors(120)[1:120], axes=T, legend=T)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

# ---------------------- salinity ----------------------- #
zmin <- Inf
zmax <- 0
for(i in 1:nlayers(sal_year.gsa)){
  aux.min <- min(values(sal_year.gsa[[i]]), na.rm = T)
  aux.max <- max(values(sal_year.gsa[[i]]), na.rm = T)
  
  zmax <- ifelse(aux.max > zmax, aux.max, zmax)
  zmin <- ifelse(aux.min < zmin, aux.min, zmin)
}

par(mfrow=c(5,5), mar = c(1,1,2,2.2))
for(i in 1:nlayers(sal_year.gsa)){
  if(i != nlayers(sal_year.gsa)){
    plot(sal_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(25)[1:25], 
         main= paste0("Año ", 1997+i), axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(sal_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= paste0("Año ", 1997+i), axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# --------------------- temperature --------------------- #
zmin <- Inf
zmax <- 0
for(i in 1:nlayers(tem_year.gsa)){
  aux.min <- min(values(tem_year.gsa[[i]]), na.rm = T)
  aux.max <- max(values(tem_year.gsa[[i]]), na.rm = T)
  
  zmax <- ifelse(aux.max > zmax, aux.max, zmax)
  zmin <- ifelse(aux.min < zmin, aux.min, zmin)
}

par(mfrow=c(5,5), mar = c(1,1,2,2.2))
for(i in 1:nlayers(tem_year.gsa)){
  if(i != nlayers(tem_year.gsa)){
    plot(tem_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(25)[1:25], 
         main= paste0("Año ", 1997+i), axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(tem_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= paste0("Año ", 1997+i), axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# ------------------------- npp ------------------------- #
zmin <- Inf
zmax <- 0
for(i in 1:nlayers(npp_year.gsa)){
  aux.min <- min(values(npp_year.gsa[[i]]), na.rm = T)
  aux.max <- max(values(npp_year.gsa[[i]]), na.rm = T)
  
  zmax <- ifelse(aux.max > zmax, aux.max, zmax)
  zmin <- ifelse(aux.min < zmin, aux.min, zmin)
}

par(mfrow=c(5,5), mar = c(1,1,2,2.2))
for(i in 1:nlayers(npp_year.gsa)){
  if(i != nlayers(npp_year.gsa)){
    plot(npp_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(25)[1:25], 
         main= paste0("Año ", 1997+i), axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(npp_year.gsa[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= paste0("Año ", 1997+i), axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}




#################################
# 4. DELTA LOG-NORMAL BRT MODEL # 
#################################

# Delta log-normal boosted regression trees split zero-inflated, long-tailed distribution data
# (most trawls caught nothing, and very few caught many specimens) into zero/non-zero catches,
# converted to 0 or 1, and log-normalized non-zero catches (i.e. biomass).

# 4.1. Models for probability of ocurrence (Binomial BRT models). 
# ===============================================================

## 4.1.1. Models with ONLY environmental variables
brt.ber.anx <- gbm.step(data= sample_anchovy, gbm.x = c(1,5:8), gbm.y = 9, 
                        family = "bernoulli", silent=T)
brt.ber.sar <- gbm.step(data= sample_sardine, gbm.x = c(1,5:8), gbm.y = 9,
                        family = "bernoulli", silent=T)
brt.ber.sardinella <- gbm.step(data= sample_sardinella, gbm.x = c(1,5:8), 
                               gbm.y = 9, family = "bernoulli", silent=T)


## 4.1.2. Autocovariate calculation from the residuals of the previous models.
coords <- data.frame("lon" = sample_anchovy$lon, "lat" = sample_anchovy$lat)

# ------------------------- anchovy ------------------------- #
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.ber.anx))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.ber.anx <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T) # focal calculation neighbors of 1st order
sample_anchovy$autocovariate <- extract(spatial.ber.anx, coords) # We add the autocovariate to the dataframe

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_anchovy, coords = c("lon", "lat"),
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.ber.anx, gs)
spatial.ber.anx <- mask(idw, depth)


# ------------------------- sardine ------------------------- #
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.ber.sar))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.ber.sar <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T)
sample_sardine$autocovariate <- extract(spatial.ber.sar, coords)

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_sardine, coords = c("lon", "lat"), 
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.ber.sar, gs)
spatial.ber.sar <- mask(idw, depth)


# ------------------------- sardinella ------------------------- #
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.ber.sardinella))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.ber.sardinella <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T)
sample_sardinella$autocovariate <- extract(spatial.ber.sardinella, coords)

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_sardinella, coords = c("lon", "lat"), 
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.ber.sardinella, gs)
spatial.ber.sardinella <- mask(idw, depth)

rm(list=c("coords", "rast", "xy_residuals", "sf.data", "gs", "idw", 
          "brt.ber.anx", "brt.ber.sar", "brt.ber.sardinella"))


## 4.1.3. Models with the environmental data + the autocovariate.
brt.ber.anx <- gbm.step(data= sample_anchovy, gbm.x = c(1,5:8,10), gbm.y = 9,
                        family = "bernoulli", silent=T)
brt.ber.sar <- gbm.step(data= sample_sardine, gbm.x = c(1,5:8,10), gbm.y = 9,
                        family = "bernoulli", silent=T)
brt.ber.sardinella <- gbm.step(data= sample_sardinella, gbm.x = c(1,5:8,10), 
                               gbm.y = 9, family = "bernoulli", silent=T)


## 4.1.4. Deviance explained by the model (performance evaluation of the models)
devexpl.ber.anx <- ((brt.ber.anx$self.statistics$null-brt.ber.anx$self.statistics$resid)/brt.ber.anx$self.statistics$null)*100
devexpl.ber.sar <- ((brt.ber.sar$self.statistics$null-brt.ber.sar$self.statistics$resid)/brt.ber.sar$self.statistics$null)*100
devexpl.ber.sardinella <- ((brt.ber.sardinella$self.statistics$null-brt.ber.sardinella$self.statistics$resid)/brt.ber.sardinella$self.statistics$null)*100
devexpl.ber.anx; devexpl.ber.sar; devexpl.ber.sardinella


## 4.1.5. Most influential variables in each model.
brt.ber.anx$contributions
brt.ber.sar$contributions
brt.ber.sardinella$contributions


## 4.1.6. We save results in a new workspace.
dir.create("./outputs") # it creates a folder for output files in our current directory
save.image(file = "results.RData")
load("outputs/results.RData")# run to load the workspace



# 4.2. Models for the non-zero log(biomass) observations (Gaussian BRT models).
# =============================================================================

## 4.2.1. We log-transform the biomasses.
sample_anchovy_bio$logbio <- log(sample_anchovy_bio$Biomasa)
sample_sardine_bio$logbio <- log(sample_sardine_bio$Biomasa)
sample_sardinella_bio$logbio <- log(sample_sardinella_bio$Biomasa)


## 4.2.2.  Models with ONLY environmental variables.
brt.gau.anx <- gbm.step(data= sample_anchovy_bio, gbm.x = c(1,5:8), gbm.y = 9,
                        family = "gaussian", silent=T)
brt.gau.sar <- gbm.step(data= sample_sardine_bio, gbm.x = c(1,5:8), gbm.y = 9,
                        family = "gaussian", silent=T)
brt.gau.sardinella <- gbm.step(data= sample_sardinella_bio, gbm.x = c(1,5:8), 
                               gbm.y = 9, family = "gaussian")



## 4.2.3. Autocovariate calculation from the residuals of the previous models.

# ------------------------- anchovy ------------------------- #
coords <- data.frame("lon" = sample_anchovy_bio$lon, "lat" = sample_anchovy_bio$lat)
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.gau.anx))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.gau.anx <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T) # focal calculation neighbors of 1st order
sample_anchovy_bio$autocovariate <- extract(spatial.gau.anx, coords) # We add the autocovariate to the dataframe

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_anchovy_bio, coords = c("lon", "lat"), 
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.gau.anx, gs)
spatial.gau.anx <- mask(idw, depth)


# ------------------------- sardine ------------------------- #
coords <- data.frame("lon" = sample_sardine_bio$lon, "lat" = sample_sardine_bio$lat)
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.gau.sar))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.gau.sar <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T)
sample_sardine_bio$autocovariate <- extract(spatial.gau.sar, coords)

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_sardine_bio, coords = c("lon", "lat"), 
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.gau.sar, gs)
spatial.gau.sar <- mask(idw, depth)


# ------------------------- sardinella ------------------------- #
coords <- data.frame("lon" = sample_sardinella_bio$lon, 
                     "lat" = sample_sardinella_bio$lat)
rast <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
               ymn = extent(depth)[3], ymx = extent(depth)[4], 
               nrows = dim(depth)[1], ncols = dim(depth)[2]) 
xy_residuals <- cbind(coords, resid(brt.gau.sardinella))
rast[cellFromXY(rast, xy_residuals)] <- xy_residuals[,3]        
spatial.gau.sardinella <- raster::focal(rast, w = matrix(1, nr= 9, nc=9), na.rm =T)
sample_sardinella_bio$autocovariate <- extract(spatial.gau.sardinella, coords)

## IDW interpolation of the autocovariate:
sf.data <- as(st_as_sf(sample_sardinella_bio, coords = c("lon", "lat"),
                       crs= "+proj=longlat +datum=WGS84"), "Spatial")
gs <- gstat::gstat(formula = autocovariate ~ 1, data = sf.data)
idw <- interpolate(spatial.gau.sardinella, gs)
spatial.gau.sardinella <- mask(idw, depth)

rm(list=c("coords", "rast", "xy_residuals", "sf.data", "gs", "idw", 
          "brt.gau.anx", "brt.gau.sar", "brt.gau.sardinella"))


## 4.2.4. Models with the environmental data + the autocovariate.
brt.gau.anx <- gbm.step(data= sample_anchovy_bio, gbm.x = c(1,5:8,10), gbm.y = 9,
                        family = "gaussian", silent=T)
brt.gau.sar <- gbm.step(data= sample_sardine_bio, gbm.x = c(1,5:8,10), gbm.y = 9,
                        family = "gaussian", silent=T)
brt.gau.sardinella <- gbm.step(data= sample_sardinella_bio, gbm.x = c(1,5:8,10),
                               gbm.y = 9, family = "gaussian", silent=T)

## 4.2.5. Deviance explained by the model (performance evaluation of the models)
devexpl.gau.anx <- ((brt.gau.anx$self.statistics$null-brt.gau.anx$self.statistics$resid)/brt.gau.anx$self.statistics$null)*100
devexpl.gau.sar <- ((brt.gau.sar$self.statistics$null-brt.gau.sar$self.statistics$resid)/brt.gau.sar$self.statistics$null)*100
devexpl.gau.sardinella <- ((brt.gau.sardinella$self.statistics$null-brt.gau.sardinella$self.statistics$resid)/brt.gau.sardinella$self.statistics$null)*100
devexpl.gau.anx; devexpl.gau.sar; devexpl.gau.sardinella 


## 4.2.6. Most influential variables in each model.
brt.gau.anx$contributions
brt.gau.sar$contributions
brt.gau.sardinella$contributions

## 4.2.7. We save results in a the workspace.
save.image(file = "results.RData")



#####################################
# 5. PROVES AUTOCORRELACIÓ ESPACIAL #
#####################################

# 5.1. Occurrence models.
# =======================
coords <- as.matrix(cbind(sample_anchovy$lon, sample_anchovy$lat))
nb <- dnearneigh(as.matrix(coords), 0, 20)
listw <- nb2listw(nb)

nsim <- 999
set.seed(1234)

MoranI.anx.ber <- moran.mc(residuals(brt.ber.anx), listw=listw, nsim=nsim) 
MoranI.sar.ber <- moran.mc(residuals(brt.ber.sar), listw=listw, nsim=nsim) 
MoranI.sardinella.ber <- moran.mc(residuals(brt.ber.sardinella), 
                                  listw=listw, nsim=nsim) 
MoranI.anx.ber; MoranI.sar.ber; MoranI.sardinella.ber




# 5.2. Biomass models.
# ====================
nsim <- 999
set.seed(1234)

coords <- as.matrix(cbind(sample_anchovy_bio$lon, sample_anchovy_bio$lat))
nb <- dnearneigh(as.matrix(coords), 0, 20)
listw <- nb2listw(nb)
MoranI.anx.gau <- moran.mc(residuals(brt.gau.anx), listw=listw, nsim=nsim)

coords <- as.matrix(cbind(sample_sardine_bio$lon, sample_sardine_bio$lat))
nb <- dnearneigh(as.matrix(coords), 0, 20)
listw <- nb2listw(nb)
MoranI.sar.gau  <- moran.mc(residuals(brt.gau.sar), listw=listw, nsim=nsim)

coords <- as.matrix(cbind(sample_sardinella_bio$lon, sample_sardinella_bio$lat))
nb <- dnearneigh(as.matrix(coords), 0, 20)
listw <- nb2listw(nb)
MoranI.sardinella.gau  <- moran.mc(residuals(brt.gau.sardinella), listw=listw, nsim=nsim)
MoranI.anx.gau ; MoranI.sar.gau ; MoranI.sardinella.gau 

remove(list=c("coords", "nb", "listw", "nsim"))

## Save workspace.
save.image(file = "results.RData")





#########################
# 6. MAKING PREDICTIONS # 
#########################

# 6.1. Function to combine occurrence and biomass predictions.
# ============================================================

pred <- function(model1, model2, spatial1, spatial2, depth){
  
  year <- raster(xmn = extent(depth)[1], xmx = extent(depth)[2], 
                 ymn = extent(depth)[3], ymx = extent(depth)[4], 
                 nrows = dim(depth)[1], ncols = dim(depth)[2]) 
  
  for(i in 1:length(sal_year)){
    # 0. Create raster with the year
    values(year) <- 1997+i
    
    # 1. Prediction of ocurrence.
    predictors1 <- stack(year, sal_year[[i]], tem_year[[i]],
                         npp_year[[i]], depth, spatial1)
    names(predictors1) <- c("year", "SSS", "SST", "NPP", "Depth", "autocovariate")
    pred.ber <- predict(predictors1, model1, type="response", 
                        n.trees=model1$n.trees,
                        shrinkage= 0.01, distribution="bernoulli")
    if(i==1){
      pred.ocurrence <- pred.ber
    } else{
      pred.ocurrence <- stack(pred.ocurrence, pred.ber)
    }
    # 2. Prediction of log-biomass
    predictors2 <- stack(year, sal_year[[i]], tem_year[[i]], 
                         npp_year[[i]], depth, spatial2)
    names(predictors2) <- c("year", "SSS", "SST", "NPP", "Depth", "autocovariate")
    pred.gau <- predict(predictors2, model2, type="response", 
                        n.trees=model2$n.trees,
                        shrinkage= 0.01, distribution="gaussian")
    
    # 3. Reverse log-transformation
    values(pred.gau) <- exp(values(pred.gau))
    
    # 4. Combine results.
    if(i==1){
      pred <- pred.ber*pred.gau
    } else{
      pred <- stack(pred, pred.ber*pred.gau)
    }
  }
  return(list(pred.ocurrence, pred))
}


# 6.2. We calculate predictions for each species.
# ===============================================

##  6.2.1. First we trim the bathymetry to <= -850.
matrix <- cbind(coordinates(depth), depth=raster::getValues(depth))
matrix <- as.data.frame(matrix[!is.na(matrix[,3]),])
matrix <- subset(matrix, matrix$depth  >= -850)
depth2 <- rasterFromXYZ(matrix)
rm(matrix)

##  6.2.2. Now we predict:
pred.anx <- pred(model1 = brt.ber.anx,
                 model2 = brt.gau.anx,
                 spatial1 = spatial.ber.anx, 
                 spatial2 = spatial.gau.anx, depth2)
pred.ber.anx <- pred.anx[[1]]
pred.anx <- pred.anx[[2]]

pred.sar <- pred(model1 = brt.ber.sar,
                 model2 = brt.gau.sar,
                 spatial1 = spatial.ber.sar, 
                 spatial2 = spatial.gau.sar, depth2)
pred.ber.sar <- pred.sar[[1]]
pred.sar <- pred.sar[[2]]

pred.sardinella <- pred(model1 = brt.ber.sardinella, 
                        model2 = brt.gau.sardinella,
                        spatial1 = spatial.ber.sardinella, 
                        spatial2 = spatial.gau.sardinella, depth2)
pred.ber.sardinella <- pred.sardinella[[1]]
pred.sardinella <- pred.sardinella[[2]]

names(pred.anx) <- names(pred.sar) <- names(pred.sardinella) <- paste0("Year:", 1998:2019)
names(pred.ber.anx) <- names(pred.ber.sar) <- names(pred.ber.sardinella) <- paste0("Year:", 1998:2019)




############################################
# 7. CUT MAPS FOR GSA06 and THE DELTA ZEPA # 
############################################

# ------------------------- anchovy ------------------------- #
for (i in 1:nlayers(pred.anx)){
  if(i == 1){
    pred.anx.gsa <- mask(pred.anx[[i]], gsa06)
    pred.anx.zepa <- mask(pred.anx[[i]], zepa)
    pred.anx.ber.gsa <- mask(pred.ber.anx[[i]], gsa06)
    pred.anx.ber.zepa <- mask(pred.ber.anx[[i]], zepa)
  } else{
    pred.anx.gsa <- stack(pred.anx.gsa, mask(pred.anx[[i]], gsa06))
    pred.anx.zepa <- stack(pred.anx.zepa, mask(pred.anx[[i]], zepa))
    pred.anx.ber.gsa <- stack(pred.anx.ber.gsa, 
                              mask(pred.ber.anx[[i]], gsa06))
    pred.anx.ber.zepa <- stack(pred.anx.ber.zepa, 
                               mask(pred.ber.anx[[i]], zepa))
  }
}


# ------------------------- sardine ------------------------- #
for (i in 1:nlayers(pred.sar)){
  if(i == 1){
    pred.sar.gsa <- mask(pred.sar[[i]], gsa06)
    pred.sar.zepa <- mask(pred.sar[[i]], zepa)
    pred.sar.ber.gsa <- mask(pred.ber.sar[[i]], gsa06)
    pred.sar.ber.zepa <- mask(pred.ber.sar[[i]], zepa)
  } else{
    pred.sar.gsa <- stack(pred.sar.gsa, mask(pred.sar[[i]], gsa06))
    pred.sar.zepa <- stack(pred.sar.zepa, mask(pred.sar[[i]], zepa))
    pred.sar.ber.gsa <- stack(pred.sar.ber.gsa,
                              mask(pred.ber.sar[[i]], gsa06))
    pred.sar.ber.zepa <- stack(pred.sar.ber.gsa,
                               mask(pred.ber.sar[[i]], zepa))
  }
}


# ------------------------- sardinella ------------------------- #
for (i in 1:nlayers(pred.sardinella)){
  if(i == 1){
    pred.sardinella.gsa <- mask(pred.sardinella[[i]], gsa06)
    pred.sardinella.zepa <- mask(pred.sardinella[[i]], zepa)
    pred.sardinella.ber.gsa <- mask(pred.ber.sardinella[[i]], gsa06)
    pred.sardinella.ber.zepa <- mask(pred.ber.sardinella[[i]], zepa)
  } else{
    pred.sardinella.gsa <- stack(pred.sardinella.gsa, 
                                 mask(pred.sardinella[[i]], gsa06))
    pred.sardinella.zepa <- stack(pred.sardinella.zepa, 
                                  mask(pred.sardinella[[i]], zepa))
    pred.sardinella.ber.gsa <- stack(pred.sardinella.ber.gsa, 
                                     mask(pred.ber.sardinella[[i]], gsa06))
    pred.sardinella.ber.zepa <- stack(pred.sardinella.ber.zepa, 
                                      mask(pred.ber.sardinella[[i]], zepa))
  }
}

pred.anx <- pred.anx.gsa
pred.sar <- pred.sar.gsa
pred.sardinella <- pred.sardinella.gsa
pred.ber.anx <- pred.anx.ber.gsa
pred.ber.sar <- pred.sar.ber.gsa
pred.ber.sardinella <- pred.sardinella.ber.gsa
rm(list=c("i", "pred.anx.gsa", "pred.sar.gsa", "pred.sardinella.gsa",
          "pred.anx.ber.gsa", "pred.sar.ber.gsa", "pred.sardinella.ber.gsa"))

## Save workspace.
save.image(file = "results.RData")




##################################################################
# 8. CONSTRUCTION OF BIOMASS PREDICTION MAPS BY SPECIES AND YEAR # 
##################################################################

# 8.1. Evolution over the years.
# =============================

# ------------------------- anchovy ------------------------- #
zmin <- min(values(pred.anx), na.rm = T)
zmax <- max(values(pred.anx), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.anx)){
  if(i != nlayers(pred.anx)){
    plot(pred.anx[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.anx)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.anx[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.anx)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# ------------------------- sardine ------------------------- #
zmin <- min(values(pred.sar), na.rm = T)
zmax <- max(values(pred.sar), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.sar)){
  if(i != nlayers(pred.sar)){
    plot(pred.sar[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.sar)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.sar[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.sar)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# ------------------------- sardinella ------------------------- #
zmin <- min(values(pred.sardinella), na.rm = T)
zmax <- max(values(pred.sardinella), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.sardinella)){
  if(i != nlayers(pred.sardinella)){
    plot(pred.sardinella[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.sardinella)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.sardinella[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.sardinella)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}


# 8.2. Mean prediction (and SD).
# ==============================

# ------------------------- anchovy ------------------------- #
mean.pred.anx <- mean(pred.anx)
mean.zmax <- max(values(mean.pred.anx), na.rm = T)
sd.pred.anx <-calc(pred.anx, sd)
sd.zmax <- max(values(sd.pred.anx), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.anx,zlim=c(0, mean.zmax), col=tim.colors(120),
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.anx, zlim=c(0, sd.zmax), col=tim.colors(120), 
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)


# ------------------------- sardine ------------------------- #
mean.pred.sar <- mean(pred.sar)
mean.zmax <- max(values(mean.pred.sar), na.rm = T)
sd.pred.sar <-calc(pred.sar, sd)
sd.zmax <- max(values(sd.pred.sar), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.sar,zlim=c(0, mean.zmax), col=tim.colors(120)[1:120], 
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.sar, zlim=c(0, sd.zmax), col=tim.colors(120)[1:120], 
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)


# ------------------------- sardinella ------------------------- #
mean.pred.sardinella <- mean(pred.sardinella)
mean.zmax <- max(values(mean.pred.sardinella), na.rm = T)
sd.pred.sardinella <-calc(pred.sardinella, sd)
sd.zmax <- max(values(sd.pred.sardinella), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.sardinella,zlim=c(0, mean.zmax), col=tim.colors(120)[1:120], 
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.sardinella, zlim=c(0, sd.zmax), col=tim.colors(120)[1:120], 
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)


# 8.3. Maps of the residual autocovariate (RAC). 
# =============================================
par(mfrow=c(1,3))
plot(spatial.ber.anx, main = "Anchovy")
plot(spatial.ber.sar, main = "Sardine")
plot(spatial.ber.sardinella, main = "Alacha")

par(mfrow=c(1,3))
plot(spatial.gau.anx, main = "Anchovy")
plot(spatial.gau.sar, main = "Sardine")
plot(spatial.gau.sardinella, main = "Alacha")




#####################################################################
# 9. CONSTRUCTION OF OCCURRENCE PREDICTION MAPS BY SPECIES AND YEAR # 
#####################################################################

# 9.1. Evolution over the years.. 
# ===============================

# ------------------------- anchovy ------------------------- #
zmin <- min(values(pred.ber.anx), na.rm = T)
zmax <- max(values(pred.ber.anx), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.ber.anx)){
  if(i != nlayers(pred.ber.anx)){
    plot(pred.ber.anx[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.anx)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.ber.anx[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.anx)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# ------------------------- sardine ------------------------- #
zmin <- min(values(pred.ber.sar), na.rm = T)
zmax <- max(values(pred.ber.sar), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.ber.sar)){
  if(i != nlayers(pred.ber.sar)){
    plot(pred.ber.sar[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.sar)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.ber.sar[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.sar)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}

# ------------------------- sardinella ------------------------- #
zmin <- min(values(pred.ber.sardinella), na.rm = T)
zmax <- max(values(pred.ber.sardinella), na.rm = T)
par(mfrow=c(5,5), mar = c(1,1,1,1))
for(i in 1:nlayers(pred.ber.sardinella)){
  if(i != nlayers(pred.ber.sardinella)){
    plot(pred.ber.sardinella[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.sardinella)[[i]], axes=F, legend=F)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }else{
    plot(pred.ber.sardinella[[i]],zlim=c(zmin, zmax), col=tim.colors(120)[1:120], 
         main= names(pred.ber.sardinella)[[i]], axes=F, legend=TRUE)
    plot(spain, col='dark grey',add=T)
    plot(france, col='dark grey',add=T)
  }
}


# 9.2. Mean prediction (and SD). 
# ==============================

# ------------------------- anchovy ------------------------- #
mean.pred.ber.anx <- mean(pred.ber.anx)
mean.zmax <- max(values(mean.pred.ber.anx), na.rm = T)
sd.pred.ber.anx <-calc(pred.ber.anx, sd)
sd.zmax <- max(values(sd.pred.ber.anx), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.ber.anx,zlim=c(0, mean.zmax), col=tim.colors(120),
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.ber.anx, zlim=c(0, sd.zmax), col=tim.colors(120),
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)


# ------------------------- sardine ------------------------- #
mean.pred.ber.sar <- mean(pred.ber.sar)
mean.zmax <- max(values(mean.pred.ber.sar), na.rm = T)
sd.pred.ber.sar <-calc(pred.ber.sar, sd)
sd.zmax <- max(values(sd.pred.ber.sar), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.ber.sar,zlim=c(0, mean.zmax), col=tim.colors(120)[1:120], 
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.ber.sar, zlim=c(0, sd.zmax), col=tim.colors(120)[1:120],
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)


# ------------------------- sardinella ------------------------- #
mean.pred.ber.sardinella <- mean(pred.ber.sardinella)
mean.zmax <- max(values(mean.pred.ber.sardinella), na.rm = T)
sd.pred.ber.sardinella <-calc(pred.ber.sardinella, sd)
sd.zmax <- max(values(sd.pred.ber.sardinella), na.rm = T)

par(mfrow=c(1,2), mar= c(2,1.5,3,6))
plot(mean.pred.ber.sardinella,zlim=c(0, mean.zmax), col=tim.colors(120)[1:120],
     main= "Media", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)

plot(sd.pred.ber.sardinella, zlim=c(0, sd.zmax), col=tim.colors(120)[1:120],
     main= "Desviación típica", axes=F)
plot(spain, col='dark grey',add=T)
plot(france, col='dark grey',add=T)



#############################
# 10. SAVING PREDICION MAPS # 
#############################

# ------------------------- anchovy ------------------------- #

## OCURRENCE
writeRaster(mean.pred.ber.anx, "prediccions/anchovy/mean_ocurrence.tif", overwrite=T)
writeRaster(sd.pred.ber.anx, "prediccions/anchovy/sd_ocurrence.tif", overwrite=T)

for(i in 1:nlayers(pred.ber.anx)){
  path <- paste0("prediccions/anchovy/", 
                 substr(names(pred.ber.anx)[[i]], 6, 9), "_ocurrence.tif")
  writeRaster(pred.ber.anx[[i]], path, overwrite=T)
}

## BIOMASS
writeRaster(mean.pred.anx, "prediccions/anchovy/mean_biomass.tif", overwrite=T)
writeRaster(sd.pred.anx, "prediccions/anchovy/sd_biomass.tif", overwrite=T)

for(i in 1:nlayers(pred.anx)){
  path <- paste0("prediccions/anchovy/", 
                 substr(names(pred.anx)[[i]], 6, 9), "_biomass.tif")
  writeRaster(pred.anx[[i]], path, overwrite=T)
}

# ------------------------- sardine ------------------------- #

## OCURRENCE
writeRaster(mean.pred.ber.sar, "prediccions/sardine/mean_ocurrence.tif", overwrite=T)
writeRaster(sd.pred.ber.sar, "prediccions/sardine/sd_ocurrence.tif", overwrite=T)

for(i in 1:nlayers(pred.ber.sar)){
  path <- paste0("prediccions/sardine/", 
                 substr(names(pred.ber.sar)[[i]], 6, 9), "_ocurrence.tif")
  writeRaster(pred.ber.sar[[i]], path, overwrite=T)
}

## BIOMASS
writeRaster(mean.pred.sar, "prediccions/sardine/mean_biomass.tif", overwrite=T)
writeRaster(sd.pred.sar, "prediccions/sardine/sd_biomass.tif", overwrite=T)

for(i in 1:nlayers(pred.sar)){
  path <- paste0("prediccions/sardine/", 
                 substr(names(pred.sar)[[i]], 6, 9), "_biomass.tif")
  writeRaster(pred.sar[[i]], path, overwrite=T)
}



# ------------------------- sardinella ------------------------- #

## OCURRENCE
writeRaster(mean.pred.ber.sardinella, "prediccions/sardinella/mean_ocurrence.tif", overwrite=T)
writeRaster(sd.pred.ber.sardinella, "prediccions/sardinella/sd_ocurrence.tif", overwrite=T)

for(i in 1:nlayers(pred.ber.sardinella)){
  path <- paste0("prediccions/sardinella/", 
                 substr(names(pred.ber.sardinella)[[i]], 6, 9), "_ocurrence.tif")
  writeRaster(pred.ber.sardinella[[i]], path, overwrite=T)
}

## BIOMASS
writeRaster(mean.pred.sardinella, "prediccions/sardinella/mean_biomass.tif", overwrite=T)
writeRaster(sd.pred.sardinella, "prediccions/sardinella/sd_biomass.tif", overwrite=T)

for(i in 1:nlayers(pred.sardinella)){
  path <- paste0("prediccions/sardinella/", 
                 substr(names(pred.sardinella)[[i]], 6, 9), "_biomass.tif")
  writeRaster(pred.sardinella[[i]], path, overwrite=T)
}




############################
# 11. SENSITIVITY ANALYSIS # 
############################

# Here we create density plots of the log(biomass) for each year.

# ------------------------- anchovy ------------------------- #
df <- as.data.frame(pred.anx)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_anx_gsa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) + 
  labs(title = "Anchoa", tag="A") +
  theme_minimal() + labs(y="Densidad", x="log-biomasa") + theme(legend.position = "none")

df <- as.data.frame(pred.anx.zepa)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_anx_zepa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) + 
  labs(title = "Anchoa", tag="A") +
  theme_minimal() + labs(y="Densidad", x="log-biomasa") + theme(legend.position = "none")
rm(df)

# ------------------------- sardine ------------------------- #

df <- as.data.frame(pred.sar)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_sar_gsa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) + 
  labs(title = "Sardina", tag="B") +
  theme_minimal() + labs(y="", x="log-biomasa") + theme(legend.position = "none")

df <- as.data.frame(pred.sar.zepa)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_sar_zepa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) + 
  labs(title = "Sardina", tag="B") +
  theme_minimal() + labs(y="", x="log-biomasa") + theme(legend.position = "none")

rm(df)


# ------------------------- sardinella ------------------------- #

df <- as.data.frame(pred.sardinella)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_sardinella_gsa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) +
  labs(title = "Alacha", tag="C") +
  theme_minimal() + labs(y="", x="log-biomasa") + theme(legend.position = "none")


df <- as.data.frame(pred.sardinella.zepa)
df <- df[complete.cases(df), ]
colnames(df) <- 1998:2019
df <- data.frame(stack(df))
colnames(df) <- c("Biomass", "Year")


plot_sardinella_zepa <- ggplot(df, aes(x = log(Biomass), y = Year)) +
  geom_density_ridges(rel_min_height = 0.005, scale=3) + 
  labs(title = "Alacha", tag="C") +
  theme_minimal() + labs(y="", x="log-biomasa") + theme(legend.position = "none")
rm(df)


grid.arrange(plot_anx_gsa, plot_sar_gsa, plot_sardinella_gsa, ncol=3)
grid.arrange(plot_anx_zepa, plot_sar_zepa, plot_sardinella_zepa, ncol=3)



## Save workspace.
save.image(file = "results.RData")




#################################
# 12. FUNCTIONAL RESPONSE PLOTS # 
#################################

## OCURRENCE MODELS
par(mar=c(4,4,1,2))
gbm.plot(brt.ber.anx,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))
gbm.plot(brt.ber.sar,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))
gbm.plot(brt.ber.sardinella,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))

## BIOMASS MODELS
gbm.plot(brt.gau.anx,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))
gbm.plot(brt.gau.sar,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))
gbm.plot(brt.gau.sardinella,write.title = F, show.contrib = T, smooth=T,
         y.label = "fitted function",plot.layout=c(2,3))

############### Blog: Change analysis  ########## 
## 
## DATE CREATED: 06/19/2017
## DATE MODIFIED: 06/21/2018
## AUTHORS: Benoit Parmentier 
## Version: 1
## PROJECT: General purpose
## ISSUE: 
## TO DO: Make this a function later
##
## COMMIT: initial commit
##
## Links to investigate:
#
###################################################
#

###Loading R library and packages                                                      

library(sp) # spatial/geographfic objects and functions
library(rgdal) #GDAL/OGR binding for R with functionalities
library(spdep) #spatial analyses operations, functions etc.
library(gtools) # contains mixsort and other useful functions
library(maptools) # tools to manipulate spatial data
library(parallel) # parallel computation, part of base package no
library(rasterVis) # raster visualization operations
library(raster) # raster functionalities
library(forecast) #ARIMA forecasting
library(xts) #extension for time series object and analyses
library(zoo) # time series object and analysis
library(lubridate) # dates functionality
library(colorRamps) #contains matlab.like color palette
library(rgeos) #contains topological operations
library(sphet) #contains spreg, spatial regression modeling
library(BMS) #contains hex2bin and bin2hex, Bayesian methods
library(bitops) # function for bitwise operations
library(foreign) # import datasets from SAS, spss, stata and other sources
library(gdata) #read xls, dbf etc., not recently updated but useful
library(classInt) #methods to generate class limits
library(plyr) #data wrangling: various operations for splitting, combining data
#library(gstat) #spatial interpolation and kriging methods
library(readxl) #functionalities to read in excel type data
library(psych) #pca/eigenvector decomposition functionalities
library(sf) # spatial objects classes
library(plotrix) #various graphic functions e.g. draw.circle

###### Functions used in this script


############################################################################
#####  Parameters and argument set up ###########

out_suffix <- "change_" #output suffix for the files and ouptut folder #param 12

in_dir <- "/nfs/bparmentier-data/Data/blogs/blog1_basic_change_analysis/data/"
out_dir <- "/nfs/bparmentier-data/Data/blogs/blog1_basic_change_analysis/outputs"

file_format <- ".tif" #PARAM5
NA_flag_val <- -9999 #PARAM7
create_out_dir_param=TRUE #PARAM9

#new_strata_rita_10282017.shp
#nlcd_2006_RITA.tif
#nlcd_legend.txt
#df_modis_band_info.txt

infile_reflectance_date1 <- "mosaiced_MOD09A1_A2005265__006_reflectance_masked_RITA_reg_1km.tif"
infile_reflectance_date2 <- "mosaiced_MOD09A1_A2005273__006_reflectance_masked_RITA_reg_1km.tif"

###############################################
##### PART III: Band combination: Indices and thresholding for flood mapping ##############

  
###### Read in MOD09 reflectance images before and after Hurrican Rita.
r_before <- brick(file.path(in_dir,infile_reflectance_date1)) # Before RITA, Sept. 22, 2005.
r_after <- brick(file.path(in_dir,infile_reflectance_date2)) # After RITA, Sept 30, 2005.

names(r_before) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")
names(r_after) <- c("Red","NIR","Blue","Green","SWIR1","SWIR2","SWIR3")

plot(r_before)
plot(r_after)

#### Generate Color composites
### True color composite

# 4 figures arranged in 2 rows and 2 columns
#attach(mtcars)
par(mfrow=c(2,1))

plotRGB(r_before,
        r=1,
        g=4,
        b=3,
        scale=0.6,
        stretch="hist",
        main="before")

plotRGB(r_after,
        r=1,
        g=4,
        b=3,
        scale=0.6,
        stretch="hist",
        main="after")

### False color composite:

plotRGB(r_before,
        r=2,
        g=1,
        b=4,
        scale=0.6,
        stretch="hist")

plotRGB(r_after,
        r=2,
        g=1,
        b=4,
        scale=0.6,
        stretch="hist")

### NIR experiment with threshold to  map water/flooding

############## Generating indices based on raster algebra of original bands
## Let's generate a series of indices, we list a few possibility from the literature.
#1) NDVI = (NIR - Red)/(NIR+Red)
#2) NDWI = (Green - NIR)/(Green + NIR)
#3) MNDWI = Green - SWIR2 / Green + SWIR2
#4) NDWI2 (LSWIB5) =  (NIR - SWIR1)/(NIR + SWIR1)
#5) LSWI (LSWIB5) =  (NIR - SWIR2)/(NIR + SWIR2)

names(r_before)
r_before_NDVI <- (r_before$NIR - r_before$Red) / (r_before$NIR + r_before$Red)
r_after_NDVI <- subset(r_after,"NIR") - subset(r_after,"Red")/(subset(r_after,"NIR") + subset(r_after,"Red"))

r_before_NDWI <- (r_before$Red - r_before$SWIR2) / (r_before$Red + r_before$SWIR2)
r_after_NDWI <- (r_after$Red - r_after$SWIR2) / (r_after$Red + r_after$SWIR2)

plot(r_before_NDVI,zlim=c(-1,1),col=matlab.like(255))
plot(r_after_NDVI,zlim=c(-1,1),col=matlab.like2(255))

plot(r_before_NDWI,zlim=c(-1,1),col=matlab.like(255))
plot(r_after_NDWI,zlim=c(-1,1),col=matlab.like2(255))

##### How do we map change related to event?

### Lower NIR often correlates to areas with high water fraction or inundated:

#1) Simple thresholding
#2) Difference 
#3) PCA

### 1) Simple thresholding

r_rec_NDWI_before <- subset(r_before,"NIR") < 0.2
r_rec_NDWI_after <- subset(r_after,"NIR") < 0.2

### THis is suggesting flooding!!!
plot(r_rec_NDWI_before, main="Before RITA, NDWI > 0.2")
plot(r_rec_NDWI_after, main="After RITA, NDWI > 0.2")

#Compare to actual flooding data
freq_fema_zones <- as.data.frame(freq(r_ref))
xtab_threshold <- crosstab(r_ref,r_rec_NIR_after,long=T)

## % overlap between the flooded area and values below 0.2 in NIR
(xtab_threshold[5,3]/freq_fema_zones[2,2])*100 #agreement with FEMA flooded area in %.

### 2) Difference

r_diff <- r_after_NDWI - r_before_NDWI

plot(r_diff)
histogram(r_diff)

r_std <- (r_diff - cellStats(r_diff,"mean"))/cellStats(r_diff,"sd")

plot(r_std)
hist(r_std)
abline(v=1.96,col="red",tly=1)
abline(v=-1.96,col="red",tly=1)

r_change_pos <- r_std > 1.96
r_change_neg <- r_std < -1.96
plot(r_change_pos)
plot(r_change_neg)

### 3) PCA

#Correlate long term mean to PC!
cor_mat_layerstats <- layerStats(stack(r_after_NDWI,r_before_NDWI), 'pearson', na.rm=T)
cor_matrix <- cor_mat_layerstats$`pearson correlation coefficient`
class(cor_matrix)
print(cor_matrix) #note size is 7x7

pca_mod <- principal(cor_matrix,nfactors=2,rotate="none")
class(pca_mod$loadings)
print(pca_mod$loadings)

df_loadings <- pca_mod$loadings
head(df_loadings)

### Generate scores from eigenvectors
### Using predict function: this is recommended for raster.
r_pca <- predict(stack(r_after_NDWI,r_before_NDWI), pca_mod, index=1:2,filename="pc_scores.tif",overwrite=T) # fast
plot(r_pca,y=2,zlim=c(-2,2))
plot(r_pca,y=1,zlim=c(-2,2))

histogram(r_pca)
plot(r_pca$pc_scores.2 > 1.96)
plot(r_pca$pc_scores.2 < -1.96)

plot(r_pca$pc_scores.2 > 1.645)
plot(r_pca$pc_scores.2 < -1.645)

##### Plot on the loading space using TCBI and TCWI as supplementary variables
var_labels <- rownames(loadings_df)

plot(loadings_df[,1],loadings_df[,2],
     type="p",
     pch = 20,
     col ="blue",
     xlab=names(loadings_df)[1],
     ylab=names(loadings_df)[2],
     ylim=c(-1,1),
     xlim=c(-1,1),
     axes = FALSE,
     cex.lab = 1.2)

points(cor_r_TC$pc_scores.1[1],cor_r_TC$pc_scores.2[1],col="red")
points(cor_r_TC$pc_scores.1[2],cor_r_TC$pc_scores.2[2],col="green")

axis(1, at=seq(-1,1,0.2),cex=1.2)
axis(2, las=1,at=seq(-1,1,0.2),cex=1.2) # "1' for side=below, the axis is drawned  on the right at location 0 and 1
box()    #This draws a box...

###################################   End of script ###########################################

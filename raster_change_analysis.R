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


###################################   End of script ###########################################
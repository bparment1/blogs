
############################################################################
#####  Parameters and argument set up ###########

#Note: https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20170008750.pdf
out_suffix <- "pca_test_10192018" #output suffix for the files and ouptut folder

in_dir <- "/nfs/bparmentier-data/Data/blogs/blog2_Time_series_Analysis_PCA/data"
out_dir <- "/nfs/bparmentier-data/Data/blogs/blog2_Time_series_Analysis_PCA/data/outputs"

file_format <- ".tif" #PARAM 4

pattern_str <- "ndvi3g_geo.*.201.*.tif" #pattern for hte GIMMS dataset
infile_reg_outline <- "africa_roi.tif" #region of interest

r_ref <- raster(file.path(in_dir,infile_reg_outline))
lf <- list.files(pattern=pattern_str,
                 path=in_dir,
                 full.names = T)
r_stack <- stack(lf)
plot(r_stack,y=1)
plot(r_ref)
r_NDVI_s <- crop(r_stack,r_ref)

r_mask <- r_NDVI_s < -10000

r_NDVI_s <- mask(r_NDVI_s,r_mask,maskvalue=1)
plot(r_NDVI_s,y=1)


#### Generate Color composites
### True color composite

#Correlate long term mean to PC!
cor_mat_layerstats <- layerStats(r_NDVI_s, 'pearson', na.rm=T)
#cov_mat_layerstats <- layerStats(r_NDVI_s, 'cov', na.rm=T)

cor_matrix <- cor_mat_layerstats$`pearson correlation coefficient`
#cov_matrix <- cov_mat_layerstats$`covariance`


class(cor_matrix)
dim(cor_matrix)
#View(cor_matrix)
image(cor_matrix) #visualize the matrix

pca_mod <-principal(cor_matrix,nfactors=3,rotate="none")
#pca_mod <-principal(cov_matrix,nfactors=3,rotate="none")

#pca_mod <-principal(cov_matrix,nfactors=3,rotate="none")

#class(pca_mod$loadings)
#str(pca_mod$loadings)
#plot(pca_mod$loadings[,1],type="b",
#     xlab="time steps",
#     ylab="PC loadings",
#     ylim=c(-1,1),
#     col="blue")
#lines(pca_mod$loadings[,2],type="b",col="red")
#lines(pca_mod$loadings[,3],type="b",col="black")
#title("Loadings for the first three components using T-mode")

##Make this a time series
loadings_df <- as.data.frame(pca_mod$loadings[,1:3])
list_file_names <- (str_split(names(r_NDVI_s),"_"))
list_file_names

extract_dates <- function(filename_val){
  n_string <- length(filename_val)
  n_start <- n_string - 3 + 1
  
  date_tmp <- filename_val[n_start:n_string]
  date_tmp <- paste(date_tmp,collapse=".")
  date_queried <- date_tmp
  #date_queried <- as.Date(date_tmp ,
  #                        format = "%Y.%m.%d")
  return(date_queried)
}

undebug(extract_dates)

dates_val <- extract_dates(list_file_names[[1]])

dates_val <- unlist(lapply(list_file_names,FUN=extract_dates))

dates_val <- as.Date(dates_val ,format = "%Y.%m.%d")

#dates_val <- 1:24
pca_loadings_dz <- zoo(loadings_df,dates_val) #create zoo object from data.frame and date sequence object
#?plot.zoo to find out about zoo time series plotting of indexes
plot(pca_loadings_dz,
     type="b",
     plot.type="single",
     col=c("blue","red","black"),
     xlab="time steps",
     ylab="PC loadings",
     ylim=c(-1,1))
title("Loadings for the first three components using T-mode")
names_vals <- c("pc1","pc2","pc3")
legend("bottomright",legend=names_vals,
       pt.cex=0.8,cex=1.1,col=c("blue","red","black"),
       lty=c(1,1), # set legend symbol as lines
       pch=1, #add circle symbol to line
       lwd=c(1,1),bty="n")

## Add scree plot
#plot(pca_mod$values,main="Scree plot: Variance explained")

###################

plot(Re(fft(pca_mod$loadings[,2])))

#do a AR
acf(pca_mod$loadings[,2])
acf(pca_mod$loadings[,3])
ccf(pca_mod$loadings[,2],pca_mod$loadings[,3])
ccf(pca_mod$loadings[,2],pca_mod$loadings[,2])

pacf(pca_mod$loadings[,2])
# calculate fft of data

test <- fft(pca_mod$loadings[,3])


# extract magnitudes and phases
magn <- Mod(test) # sqrt(Re(test)*Re(test)+Im(test)*Im(test))
phase <- Arg(test) # atan(Im(test)/Re(test))

# select only first half of vectors
magn.1 <- magn[1:(length(magn)/2)]
#phase.1 <- Arg(test)[1:(length(test)/2)]

# plot various vectors

# plot magnitudes as analyses by R
#x11()
plot(magn,type="l")

# plot first half of magnitude vector
#x11()
plot(magn.1,type="l")

# generate x-axis with frequencies
x.axis <- 1:length(magn.1)/time

# plot magnitudes against frequencies
x11()
plot(x=x.axis,y=magn.1,type="l")

spectrum((pca_mod$loadings[,2]))

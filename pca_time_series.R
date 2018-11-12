list_file_names <- (str_split(names(r_NDVI_s),"_"))
list_file_names

extract_dates <- function(filename_val){
  n_string <- length(filename_val)
  n_start <- n_string - 3 + 1
  
  date_tmp <- filename_val[n_start:n_string]
  dates_queried <- paste(date_tmp,collapse=".")
  date_queried <- as.Date(date_tmp ,
                           format = "%Y.%m.%d")
  return(dates_queried)
}

dates_val <- unlist(lapply(list_file_names,FUN=extract_dates))

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
legend("topright",legend=names_vals,
       pt.cex=0.8,cex=1.1,col=c("blue","red","black"),
       lty=c(1,1), # set legend symbol as lines
       pch=1, #add circle symbol to line
       lwd=c(1,1),bty="n")

## Add scree plot
plot(pca_mod$values,main="Scree plot: Variance explained")

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

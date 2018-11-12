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

date_tmp <- list_file_names[[1]][n_start:n_string]
dates_queried <- paste(date_tmp,collapse=".")

month_vals <- month_range[1]:month_range[2]
month_vals <- sprintf("%02d", month_vals)

# print today's date
today <- Sys.Date()
format(today, format="%B %d %Y")
"June 20 2007"

#list_dates <-unlist(lapply(month_vals,function(x){paste0(x,".",c("01","15"))}))
list_dates <-unlist(lapply(month_vals,function(x){paste0(x,"_",c("01","15"))}))

dates_queried <- format(date_tmp,"%Y.%m.%d") #formatting queried dates
#dates_val <- format(ll[-13],"%Y_%m_%d") #formatting queried dates
dates_val <- paste(year_val,list_dates,sep="_") #output dates
names()

plot(pca_mod$loadings[,1],type="b",
     xlab="time steps",
     ylab="PC loadings",
     ylim=c(-1,1),
     col="blue")
lines(pca_mod$loadings[,2],type="b",col="red")
lines(pca_mod$loadings[,3],type="b",col="black")
title("Loadings for the first three components using T-mode")

r_test <- r_NDVI_s < -10000


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

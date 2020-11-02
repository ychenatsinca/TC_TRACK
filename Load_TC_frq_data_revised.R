# load TC return frequency data 
# date : 2020-10-13
#

fun_read_nc <- function(arg1) {
  #load  ncdf library
  #
  library(ncdf4)
  #arg1: filepath fo the nc file from James multilayer output
  print(paste("arg1: for reading file path ;", arg1))
  # open the read in file and copy the variables to the dataframe for analysis
  input_nc <- nc_open(arg1)
  #str(input <- nc)
  #   print(i)
  
  result <- list()
  for (i in 1:input_nc$ndims ) {
    #   print(i)
    # store each  variable with respect of the dim <- name
    result[[input_nc$dim[[i]]$name]] <- input_nc$dim[[i]]$vals
  }
  
  for (i in 2:length(input_nc$var) ) {
    # store each variable with respect of the var <- name
    result[[input_nc$var[[i]]$name]] <- ncvar_get(input_nc,input_nc$var[[i]]$name)
  }
  
  nc_close(input_nc)
  # show result structure
  print(str(result))
  # export the datatable
  return(result)
}  

lai <- fun_read_nc(arg1="/lfs/home/ychen/LAI_STUDY_EAsia/LAI_DATA/c_gls_LAI_201308310000_VGT_V2_EAST_ASIA.nc")
#lai <- fun_read_nc(arg1="H:\\LAI_STUDY_EAsia\\LAI_DATA\\c_gls_LAI_201308310000_VGT_V2_EAST_ASIA.nc")

# get dimensions from grid a for working
#mask_$LANDMASK
nav_lon <- lai$lon
nav_lat <- lai$lat
nx <- length(lai$lon)
ny <- length(lai$lat)


#get land mask 
land.mask <- lai$LAI
land.mask[ is.na(lai$LAI) == FALSE ] <- 1

# load libraries and  TC_track data from the selected year 
library("plyr")
library("fields")
library("rgdal")
library("maptools")
library("raster")

#load the coastlines data by readOGR function from sp package 
coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")


#  Totoal 75 years TC occurance data 


# create idex of year

ld_go <- T
if (ld_go) {

yr.id <- seq(1999, 2018,1) - 1999 + 1
count.1999.to.2018.2d <- array( 0, dim=c(6722,6722)) 
count.1999.to.2018.3d <- array( 0, dim=c(6722,6722)) 
count.1999.to.2018.4d <- array( 0, dim=c(6722,6722)) 
count.1999.to.2018.5d <- array( 0, dim=c(6722,6722)) 
#count.1999.to.2019.7d <- array( 0, dim=c(6722,6722)) 
load("2D_8_60_1999_2018_.rda")
#sum up all arraies in ann.land.frq 
# sum up TC occurence 
for ( it in yr.id ) { 
 print(paste("working on year:",it+1999-1,"for 2D") )
      count.1999.to.2018.2d <-  count.1999.to.2018.2d + tc.occ.avg[,,it]
    }    

#count annual occurance of TC 
tc.ave.occ.2d <- count.1999.to.2018.2d / as.numeric(length(yr.id))
tc.ave.occ.2d[tc.ave.occ.2d<=0.1] <- 0
} #lf_ld_go

#copy information to arr.for.plot 
arr.for.plot <- tc.ave.occ.2d 

#load("TC_return_frq_1945to2019_7D.rda")
#sum up all arraies in ann.land.frq 
# sum up TC occurence 
#for ( it in yr.id ) { 
#print(paste("working on file:",it,"for 7D") )
#      count.1999.to.2019.7d <-  count.1999.to.2019.7d + ann.land.frq[,,it]
#    }    
#count annual occurance of TC 
#tc.ave.occ.7d <- count.1999.to.2019.7d / as.numeric(length(yr.id))
#tc.ave.occ.7d[tc.ave.occ.5d<=0.01] <- 0



#load forest cover type 

# 2. ESA land cover map from 2015, the domain shold be also matched the LAI map over WP region

load(paste("/lfs/home/ychen/LAI_STUDY_EAsia/LANDCOVER_DATA/2015.esa.landcover.east.asia.rda",sep=""))
# variable name: esa.lc
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- NA
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 40 ) ] <- 1
esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3

# create the LC mask
lc.mask <- esa.lc
# agricuture mask
lc.agr.mask <- lc.mask
lc.agr.mask[lc.mask!=1] <- NA
lc.agr.mask[lc.mask==1] <- 1
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- NA
lc.for.mask[lc.mask==2] <- 1
# other mask
lc.oth.mask <- lc.mask
lc.oth.mask[lc.mask!=3] <- NA
lc.oth.mask[lc.mask==3] <- 1



# calculate the latitude band 

lat.band.for.2d <- array(0, dim=c(60) )
lat.band.for.3d <- array(0, dim=c(60) )
lat.band.for.4d <- array(0, dim=c(60) )


#find the array index which the same latitude zone

tc.ave.occ.2d[tc.ave.occ.2d >= 0.1] <- 1
tc.ave.occ.2d[tc.ave.occ.2d < 0.1] <- 0

#tc.ave.occ.3d[tc.ave.occ.3d >= 0.05] <- 1
#tc.ave.occ.3d[tc.ave.occ.3d < 0.05] <- 0

#tc.ave.occ.5d[tc.ave.occ.5d >= 0.05] <- 1
#tc.ave.occ.5d[tc.ave.occ.5d < 0.05] <- 0

#tc.ave.occ.7d[tc.ave.occ.7d >= 0.05] <- 1
#tc.ave.occ.7d[tc.ave.occ.7d < 0.05] <- 0


for ( iy in 1:60 ) {

  y1 <- (iy-1)
  y2 <- (iy)
  y.id <- which( (lai$lat >= y1) &  (lai$lat < y2) )
  
#  lat.band.for[iy] <- sum(lc.for.mask[,y.id],  na.rm=T) 
  lat.band.for.2d[iy] <- sum( lc.for.mask[,y.id]*tc.ave.occ.2d[,y.id], na.rm=T   )
#  lat.band.for.2d[iy] <- sum( lc.for.mask[,y.id]*tc.ave.occ.2d[,y.id], na.rm=T   )
#  lat.band.for.3d[iy] <- sum( lc.for.mask[,y.id]*tc.ave.occ.3d[,y.id], na.rm=T   )
#  lat.band.for.5d[iy] <- sum( lc.for.mask[,y.id]*tc.ave.occ.5d[,y.id], na.rm=T   )
#  lat.band.for.7d[iy] <- sum( lc.for.mask[,y.id]*tc.ave.occ.7d[,y.id], na.rm=T   )
  
}#end for 
 par(mar=c(5,6,5,2))
  
 plot(y=seq(1,60,1), x=lat.band.for.2d, type="l",col="black",lty="dotdash",lwd=1,
      xlab="Potentail Affected area (km^2)", ylab="Latitude", cex.lab=1.5, xlim=c(0,5e+5) )
# lines(y=seq(1,60,1), x=lat.band.for.1d, col="orange",lty="solid",lwd=3)
# lines(y=seq(1,60,1), x=lat.band.for.3d, col="black",lty="dashed",lwd=1)
# lines(y=seq(1,60,1), x=lat.band.for.2d, col="black",lty="dotted",lwd=1)
 
# legend("bottomright", inset=0.08, title="Along Track Diameter",lwd=c(3,1,1,1),
#        c("1D","2D","3D","5D"), lty=c("solid","dotted","dashed","dotdash"),
#        col=c("orange","black","black","black"),horiz=TRUE, cex=0.8)

dev.new()

obj.go <- TRUE
if (obj.go) {
  library(raster)
  #convert array to raster
  
  land.frq.raster <-  raster( x=t(arr.for.plot*lc.for.mask), 
                              xmn=lai$lon[1],  xmx=lai$lon[nx],
                              ymn=lai$lat[ny], ymx=lai$lat[1], 
                              crs=CRS("+proj=longlat +datum=WGS84"))
  par(mar=c(5,6,5,2))
  ##assign the color palette
  #my.color<- rev(colorRampPalette(c(terrain.colors(27),"lightgray" ))(28))
  #my.breaks<- seq(1,100,length.out = 100)
 # my.color<- colorRampPalette(c("lightgray","orange","yellow","forestgreen","forestgreen"))(49)
  my.color<- colorRampPalette(c("lightgray","blue","cyan","green","yellow","red","red"))(15)
  my.breaks<- round(seq(0, 8.0, length.out = 15), digits=1)
  plot(land.frq.raster, ylim=c(0,60),xlim=c(90,150), legend=TRUE,
         col=my.color,breaks=my.breaks, 
       xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.main=2.0,
       main=paste("Annual TC Occurence ",sep="") )
  plot(coastlines, add=T,lwd=2)
  #load track data for specific event 
  #Morakot 2009
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2009/bwp092009.txt",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dashed",lwd=1.5)
  text(y=track.tmp$LatNS[1], x=track.tmp$LonEW[1], "Morakot (2009)",cex=0.8)
  #Megi 2010
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2010/bwp152010.txt",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dashed",lwd=1.5)
  text(y=track.tmp$LatNS[1], x=track.tmp$LonEW[1], "Megi (2010)",cex=0.8)
  #Haiyan 2013
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2013/bwp312013.dat",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dashed",lwd=1.5)
  text(y=5.0, x=140, "Haiyan (2013)",cex=0.8)
  #Meigi 2016
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2016/bwp202016.dat",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dashed",lwd=1.5)
  text(y= 15, x=145,"Meigi (2016)",cex=0.8)


  #plot(track.mask.raster, col = ifelse(!is.na(track.mask.raster), "red", "white"), add=T)
  #plot(track.mask.raster, col = ifelse(!is.na(track.mask.raster), "red", "white"), add=T)
  #plot(coastlines, add=T)
  # add path line by each event
  #for (ieve in 1:nevents) {
  #    #plot the lines
  #    lines ( y=tc.track$LatNS[as.integer(tc.track$CY)==ieve],
  ##	    x=tc.track$LonEW[as.integer(tc.track$CY)==ieve],
  #	    col="gray",lwd=0.1,lty="longdash" )
  #}
  #points(y=tc.track$LatNS,x=tc.track$LonEW,cex=0.2, pch=1,col="red")
  
} # end of obj.go

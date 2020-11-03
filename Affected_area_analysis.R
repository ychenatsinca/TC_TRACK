#
# This script load the track mask, rainfal, wind data, and land cover data to calculate the  affected area based on different criteria  
# Author: Yi-Ying Chen
# First Date: 2020-11-02
# Revised :

fun.tc.mask <- function (run_case=1)
{
# load the input argument for the case of irun   
run_case <- run_case
print( paste("run case id:",run_case) )

# set working dir path
wrk.dir <- c("/lfs/home/ychen/LAI_STUDY_EAsia/")

source("./src_function_ychen_ncdf4.R")
# function for reading the netCDF file
 

#setup the 12 possible cases

runs <- array(0, dim=c(12,3))  
#diameter, maximum wind, accumulative rainfall
runs[1,] <- c("2D",10,60)
runs[2,] <- c("3D",12,80)
runs[3,] <- c("4D",14,100)
runs[4,] <- c("5D",16,120)
#
runs[5,] <- c("2D",10,0)
runs[6,] <- c("3D",12,0)
runs[7,] <- c("4D",14,0)
runs[8,] <- c("5D",16,0)
#
runs[9,]  <- c("2D",0,60)
runs[10,] <- c("3D",0,80)
runs[11,] <- c("4D",0,100)
runs[12,] <- c("5D",0,120)


for (irun in run_case:run_case ) {

# set up the working years
yr.st <- 1999 
yr.nd <- 2018
nyrs <- (yr.nd-yr.st)+1 
wrk.yr <- seq(yr.st, yr.nd, 1)

#decleare TC occurency array
tc.occ.avg <- array(0, dim=c(6722,6722,nyrs))  


# load libraries and  TC_track data from the selected year 
library("plyr")
library("fields")
library("rgdal")
library("maptools")
library("ncdf4")
library("raster")

#load coastlines
coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")

# year loop 
for (iyr in  1:nyrs ) {   

#read landcover data 
# ESA land cover map at 1km, the domain should be matched the LAI map over study area

if (wrk.yr[iyr] > 1999 ) { 
    load(paste(wrk.dir,"/LANDCOVER_DATA/",as.character(wrk.yr[iyr]), ".esa.landcover.east.asia.rda",sep=""))
}else{
    load(paste(wrk.dir,"/LANDCOVER_DATA/1999.esa.landcover.east.asia.rda",sep=""))
}#end if

# variable name: esa.lc  
# group the types to croplands(type 1), forests(type 2), others(type 3)
# see http://http://maps.elie.ucl.ac.be/CCI/viewer/download/ESACCI-LC-Ph2-PUGv2_2.0.pdf   # Page-30 
# set water(210) to NA
esa.lc[ (esa.lc == 210 ) ] <- 0 
esa.lc[ (esa.lc >= 10)  & (esa.lc <= 40 ) ] <- 1
esa.lc[ (esa.lc >= 50)  & (esa.lc <= 120) ] <- 2
#esa.lc[ (esa.lc == 160) | (esa.lc == 170) ] <- 2
esa.lc[ (esa.lc > 2) ]  <- 3 

# create the LC mask
lc.mask <- esa.lc
# forest mask
lc.for.mask <- lc.mask
lc.for.mask[lc.mask!=2] <- 0
lc.for.mask[is.nan(lc.mask)] <- 0
lc.for.mask[lc.mask==2] <- 1

# search lai files/dates for the working year
lai.files.yr <- list.files(path=paste(wrk.dir,"/LAI_DATA/",sep=""), pattern=paste("c_gls_LAI_",as.character(wrk.yr[iyr]),sep="") ) 
wrk.date <- substr(lai.files.yr,start=11,stop=20)
#

# 10 days loop for reading wind speed, rainfall, and TC mask  
for (idate in 1:length(wrk.date) )  {
#
#for (idate in 1:1 )  {

#read tc mask
tc.mask <- fun_read_nc(arg1=paste(wrk.dir,"TRACK_DATA_",runs[irun,1],"/tc_mask_pot_",wrk.date[idate],".nc",sep="") ) 
tc.mask$tc_pot[tc.mask$tc_pot>=1] <- 1
tc.mask$tc_pot[tc.mask$tc_pot!=1] <- 0
tc.mask$tc_pot[is.na(tc.mask$tc_pot)] <- 0

nx <- length(tc.mask$lon) 
ny <- length(tc.mask$lat)
 


tc.pixels <- length(which((tc.mask$tc_pot>=1)))
#check availbe pixel for analysis 
if ( tc.pixels > 1000   ){
  print( paste("Total TC pixeld inthe mask:",tc.pixels,sep=""))
  #read tc 10day maximum wind speed
  tc.mws <- fun_read_nc(arg1=paste(wrk.dir,"/WIND_DATA","/max_wind_",wrk.date[idate],".nc",sep="") ) 
  tc.mws$max_wind[tc.mws$max_wind < as.numeric(runs[irun,2])]  <- 0
  tc.mws$max_wind[tc.mws$max_wind >= as.numeric(runs[irun,2])] <- 1
  tc.mws$max_wind[is.na(tc.mws$max_wind)] <- 0


  #read tc 10day accumulative rainfall
  tc.acf <- fun_read_nc(arg1=paste(wrk.dir,"/PRECIP_DATA","/acc_rainf_",wrk.date[idate],".nc",sep="") )
  tc.acf$acc_rainf[tc.acf$acc_rainf < as.numeric(runs[irun,3])]  <- 0
  tc.acf$acc_rainf[tc.acf$acc_rainf >= as.numeric(runs[irun,3])] <- 1
  tc.acf$acc_rainf[is.na(tc.acf$acc_rainf)] <- 0

  #conditional critiria for determining the affected area 
  tmp.arr <- array(0,dim=c(nx,ny))
  # wind only condition 
  if ( as.numeric(runs[irun,2]) == 0. )  tmp.arr[ (tc.acf$acc_rainf >=1) & (tc.mask$tc_pot >= 1)  ] <- 1  
  # rainfall only
  if ( as.numeric(runs[irun,3]) == 0. )  tmp.arr[ (tc.mws$max_wind >=1) & (tc.mask$tc_pot >= 1)  ] <- 1   
  #combine condition
  if ( (as.numeric(runs[irun,3]) != 0.) & (as.numeric(runs[irun,2]) != 0.) ) {
     tmp.arr[ ((tc.acf$acc_rainf >=1)|(tc.mws$max_wind >=1)) & (tc.mask$tc_pot >= 1)  ] <- 1  
  } 
 #calculate the tc occurance for 10 days
  tc.occ.avg[,,iyr] <- tmp.arr + tc.occ.avg[,,iyr]

} else {
 #eo to the next step for the it loop because of s,all sample size pixels below 20   
  print( "Go to next it/time step! Not enough TC pixels for analysising.. ") 
  next
  # exit it loop
} # end if for checking pixels 

 
} #end of 10days loop 

 tc.occ.avg[,,iyr] <- tc.occ.avg[,,iyr] * lc.for.mask 

 #output the plot for annual TC Occueance
 tc.occ.yr <- array(0,dim=c(nx,ny))
 tc.occ.yr <- raster( x=t(tc.occ.avg[,,iyr]), 
			     xmn=tc.acf$lon[1],  xmx=tc.acf$lon[nx],
			     ymn=tc.acf$lat[ny], ymx=tc.acf$lat[1], 
			     crs=CRS("+proj=longlat +datum=WGS84"))
 
 bitmap(file=paste("./png_files/TC_frq",runs[irun,1],runs[irun,2],runs[irun,3],wrk.yr[iyr],".png",sep="_"),type="png16m",
        width = 4, height = 4, units = "in", res= 300 )
 layout(matrix(data=seq(1,1,1),nrow=1, ncol=1, byrow=TRUE))
 plot(tc.occ.yr,main=paste("Annual TC Occurence:", wrk.yr[iyr],"for runs:",runs[irun,1],runs[irun,2],runs[irun,3],sep="_") )
 plot(coastlines,add=T)
 dev.off()
 # show in window
# plot(tc.occ.yr,main=paste("Annual TC Occurence:", wrk.yr[iyr],"for runs:",runs[irun,1],runs[irun,2],runs[irun,3],sep="_") )
# plot(coastlines,add=T)

 
} #end of  year-loop 

# save the tc.occ.avg array for each run

save(tc.occ.avg, file= paste(runs[irun,1],runs[irun,2],runs[irun,3],yr.st,yr.nd,".rda",sep="_"))  


} #end of run-case-loop


}#end of function 
#========== Set the script to Auto RUN=========================== 
# Start convert the file by fun_compress with specific argunment
# If you want to submit the file by queue to obelix, you need to apply the 
# following lines, which allowed this R-script can be called under the shell script
# with the arggunmets sending by specific batch jobs  
#
args<-commandArgs(TRUE)
print(args)
fun.tc.mask(args[1]) 
#
#========== End ================================================




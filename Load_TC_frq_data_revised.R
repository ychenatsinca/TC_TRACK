# load TC return frequency data 
# first date : 2020-10-13
# revised: 2020-11-06

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
#coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_110m_coastline/ne_110m_coastline.shp")
coastlines <- readOGR("/lfs/home/ychen/GIS/Coastline/ne_10m_coastline/ne_10m_coastline.shp")

# load forest cover type 
# ESA land cover map from 2015, the domain shold be also matched the LAI map over WP region
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

#load the TC occurence
yr.id <- seq(1999, 2018,1) - 1999 + 1

run.name <- c("4D_12_100_1999_2018_.rda", "4D_0_100_1999_2018_.rda", "4D_12_0_1999_2018_.rda",
              "3D_10_80_1999_2018_.rda", "3D_0_80_1999_2018_.rda", "3D_10_0_1999_2018_.rda",
              "2D_8_60_1999_2018_.rda", "2D_0_60_1999_2018_.rda", "2D_8_0_1999_2018_.rda")
              
 led.txt <- c("4D-co","4D-rf","4D-w",
              "3D-co","3D-rf","3D-w",
              "2D-co","2D-rf","2D-w")
               
run.col <- c("forestgreen","forestgreen","forestgreen",
            "blue","blue","blue",
            "black","black","black")

run.lty <- rep(c(1,5,3),3)  
#
nruns = 9
nyears = 20
nlat = 60
#load the TC Rda data
ld_go <- F 
if (ld_go) {
# calculate the latitude band 
lat.band.runs <- array(0, dim=c(nlat,nruns) )
tot.aff.runs <- array(0, dim=c(nyears,nruns) ) 

tot.aff.runs.1yr <- array(0, dim=c(nruns) )
tot.aff.runs.2yr <- array(0, dim=c(nruns) )
tot.aff.runs.5yr <- array(0, dim=c(nruns) )


for (irun in 1:nruns) {
#load the data for iruns
load( file=run.name[irun] )
#decleare array for counting the occurence
count.1999.to.2018 <- array( 0., dim=c(6722,6722)) 
tmp.arr <- array(0, dim=c(6722,6722))
#sum up all arraies in ann.land.frq 
#sum up TC occurence 
for ( it in 1:nyears ) { 
      print(paste("working on year:",it+1999-1) )
      #tc.occ.avg[,,it][ tc.occ.avg[,,it]>=1] <- 1.0
 
      count.1999.to.2018 <-  count.1999.to.2018 + tc.occ.avg[,,it]
    
      #calculate the total area
      #skip longitude <100E
      tmp.arr <- tc.occ.avg[,,it]*lc.for.mask
      x.id1 = 1  #90E
      x.id2 = 1200 #100E
      tmp.arr [x.id1:x.id2,] <- 0
      #skip logitude <110E & latitude > 30
      #x.id1 =1; x.id2 = 1200*1.5
      #y.id1 =1; y.id2 = 1200*3.0 
      #tmp.arr[x.id1:x.id2,y.id1:y.id2]  <- 0
      #count numbers

      #all events 
      tmp.arr[tmp.arr >= 0.1] <- 1
      tmp.arr[tmp.arr < 0.1] <- 0
      #plot(tmp.arr)
      #
      tot.aff.runs[it,irun] <- sum(tmp.arr,na.rm=T) 


}    

#count annual occurance of TC (divided by average event numbers from 1999 to 2018, 29 )  
tc.ave.occ <- count.1999.to.2018 / 20

#filter out 10 yrs events 
tc.ave.occ[tc.ave.occ<0.1] <- 0

#skip longitude <100E
x.id1 = 1  #90E
x.id2 = 1200 #100E
tc.ave.occ[x.id1:x.id2,] <- 0
#skip logitude <110E & latitude > 30
#x.id1 =1; x.id2 = 1200*1.5
#y.id1 =1; y.id2 = 1200*3.0
#tc.ave.occ[x.id1:x.id2,y.id1:y.id2] <- 0
  
# calculate the return frequency 
tc.occ.1yr <- tc.ave.occ
tc.occ.1yr[tc.occ.1yr <  1] <- NA
tc.occ.1yr[tc.occ.1yr >= 1] <- 1.0 
image(tc.occ.1yr[,6722:1])

# pixels with 2yr return frequency
tc.occ.2yr <- tc.ave.occ  
tc.occ.2yr[tc.occ.2yr <  0.5] <- NA
tc.occ.2yr[tc.occ.2yr >= 1.0] <- NA
tc.occ.2yr[tc.occ.2yr >= 0.5] <- 1.0
image(tc.occ.2yr[,6722:1])

# pixels with 5yr return frequency
tc.occ.5yr <- tc.ave.occ  
tc.occ.5yr[tc.occ.5yr >= 0.2] <- NA
tc.occ.5yr[tc.occ.5yr == 0.0] <- NA
tc.occ.5yr[((tc.occ.5yr <=  0.2) & (tc.occ.5yr >  0.0)) ] <- 1.0
image.plot(tc.occ.5yr[,6722:1],legend=T)

#calculate the pixels for difference returen frequency
tot.aff.runs.1yr[irun] <- sum(tc.occ.1yr,na.rm=T)
tot.aff.runs.2yr[irun] <- sum(tc.occ.2yr,na.rm=T)
tot.aff.runs.5yr[irun] <- sum(tc.occ.5yr,na.rm=T)


print(paste("1yr.occ:",tot.aff.runs.1yr[irun]/10000.,
            "; 2yr.occ:",tot.aff.runs.2yr[irun]/10000.,
            "; 5yr.occ:",tot.aff.runs.5yr[irun]/10000.,
            "Unit in [Mha2]",sep="")) 


#copy information to arr.for.plot and accumulate
arr.for.plot <- tc.ave.occ + arr.for.plot 


#find the array index which the same latitude zone
tc.ave.occ[tc.ave.occ > 0.1] <- 1
tc.ave.occ[tc.ave.occ <= 0.1] <- 0

# calculate the latitude zonal band 
for ( iy in 1:nlat ) {

  y1 <- (iy)
  y2 <- (iy)+1
  y.id <- which( (lai$lat >= y1) &  (lai$lat < y2) )
  #  lat.band calculation for irun: 
  lat.band.runs[iy,irun] <- sum( lc.for.mask[,y.id]*tc.ave.occ[,y.id], na.rm=T   )
  
}#end of iy 

}#end of iruns

#average the arr.for.plot
arr.for.plot <- arr.for.plot/as.numeric(nruns)
arr.for.plot[arr.for.plot<=0.1] <- 0

}#end of ld_go

run.lwd <- rep( c(1.,1.,1.),3)

obj.go <- T
if (obj.go) {
  pdf(file="TC_Occ_Fig1.pdf",width=8, height=8)
  library(raster)
  #convert array to raster
  # set a layout for spactail TC occurenc (a), latitudal mean (b), time series plot (c). 
  #layoutsetting for combine plot 
 
  arr.for.plot[arr.for.plot>=4.5] <- 4.9
  land.frq.raster <-  raster( x=t(arr.for.plot*lc.for.mask), 
                              xmn=lai$lon[1],  xmx=lai$lon[nx],
                              ymn=lai$lat[ny], ymx=lai$lat[1], 
                              crs=CRS("+proj=longlat +datum=WGS84"))
  par(oma = c(0.1, 0.1, 0.1, 0.1),mgp=c(2.5,1,0))
# par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
  #layout( matrix(c(3,2,1,1), nrow=2, ncol=2, byrow = TRUE), width=c(1,1),height=c(1,1) )
  #op = par(no.readonly = TRUE)

#subplot (a) left: spatial distribution map of TC Occ.
  my.color<- colorRampPalette(c("gray","cyan","blue","forestgreen","yellow","orange","red","red"))(12)
#  my.breaks<- round(seq(0, 6.0, length.out = 13), digits=1)
  my.breaks<- c(-.5,0.05,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0)
  my.breaks.txt<- c("0","0.05","0.5","1.0","1.5","2.0","2.5","3.0","3.5","4.0","4.5",">5")
  # set subplot margin
  par(fig=c(0.0,0.7,0.2,1.0), mar = c(4, 4.5, 1.5 ,0.1), usr=c(100,150,0,60)) #no buffer in x-y axis 
  plot(coastlines,col="white",lwd=0,ylim=c(0,60),xlim=c(100,150))
  par(fig=c(0.0,0.7,0.2,1.0), mar = c(4, 4.5, 1.5 ,0.1), usr=c(100,150,0,60)) #no buffer in x-y axis 
  plot(land.frq.raster,  legend=FALSE,ylim=c(0,60),xlim=c(90,150),
         col=my.color,breaks=my.breaks, box=T, 
       xlab="Longitude", ylab="Latitude", cex.lab=1., cex.main=2.0,
       main=paste("",sep=""),add=T )
  #plot legend 
  par(fig=c(0.0,0.7,0.2,1.0), mar = c(4, 4.5, 1.5 ,0.1), usr=c(100,150,0,60)) #no buffer in x-y axis 
  plot(land.frq.raster, legend.only=TRUE,ylim=c(0,60),xlim=c(100,150), col=my.color,breaks=my.breaks, 
       smallplot=c(0.22,.23, 0.63, 0.88), add=T,
       axis.args=list(at=my.breaks, labels=my.breaks.txt, cex.axis=0.8),
       legend.args=list(text= expression("TC Occurrence (" * yr^-1 * ")") ,side=2, line=.1, cex=0.8) )
  #add coastaline
  par(fig=c(0.0,0.7,0.2,1.0), mar = c(4, 4.5, 1.5 ,0.1), usr=c(100,150,0,60), xpd=FALSE)  
  plot(coastlines, lwd=0.8,ylim=c(0,60),xlim=c(100,150),add=T  )
  text(x=103,y=57, label="a",cex=1.2, font=2)  
  box()
  axis(side = 2, at = seq(0,60,length.out=7), 
       labels = c("0","10","20","30","40","50","60"), tck = -0.02)
  mtext("Latitude" , side=2, line=2, cex=1.0) 

  axis(side = 1, at = seq(90,150,length.out=7), 
       labels = c("","100","110","120","130","140","150"), tck = -0.02)
  mtext("Longitude" , side=1, line=2, cex=1.0) 

  # 
  ld_go <- F
  if (ld_go) {
  for (iyr in 1999:2019) {
  print(paste("Year:",iyr,sep=""))
  #show tc track for a specific year
  track.files <-  list.files(paste("/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp",iyr,"/",sep=""))
  for ( id in 1:length(track.files) ) {
  track.tmp <- read.csv(file= paste("/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp",iyr,"/",track.files[id] ,sep="") , 
                              sep=",", header=FALSE,colClasses =c( rep("character",28)),fill=TRUE )[,c(1:9,28)]
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX","NAME")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  #show the TC name in the final step
 #
  if ( is.na(track.tmp$LatNS[length(track.tmp$LatNS)])!=T ) { 
    if ( (as.numeric(track.tmp$LatNS[length(track.tmp$LatNS)]) >= 40) ) {
       print( paste(id,":",track.tmp$NAME[length(track.tmp$NAME)-2],">40N",sep="" ))
       lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dotdash",lwd=.5,col="red")
       text(y=track.tmp$LatNS[length(track.tmp$LatNS)]+1,x=track.tmp$LonEW[length(track.tmp$LatNS)]-5,
            labels=paste(id,":",track.tmp$NAME[length(track.tmp$NAME)],sep="") ,cex=0.7,font=2)
    } #end if 
    if ( (as.numeric(track.tmp$LatNS[length(track.tmp$LatNS)]) <= 5) ) {
       print( paste(id,":",track.tmp$NAME[length(track.tmp$NAME)-2],"<5N",sep="" ))
       lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dotdash",lwd=.5,col="blue")
       text(y=track.tmp$LatNS[length(track.tmp$LatNS)]+1,x=track.tmp$LonEW[length(track.tmp$LatNS)]-5,
            labels=paste(id,":",track.tmp$NAME[length(track.tmp$NAME)],sep="") ,cex=0.7,font=2)
   } #end if 
  }#end if 

  }#end for   
  }#end for 

  }#ld_go 
  

  #add load track data for specific event 
  ld_go <- T
  if (ld_go) {
  #Morakot 2009
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2009/bwp092009.txt",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dotdash",lwd=.8)
  text(y=track.tmp$LatNS[1]+4, x=track.tmp$LonEW[1]+5, "Morakot (2009)",cex=0.8)
  #Megi 2010
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2010/bwp152010.txt",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dotdash",lwd=.8)
  text(y=track.tmp$LatNS[1]-2, x=track.tmp$LonEW[1]-2, "Megi (2010)",cex=0.8)
  #Haiyan 2013
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2013/bwp312013.dat",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS[track.tmp$LonEW<=145], x=track.tmp$LonEW[track.tmp$LonEW<=145],pch=21,lty="dotdash",lwd=.8)
  text(y=5.0, x=140, "Haiyan (2013)",cex=0.8)
  #Meigi 2016
  track.tmp <- read.csv(file="/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp2016/bwp202016.dat",
                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns
  colnames(track.tmp) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX")
  # seperate text and numbers
  track.tmp$LatNS <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LatNS) )/10.
  track.tmp$LonEW <- as.numeric( gsub("[^[:digit:]]","",track.tmp$LonEW) )/10.
  lines(y=track.tmp$LatNS, x=track.tmp$LonEW,pch=21,lty="dotdash",lwd=.8)
  text(y= 17, x=140,"Meigi (2016)",cex=0.8)
  }

#plot (a) right: latidual mean 
  par(fig=c(0.72,1.0,0.3,1.0),new=T)
  par(mar = c(0.02, 0.2, 1.5 , 1),mgp=c(3,1,0), xaxs = "i", yaxs = "i")
  plot(y=seq(1,60,1), x=lat.band.runs[,1], type="l",col="white",lty="solid",lwd=1.,
      xlab="", ylab="",bty="n",
      xaxt = "n", yaxt="n", xlim=c(0,2.0e+5),ylim=c(0,60))
  #add axis
  axis(side = 2, at = seq(0,60,length.out=7), labels = FALSE, tck = -0.05)
  par(mgp=c(0,1.0,0))
  axis(side = 1, at = seq(0,2.0e+5,length.out=5),cex=0.5, 
       labels = c("0","5","10","15","20"), tck = -0.05)
  #mtext( "Affected Area", side=1, line=-2.5, cex=1.0 ) 
  mtext( expression("Affected Area (" * Mha^2 * yr^-1 * deg^-1 * ") "), side=1, line=2.5, cex=1.0) #adj 0 for left alignment ffecteci
  #add gray area from 10 to 40 N 
  # usr <- par('usr')
   #rect(usr[1]+2e+3, 10, usr[2]+2e+3, 40, col="gray95", border=NA ) 
  # add 12 runs 
  for (irun in 1:nruns) {  
      lines(y=seq(1,60,1), x=lat.band.runs[,irun], col=run.col[irun],lty=run.lty[irun],lwd=run.lwd[irun])
  } 
  # legend("right", inset=0.08, title="", bty = "n", bg=NA,
  #      legend=led.txt, lwd=run.lwd, lty=run.lty, col=run.col, horiz=FALSE, cex=0.8)
  #text(x=0.6e+5,y=57, labels="",cex=1.2, font=2)
  box()
 
#plot(b) time series 
  par(fig=c(0.0,1.0,0.0,0.22),new=T)
  par(mar = c(3, 4.5, 0.2 , 1), xaxs = "i", yaxs = "i",  mgp=c(2.5, 1, 0))
  plot(x=seq(1999,2018,1),y=tot.aff.runs[,1],ylim=c(0,4e+6),type="n",cex.lab=1., bty="n",
       xaxt="n",yaxt="n",xlab="",ylab="",xlim=c(1997.5,2020))
 #xlab="Year", ylab=expression("TC Disturbed Area (" * km^2 *")") ) 
  axis(side = 1, at = c("1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2018","2019"), 
      labels = c("1998","2000","2002","2004","2006","2008","2010","2012","2014","2016","2018",""), tck = -0.05)
#  mtext("Year", side=1, line=2, cex=1.) 
  axis(side = 2, at = seq(0,4e+6,length.out=5), 
     labels = c("0","100","200","300","400"), tck = -0.05)
  #mtext( "Affected Area" , side=2, line=2, cex=1.) 
  mtext( expression("Total Area (" * Mha^2 * yr^-1 * ")") , side=2, line=2, cex=1.) #adj 0 for left alignment 
  
  # add bars for shadding the potentail year with events of developing extratrpical systems   
  # usr <- par('usr')
  # x.yr <- c(1999,2002,2004,2010,2013,2015,2016,2018)
  # x.wd <- 0.05 
  # for (iyr in x.yr) {
  #     #plot the bar
  #     rect((iyr-x.wd),usr[3]+0.1e+6,(iyr+x.wd),usr[4]*0.75,col="gray95",border=NA ) 
  # }
 #add more runs 
  for (irun in 1:nruns) {
      lines(x=seq(1999,2018,1), y=tot.aff.runs[,irun], col=run.col[irun],lty=run.lty[irun],lwd=run.lwd[irun])
  }
  text(x=1998.5,y=3.7e+6,labels="b",cex=1.2, font=2) #font=2 for font "face" bold

  dev.off()
} # end of obj.go





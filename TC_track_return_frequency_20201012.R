#
# This script load TC track data and create the 10 days TC affected area  
# Author: Yi-Ying Chen
# First Date: 2019-08-30
# Revised : 2020-10-13 

run.case <- c("1D")
factor <- as.numeric(substr(run.case, star=1, stop=1))  

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


#print(paste("fun_LAI_dyn works on date:", lai$time[it.tc], sep=" "))


# get dimensions from grid a for working
#mask_$LANDMASK
#nav_lon <- lai$lon
#nav_lat <- lai$lat
nx <- length(lai$lon)
ny <- length(lai$lat)

dis2center <- array(0, dim=c(nx, ny))

#100km ~one degree
radius <-  c(1.0)

#get land mask 
land.mask <- lai$LAI
land.mask[ is.na(lai$LAI) == FALSE ] <- 1

#create a inital frequency map only for land mask 
land.frq <- land.mask 
land.frq[ is.na(land.mask == 1) ] <- 0 

tmp.track.mask <- array(0,dim=c(nx,ny))
#pdf("track_mask_1999_2017_100km.pdf",width=30, height=24)
#png("track_mask.png", width = 1024, height = 1024, units = "px")
#layout(matrix(data=seq(1,20),nrow=4, ncol=5, byrow=TRUE))
layout(matrix(data=c(1),nrow=1, ncol=1, byrow=TRUE))


# load libraries and  TC_track data from the selected year 
library("plyr")
library("fields")
library("rgdal")
library("maptools")
library("raster")
#load the coastlines data by readOGR function from sp package 
#coastlines <- readOGR("ne-coastlines-10m\\ne_10m_coastline.shp")
#coastlines <- readOGR("I:\\Rscripts\\TC_tack_mask\\ne-coastlines-10m\\ne_10m_coastline.shp")

#wrk.yr <-seq(1999,1999,1)
wrk.yr <-seq(1999,2019,1)
st.yr <- 1
ed.yr <- 75
  
nyr <- (ed.yr-st.yr+1)

ann.land.frq <- array(0, dim=c(nx,ny,nyr))
#
nevents = 1
for ( kk in st.yr:ed.yr) { 

# create a table for hourly TC track 
hourly.track <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("CY","YYYYMMDDHH", "LatNS" , "LonEW","VMAX") )

    #update lai for the year	
    #lai <- fun_read_nc(arg1=paste("/lfs/home/ychen/LAI_STUDY_EAsia/LAI_DATA/c_gls_LAI_201308310000_VGT_V2_EAST_ASIA.nc",sep=""))
    print(paste("working on year:",wrk.yr[kk],sep=""))
    fname.list <-  list.files(path=paste("/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp",wrk.yr[kk],"/",sep=""), pattern = "*.txt|*.dat" )

    #fname.list <-  list.files(path=paste("H:\\LAI_STUDY_EAsia\\TRACK_DATA\\bwp",wrk.yr[kk],"\\",sep=""), pattern = "*.txt|*.dat" )
    # tc.track as a data.frame, update by wrk.yr[kk]  
    tc.track <- data.frame()


    # for the years after 2001 the TC radius is avariable 
    if ( as.integer(wrk.yr[kk]) < as.integer(2001) ) {   

      for ( i in 1:length(fname.list) ) {
        track.tmp <- read.csv(file=paste("/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp",wrk.yr[kk],"/",fname.list[i],sep="") ,
			      sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns  
#        track.tmp <- read.csv(file=paste("H:\\LAI_STUDY_EAsia\\TRACK_DATA\\bwp",wrk.yr[kk],"\\",fname.list[i],sep="") ,
#                              sep=",", header=FALSE,colClasses = "character" )[,1:9] # only import the 1:9 columns  
        #combine the table
        tc.track <- rbind( tc.track, track.tmp)
      }	#end for 
       tc.track$MSLP <- NA
       tc.track$TY   <- NA
       tc.track$WINDCODE <- NA
       tc.track$RAD   <- 0      
       tc.track$RAD1  <- 0      
  
    }else {
      # years after 2001 
    
      for ( i in 1:length(fname.list) ) {
       print(paste("reading file:",fname.list[i],sep=""))
      track.tmp <- read.csv(file=paste("/lfs/home/ychen/LAI_STUDY_EAsia/TRACK_DATA/bwp",wrk.yr[kk],"/",fname.list[i],sep="") ,
			      sep=",", header=FALSE,colClasses = "character" )[,1:14] # only import the 1:17 columns  
#      track.tmp <- read.csv(file=paste("H:\\LAI_STUDY_EAsia\\TRACK_DATA\\bwp",wrk.yr[kk],"\\",fname.list[i],sep="") ,
#                            sep=",", header=FALSE,colClasses = "character" )[,1:14] # only import the 1:17 columns  
      #combine the table
        tc.track <- rbind( tc.track, track.tmp)
      }	#end for 
 
  
    }# end if

     # assign the column names
      colnames(tc.track) <- c("BASIN","CY","YYYYMMDDHH","TECHNUM","TECH","TAU","LatNS","LonEW","VMAX","MSLP","TY","RAD","WINDCODE","RAD1") 
      # seperate text and numbers 
      tc.track$LatNS <- as.numeric( gsub("[^[:digit:]]","",tc.track$LatNS) )/10.
      tc.track$LonEW <- as.numeric( gsub("[^[:digit:]]","",tc.track$LonEW) )/10.
      # covert the unit from nmi to degree (100 km)  1.0 nmi(nutical mile) = 1.852 km
      tc.track$RAD1 <- as.numeric(tc.track$RAD1)*1.852/100. * factor
      # a minimun 100km radius has been assumped
      tc.track$RAD1[tc.track$RAD1 == 0] <- 1.0 * factor
      # subset the dataframe only for the WP basin
      tc.track <- subset(tc.track, BASIN =="WP")

      print(tc.track)

      # copy tarack information from land.mask
        track.mask <- array(0, dim=c(nx,ny) )  
        track.mask <- track.mask*land.mask
      # make a TC mask by the assumption of 50 km  diameter TC track data
      # six hourly step loop (k)
	for ( k in 1:length(tc.track$BASIN) ) { 
         #for ( k in 1:1 ) { 
	 # using linear interpolation to generate the hourly track fro six-hourly track data  
	 if ( k > 1 ) {
          n.stp = 6 
	 # do the linear interpolation from two nearby points 
	   for (iter in 1:n.stp) { 	   
	
	       del.x <- (tc.track$LonEW[k] - tc.track$LonEW[k-1])
               del.y <- (tc.track$LatNS[k] - tc.track$LatNS[k-1])

               # find TC center location in the map 
               x.center <- tc.track$LonEW[k-1]+ as.numeric(del.x/n.stp)*(iter-1)
               y.center <- tc.track$LatNS[k-1]+ as.numeric(del.y/n.stp)*(iter-1)
               
               # write the track information to the hourly table
               temp.record <- data.frame (CY=tc.track$CY[k], YYYYMMDDHH=as.integer(tc.track$YYYYMMDDHH[k])+(iter-1), 
	                                  LatNS=formatC(y.center,digits=3,format="f",width=7),
					  LonEW=formatC(x.center,digits=3,format="f",width=7),
					  VMAX= tc.track$VMAX[k], RAD=tc.track$RAD1[k])
              hourly.track <- rbind( hourly.track, temp.record )
             
              # update the information of TC radius 
               radius.now <- as.numeric(tc.track$RAD1[k])
               #print(paste("TC radius:",radius.now*100,"km",sep=""))
               # find west/upper and east/lower limits of an retangular domain  
               x1 <- x.center - radius.now 
               x2 <- x.center + radius.now
         
               y1 <- y.center - radius.now
               y2 <- y.center + radius.now

               # find the TC center lon, lat index with a fixed diameter 
               x.id <- which( (lai$lon > x1) &  (lai$lon < x2) )
               y.id <- which( (lai$lat > y1) &  (lai$lat < y2) )			 

              # find the distance to the center in the square area  
                   for (jj in y.id ) { 
	           # update the new track area  
	           dis2center[x.id,jj] <- sqrt( (lai$lon[x.id] - x.center)**2. + (lai$lat[jj] - y.center)**2. )
                   # if ( dis2center[x.id,jj]  < radius) track.mask[x.id,jj] <- 1                
                   #track.mask[x.id[which(dis2center[x.id,jj]<radius.now)],jj] <- 1
                   tmp.track.mask[x.id[which(dis2center[x.id,jj]<radius.now)],jj] <- 1
                   # count the accumulative return freq 

                   #land.frq[ x.id[which(dis2center[x.id,jj]<radius.now)],jj ] <-  land.frq[ x.id[which(dis2center[x.id,jj]<radius.now)],jj ] + 
	           #                                                           track.mask[x.id[which(dis2center[x.id,jj]<radius.now)],jj] 
		   } #end for jj
            
	   }# end for iter
           
           new.eve <- as.integer(tc.track$CY[k]) - as.integer(tc.track$CY[k-1])  
           if ( new.eve == 1 )  {
      	      nevents = nevents + 1
    	      #print(as.integer(tc.track$CY[k]   ))
              #print(as.integer(tc.track$CY[k-1] ))
	      #print(new.eve) 
	      print(paste("Event:", tc.track$CY[k], 
	      		  "Date:", tc.track$YYYYMMDDHH[k],
	      		  "VMAX:", tc.track$VMAX[k],
                          "RAD:",tc.track$RAD1[k], sep=""))
              # count the occurance of TC
              land.frq =  land.frq + tmp.track.mask
              track.mask = track.mask + tmp.track.mask
              # loop for next iteration for the new event
              # renew event track mask 
              tmp.track.mask <- array(0, dim=c(nx,ny) )  
              
            } #endif 

            }# end for k-if	   

    ld_table <- TRUE 
    if (ld_table) {
       # write out the hourly track data for a txt file
       write.table(hourly.track, file = paste("./wpb_files/wpb_",wrk.yr[kk],"_hourly.txt", sep = ""), 
                   sep = ",", dec = ".", row.names = FALSE, col.names = TRUE, quote = FALSE )
    }



#create date.id for the analyisi year reference to the LAI 10days product names to hourly step

#date.id=substr(x=list.files(path="H:\\LAI_STUDY_EAsia\\LAI_DATA\\",pattern=paste("*",wrk.yr[kk],"*",sep="")),start=11,stop=20)
date.id=substr(x=list.files(path="/lfs/home/ychen/LAI_STUDY_EAsia/LAI_DATA/",pattern=paste("*",wrk.yr[kk],"*",sep="")),start=11,stop=20)




}# end of kk-loop within a year 

#update annual land.frq to ann.land.frq array
ann.land.frq[,,kk] <- track.mask

#track.mask <- track.mask*land.mask
#lai$LAI[track.mask==1] <- 10
# plot the figure
#image.plot(lai$LAI[,ny:1], main=paste("Typhoon Radius:", radius, "in degree, Year:", wrk.yr[kk], sep=""))  
#image.plot(track.mask[,ny:1], add=T )

# get the information of subset dataset  to Taiwan  
#subset.xid <- which( (lai$lon > 119.901785) & (lai$lon < 122.1251))
#subset.yid <- which( (lai$lat > 21.83928) & (lai$lat < 25.4018))
#track.mask.tw <- track.mask[subset.xid,subset.yid]

#image.plot(track.mask.tw[,length(subset.yid):1], main = paste("Taiwan Region", wrk.yr[kk],sep="") )
n.maps <- length(date.id) 
obj.go <- FALSE
if (obj.go) {
library(raster)
#convert array to raster
#track.mask.raster <- raster( x=t(track.mask), 
#			     xmn=lai$lon[1],  xmx=lai$lon[nx],
#			     ymn=lai$lat[ny], ymx=lai$lat[1], 
#			     crs=CRS("+proj=longlat +datum=WGS84"))

land.frq.raster <-  raster( x=t(ann.land.frq[,,kk]*land.mask), 
			     xmn=lai$lon[1],  xmx=lai$lon[nx],
			     ymn=lai$lat[ny], ymx=lai$lat[1], 
			     crs=CRS("+proj=longlat +datum=WGS84"))

#lai.raster <- raster( x=t(lai$LAI), 
#			     xmn=lai$lon[1],  xmx=lai$lon[nx],
#			     ymn=lai$lat[ny], ymx=lai$lat[1], 
#			     crs=CRS("+proj=longlat +datum=WGS84"))
par(mar=c(5,8,5,0))
##assign the color palette
#my.color<- colorRampPalette(c(terrain.colors(8),"darkgray","darkgray"))(100)
#my.breaks<- seq(1,100,length.out = 100)

plot(land.frq.raster, ylim=c(0,60),xlim=c(90,150), legend=TRUE,
 #    col=my.color,breaks=my.breaks, 
     xlab="Longitude", ylab="Latitude", cex.lab=1.5, cex.main=2.0,
     main=paste(wrk.yr[kk]," TC occurance",sep="") )
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
#
#save( file= paste("/work/ychen/scripts/R/Rscripts/SPOT_VGT/west_pacific/",wrk.yr[kk],"_TC_mask_",radius*100,"km.rda",sep=""), track.mask) 

#dev.off()
} # end of kk in year loop  


go.control <- FALSE

if(go.control) {
#add ledgen for color bar, TC track, affect area 
my.color<- colorRampPalette(c(terrain.colors(12)))(100)
my.breaks<- seq(1, 10,length.out = 100)
plot.new()
# plot region [0:1; 0:1]
# color bar
plot(lai.raster, legend.only=TRUE, col=my.color, 
     legend.width=4, legend.shrink=0.75,horizontal=TRUE,
     smallplot=c(0.2,0.80, 0.7,0.9),
     legend.args=list(text='Background LAI (m2/m2)', side=3, font=2, line=1.0, cex=1.5))

par(mar = par("mar"))
# line
line.x<-seq(0.1,0.5,0.1)
line.y<-rep(0.5,5)
lines(x=line.x, y=line.y, col="gray", lwd=3,lty="dashed")
points(x=line.x, y=line.y, pch=1, col="red",cex=2)
text(x=0.6,y=0.5, label="TC Tracks", cex=2.5, adj = c(0,0) )
#plot circle
symbols (0.3, 0.2, circles=0.01, bg="gray", add=T)
text(x=0.55,y=0.2, label="Affected Area",cex=2.5,adj = c(0,0))
}
#dev.off()
# export the TC affect area return frequency to R data 
# save( file= paste("I:\\Rscripts\\TC_tack_mask\\TC_return_frq_",wrk.yr[st.yr],"to",wrk.yr[ed.yr],".rda",sep=""), ann.land.frq )

#save( file= paste("./TC_return_frq_",wrk.yr[st.yr],"to",wrk.yr[ed.yr],"_",run.case,".rda",sep=""), ann.land.frq )
# dev.off()

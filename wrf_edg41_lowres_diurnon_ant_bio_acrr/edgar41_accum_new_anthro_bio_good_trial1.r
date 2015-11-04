############################################################################
# Program to prepare emission input files for WRF-Chem Version 3.4         #    
############################################################################

rootpath<-"/global/u2/y/yashar/hopper/"
rootpath0<-"/scratch/scratchdirs/yashar/"
rpath<-paste(rootpath0,"wrf_in_out_co2tr_trial1/STILT_svn/trunk/stiltR/",sep="")
rootpath1<-"/global/u2/y/yashar/hopper/"
rootpath2<-"/scratch/scratchdirs/yashar/"
rootpath3<-"/scratch/scratchdirs/yashar/wrfoutput_input_chem-3.4_8_09_12/"
path_wrfoldin <- "/scratch/scratchdirs/yashar/wrfinput_real_ghg_temp/"
wrfoutput<-paste(rootpath2,"wrfoutput_input_chem-3.4_8_09_12/",sep="")
path_edgar<-paste(rootpath1,"WRF/WRFV3-Chem-3.4/test/em_real/",sep="")
path_edgar2<-paste(rootpath2,"PREP-CHEM-SRC-1.3.2_zl/bin/",sep="")
path_anthro<-paste(rootpath1,"WRF/WRFV3-Chem-3.4/test/em_real/",sep="")


#rqd functions
source(paste(rpath,"assignr.r",sep=""))
source(paste(rpath,"getr.r",sep=""))
source(paste(rpath,"existsr.r",sep=""))
source(paste(rpath,"lsr.r",sep=""))
source(paste(rpath,"distance.r",sep=""))
source(paste(rpath,"image.plot.fix.r",sep=""))
source(paste(rpath,"image.plot.plt.fix.r",sep=""))
source(paste(rpath,"col.grads.r",sep=""))

go<-graphics.off
pc<-F #no windows
library(RNetCDF)
require("abind",lib.loc="/global/u2/y/yashar/hopper/WRF/WRFV3-Chem-3.4/test/em_real/")
library(abind)

print("hello1")

#read wrfoutput files

  # get the domain (nest, coarse, etc.)
  domain="2" #"2" for nest

  # count through WRF output files
  final <- 27
   for (j in 1:final){
    ifelse(j <= 9, final<-paste("0",j, sep=""), final <- j)  
  date <- paste("02-",final,sep="")
  
  # get coordinates from WRF output files
   
  filename1 <- paste(wrfoutput,"wrfout_d0",domain,"_2006-",date,"_00:00:00",sep="")
  # load wrf outpur grid
  nc       <- open.nc(filename1)
  lat.wrf  <- (var.get.nc(nc,"XLAT"))
  lon.wrf  <- (var.get.nc(nc,"XLONG"))  
  close.nc(nc)

  lat.wrf0 <- t(lat.wrf[,,1])
  lon.wrf0 <- t(lon.wrf[,,1])
  
  print("hello2") 

######################################################################################
######################################################################################
# For new prep_chem_sources and wrf-chem version
######################################################################################
######################################################################################

 ############load the EDGAR data for each hour separetely###############################

  time <- 24
  edgar_tot_ch4 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],1,24)) # 1 file per each hour
  edgar_tot_co2 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],1,24)) # 1 file per each hour
  edgar_tot_co <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],1,24)) # 1 file per each hour

  xlat <-  array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],24))  # MY: Needed to add these lines
  xlong <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],24))


for(i in 1:time){
  ifelse(i <= 9, time<-paste("0",i, sep=""), time <- i)

# # Get the lat/lon values of the EDGAR CO2 source file and wrf output files (see above)
  filename_edgar <- paste(path_edgar2,"EDG41-CO2_ANT_EMISS_DIURN_ON-T-2006-",date,"-",time,"0000-g",domain,".nc",sep="")
  nc <- open.nc(filename_edgar)
  edgar_ch4 <- (var.get.nc(nc,"ch4_edgar"))
  edgar_co2 <- (var.get.nc(nc,"co2_edgar"))
  edgar_co <- (var.get.nc(nc,"co_edgar"))
  close.nc(nc)
  
  edgar_ch4 <- (edgar_ch4)
  edgar_co2 <- (edgar_co2)
  edgar_co <- (edgar_co)
  
  edgar_tot_ch4[,,1,i] <- edgar_ch4[,] 
  edgar_tot_co2[,,1,i] <- edgar_co2[,] 
  edgar_tot_co[,,1,i] <- edgar_co[,] 

  xlat[,,i]<-  lat.wrf0[,]    #MY: Needed to add these lines
  xlong[,,i] <-  lon.wrf0[,]
  }
### get the units from kg/m²s to mol/km²h

edgar_tot_ch4 <- edgar_tot_ch4*3.6e9/(0.012+4.*0.001)
edgar_tot_co2 <- edgar_tot_co2*3.6e9/(0.012+2.*0.016) # Use this when diurnal cycle is turned on in pre-processor
edgar_tot_co <- edgar_tot_co*3.6e9/(0.012+1.*0.016)
print("hello2a")

### get the times from the WRF output file

path_time <- "/global/u2/y/yashar/hopper/WRF/WRFV3-Chem-3.4/test/em_real/"
 filename_time <- paste(wrfoutput,"wrfout_d0",domain,"_2006-",date,"_00:00:00",sep="")
 
  nc <- open.nc(filename_time)
  times <- (var.get.nc(nc,"Times"))
  close.nc(nc)

times2 <- times[0:24]

print("hello2b")


####wite the netcdf file for anthropogenic emissions

nc.out_wrf  <- create.nc(paste(path_antro,"CO2_edg41_ant_diurnon_",date,"_2006_d0",domain,sep=""))
##############################################################
dim.def.nc(nc.out_wrf,"Time",unlim=TRUE)
dim.def.nc(nc.out_wrf,"DateStrLength",19)
dim.def.nc(nc.out_wrf,"south_north",dim(lat.wrf0)[1])
dim.def.nc(nc.out_wrf,"west_east",dim(lat.wrf0)[2])
dim.def.nc(nc.out_wrf,"bottom_top",1)

print("hello2c")

var.def.nc(nc.out_wrf,"E_CH4", "NC_FLOAT",c(3,2,4,0))
# var.def.nc(nc.out_wrf,"E_CH4TST", "NC_FLOAT",c(3,2,4,0))
var.def.nc(nc.out_wrf,"E_CO2", "NC_FLOAT",c(3,2,4,0))
var.def.nc(nc.out_wrf,"E_CO", "NC_FLOAT",c(3,2,4,0))
var.def.nc(nc.out_wrf,"Times","NC_CHAR",c(1,0)) 

print("hello2c1")

var.def.nc(nc.out_wrf,"XLAT", "NC_FLOAT",c(3,2,0))
att.put.nc(nc.out_wrf, "XLAT", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "XLAT", "MemoryOrder", "NC_CHAR", "XY")
att.put.nc(nc.out_wrf, "XLAT", "description", "NC_CHAR", "LATITUDE, SOUTH IS NEGATIVE")
att.put.nc(nc.out_wrf, "XLAT", "units", "NC_CHAR", "degrees_north")
att.put.nc(nc.out_wrf, "XLAT", "stagger", "NC_CHAR", " ")

print("hello2c2")

var.def.nc(nc.out_wrf,"XLONG", "NC_FLOAT",c(3,2,0))
att.put.nc(nc.out_wrf, "XLONG", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "XLONG", "MemoryOrder", "NC_CHAR", "XY")
att.put.nc(nc.out_wrf, "XLONG", "description", "NC_CHAR", "LONGITUDE, WEST IS NEGATIVE")
att.put.nc(nc.out_wrf, "XLONG", "units", "NC_CHAR", "degrees_east")
att.put.nc(nc.out_wrf, "XLONG", "stagger", "NC_CHAR", " ")

var.put.nc(nc.out_wrf, "E_CH4", edgar_tot_ch4)
att.put.nc(nc.out_wrf, "E_CH4", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "E_CH4", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "E_CH4", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "E_CH4", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "E_CH4", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "E_CH4", "units", "NC_CHAR", "mol km^-2 hr^-1")

# var.put.nc(nc.out_wrf, "E_CH4TST", wet1)
# att.put.nc(nc.out_wrf, "E_CH4TST", "FieldType", "NC_INT", 104)
# att.put.nc(nc.out_wrf, "E_CH4TST", "MemoryOrder", "NC_CHAR", "XYZ")
# att.put.nc(nc.out_wrf, "E_CH4TST", "coordinates", "NC_CHAR", "XLONG XLAT")
# att.put.nc(nc.out_wrf, "E_CH4TST", "description", "NC_CHAR", "EMISSIONS")
# att.put.nc(nc.out_wrf, "E_CH4TST", "stagger", "NC_CHAR", " ")
# att.put.nc(nc.out_wrf, "E_CH4TST", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "E_CO2", edgar_tot_co2)    
# Here, we write EDGAR CO2 emissions from the EDGAR CO2 source emissions file (generated by
# the pre-processor) to this E_CO2 array that we've created here for this new netcdf file
# (which is to be directly input to WRF-Chem in the namelist.input file) for anthropogenic
# CO2 emissions.
att.put.nc(nc.out_wrf, "E_CO2", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "E_CO2", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "E_CO2", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "E_CO2", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "E_CO2", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "E_CO2", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "E_CO", edgar_tot_co)
att.put.nc(nc.out_wrf, "E_CO", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "E_CO", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "E_CO", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "E_CO", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "E_CO", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "E_CO", "units", "NC_CHAR", "mol km^-2 hr^-1")

print("hello2d")

# Global attributes

att.put.nc(nc.out_wrf, "NC_GLOBAL", "TITLE", "NC_CHAR", "OUTPUT FROM *             PROGRAM:WRF/CHEM V3.4 MODEL")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "START_DATE", "NC_CHAR", "2006-02-01_00:00:00")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "SIMULATION_START_DATE", "NC_CHAR", "2006-02-01_00:00:00")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "WEST-EAST_GRID_DIMENSION", "NC_FLOAT", 59)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "SOUTH-NORTH_GRID_DIMENSION", "NC_FLOAT", 59)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "BOTTOM-TOP_GRID_DIMENSION", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DX", "NC_FLOAT", 30000)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DY", "NC_FLOAT", 30000)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "GRID_ID", "NC_INT", 1)                                                                            
att.put.nc(nc.out_wrf, "NC_GLOBAL", "PARENT_ID", "NC_INT", 0)   
att.put.nc(nc.out_wrf, "NC_GLOBAL", "I_PARENT_START", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "J_PARENT_START", "NC_INT", 1)  
att.put.nc(nc.out_wrf, "NC_GLOBAL", "PARENT_GRID_RATIO", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DT", "NC_FLOAT", 120)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "CEN_LAT", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "CEN_LON", "NC_FLOAT", -123.64)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "TRUELAT1", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "TRUELAT2", "NC_FLOAT", 81.)  
att.put.nc(nc.out_wrf, "NC_GLOBAL", "MOAD_CEN_LAT", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "STAND_LON", "NC_FLOAT", -123.64)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "POLE_LAT", "NC_FLOAT", 90.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "POLE_LON", "NC_FLOAT", 0.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "GMT", "NC_FLOAT", 0.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "JULYR", "NC_INT", 2006)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "JULDAY", "NC_INT", 32)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "MAP_PROJ", "NC_INT", 1) 


var.put.nc(nc.out_wrf,"Times",times2)

var.put.nc(nc.out_wrf,"XLAT",xlat)    # M.Y.: Needed to add these two lines
var.put.nc(nc.out_wrf,"XLONG",xlong)

close.nc(nc.out_wrf)

print("hello3")

# Here we use matlab scripts to fill in the required wrfinput_d0* arrays

#########################################################################################
# get biogenic emissions together for the ACRR region (in Northern CA).
#########################################################################################



#### VPRM STUFF

 # filenamewrfin <- paste(path_wrfoldin,"wrfinput_d0",domain,sep="")
 filenamewrfin <- paste(path_edgar,"wrfinput_d0",domain,sep="")
  # load wrf outpur grid (wrfinput_d01)
  nc       <- open.nc(filenamewrfin)
  evi  <- (var.get.nc(nc,"EVI"))
  # print(evi)
  evimax  <- (var.get.nc(nc,"EVIMAX"))
  evimin  <- (var.get.nc(nc,"EVIMIN"))
  lswi  <- (var.get.nc(nc,"LSWI"))
  lswimax  <- (var.get.nc(nc,"LSWIMAX"))
  lswimin  <- (var.get.nc(nc,"LSWIMIN"))
  vegfra  <- (var.get.nc(nc,"VEGFRA_VPRM"))
  t2  <- (var.get.nc(nc,"T2"))
  close.nc(nc)

### times
 path_time <- "/global/u2/y/yashar/hopper/WRF/WRFV3-Chem-3.4/test/em_real/"
 filename_time <- paste(wrfoutput,"wrfout_d0",domain,"_2006-",date,"_00:00:00",sep="")
  nc <- open.nc(filename_time)
  times <- (var.get.nc(nc,"Times"))
  # print(times)
  close.nc(nc)
  print("hello4")
# times3 <- times[0]
  times3 <- times[1]
# times3 <- times[0:24]
print(times3)

evi1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
# print(evi1)
evimax1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
evimin1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
lswi1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
lswimax1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
lswimin1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
vegfra_vprm1 <- array(0,cbind(dim(lat.wrf0)[2],dim(lat.wrf0)[1],8,1))
print("hello5")


# for (k in 0:24){
#  nc       <- open.nc(filenamewrfin)
#   evi  <- (var.get.nc(nc,"EVI"))
#   # print(evi)
#   evimax  <- (var.get.nc(nc,"EVIMAX"))
#   evimin  <- (var.get.nc(nc,"EVIMIN"))
#   lswi  <- (var.get.nc(nc,"LSWI"))
#   lswimax  <- (var.get.nc(nc,"LSWIMAX"))
#   lswimin  <- (var.get.nc(nc,"LSWIMIN"))
#   vegfra  <- (var.get.nc(nc,"VEGFRA_VPRM"))
#   t2  <- (var.get.nc(nc,"T2"))
#   close.nc(nc)

# evi <- (evi)
# evimax <- (evimax)
# evimin <- (evimin)
# lswi <- (lswi)
# lswimax <- (lswimax)
# lswimin <- (lswimin)
# vgfra <- (vegfra)

evi1[,,,1] <- evi
# evi1[,,,k] <- evi[,]
print(evi1)
print("hello5a")
evimax1[,,,1] <- evimax
evimin1[,,,1] <- evimin
lswi1[,,,1] <- lswi
lswimax1[,,,1] <- lswimax
lswimin1[,,,1] <- lswimin
vegfra_vprm1[,,,1] <- vegfra
print("hello6")
# }
nc.out_wrf  <- create.nc(paste(path_antro,"CO2_biogenic_",date,"_2006_d0",domain,sep=""))

##############################################################
dim.def.nc(nc.out_wrf,"Time",unlim=TRUE)
dim.def.nc(nc.out_wrf,"DateStrLength",19)
dim.def.nc(nc.out_wrf,"south_north",dim(lat.wrf0)[1])
dim.def.nc(nc.out_wrf,"west_east",dim(lat.wrf0)[2])
dim.def.nc(nc.out_wrf,"one",1)
# dim.def.nc(nc.out_wrf,"bottom_top",1)
dim.def.nc(nc.out_wrf,"vprm_vgcls",8)

print("hello7")

var.def.nc(nc.out_wrf,"EVI", "NC_FLOAT",c(3,2,5,0))
print("hello7a")
var.def.nc(nc.out_wrf,"EVI_MAX", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"EVI_MIN", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"LSWI", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"LSWI_MAX", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"LSWI_MIN", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"VEGFRA_VPRM", "NC_FLOAT",c(3,2,5,0))
var.def.nc(nc.out_wrf,"Times","NC_CHAR",c(1,0)) 
var.def.nc(nc.out_wrf,"XLAT", "NC_FLOAT",c(3,2,0))
var.def.nc(nc.out_wrf,"XLONG", "NC_FLOAT",c(3,2,0))

print("hello8")

var.put.nc(nc.out_wrf,"XLAT",xlat)
att.put.nc(nc.out_wrf, "XLAT", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "XLAT", "MemoryOrder", "NC_CHAR", "XY")
att.put.nc(nc.out_wrf, "XLAT", "description", "NC_CHAR", "LATITUDE, SOUTH IS NEGATIVE")
att.put.nc(nc.out_wrf, "XLAT", "units", "NC_CHAR", "degrees_north")
att.put.nc(nc.out_wrf, "XLAT", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "XLAT", "long_name", "NC_CHAR", "latitude coordinate")
att.put.nc(nc.out_wrf, "XLAT", "standard_name", "NC_CHAR", "latitude")

print("hello9")

var.put.nc(nc.out_wrf,"XLONG",xlong)
att.put.nc(nc.out_wrf, "XLONG", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "XLONG", "MemoryOrder", "NC_CHAR", "XY")
att.put.nc(nc.out_wrf, "XLONG", "description", "NC_CHAR", "LONGITUDE, WEST IS NEGATIVE")
att.put.nc(nc.out_wrf, "XLONG", "units", "NC_CHAR", "degrees_east")                     
att.put.nc(nc.out_wrf, "XLONG", "stagger", "NC_CHAR", " ")                              
att.put.nc(nc.out_wrf, "XLONG", "long_name", "NC_CHAR", "longitude coordinate")         
att.put.nc(nc.out_wrf, "XLONG", "standard_name", "NC_CHAR", "longitude")   

print("hello9a")

# var.put.nc(nc.out_wrf, "CPOOL", cpool1)
# att.put.nc(nc.out_wrf, "CPOOL", "FieldType", "NC_INT", 104)
# att.put.nc(nc.out_wrf, "CPOOL", "MemoryOrder", "NC_CHAR", "XYZ")
# att.put.nc(nc.out_wrf, "CPOOL", "coordinates", "NC_CHAR", "XLONG XLAT")
# att.put.nc(nc.out_wrf, "CPOOL", "description", "NC_CHAR", "EMISSIONS")
# att.put.nc(nc.out_wrf, "CPOOL", "stagger", "NC_CHAR", " ")
# att.put.nc(nc.out_wrf, "CPOOL", "units", "NC_CHAR", "mol km^-2 hr^-1")

# var.put.nc(nc.out_wrf, "T_ANN", tann1)
# att.put.nc(nc.out_wrf, "T_ANN", "FieldType", "NC_INT", 104)
# att.put.nc(nc.out_wrf, "T_ANN", "MemoryOrder", "NC_CHAR", "XYZ")
# att.put.nc(nc.out_wrf, "T_ANN", "coordinates", "NC_CHAR", "XLONG XLAT")
# att.put.nc(nc.out_wrf, "T_ANN", "description", "NC_CHAR", "EMISSIONS")
# att.put.nc(nc.out_wrf, "T_ANN", "stagger", "NC_CHAR", " ")
# att.put.nc(nc.out_wrf, "T_ANN", "units", "NC_CHAR", "mol km^-2 hr^-1")

# var.put.nc(nc.out_wrf, "WETMAP", wetmap1)
# att.put.nc(nc.out_wrf, "WETMAP", "FieldType", "NC_INT", 104)
# att.put.nc(nc.out_wrf, "WETMAP", "MemoryOrder", "NC_CHAR", "XYZ")
# att.put.nc(nc.out_wrf, "WETMAP", "coordinates", "NC_CHAR", "XLONG XLAT")
# att.put.nc(nc.out_wrf, "WETMAP", "description", "NC_CHAR", "EMISSIONS")
# att.put.nc(nc.out_wrf, "WETMAP", "stagger", "NC_CHAR", " ")
# att.put.nc(nc.out_wrf, "WETMAP", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "EVI", evi1)
print("hello9b")
att.put.nc(nc.out_wrf, "EVI", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "EVI", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "EVI", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "EVI", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "EVI", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "EVI", "units", "NC_CHAR", "mol km^-2 hr^-1")

print("hello9c")

var.put.nc(nc.out_wrf, "EVI_MAX", evimax1)
att.put.nc(nc.out_wrf, "EVI_MAX", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "EVI_MAX", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "EVI_MAX", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "EVI_MAX", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "EVI_MAX", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "EVI_MAX", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "EVI_MIN", evimin1)
att.put.nc(nc.out_wrf, "EVI_MIN", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "EVI_MIN", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "EVI_MIN", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "EVI_MIN", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "EVI_MIN", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "EVI_MIN", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "LSWI", lswi1)
att.put.nc(nc.out_wrf, "LSWI", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "LSWI", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "LSWI", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "LSWI", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "LSWI", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "LSWI", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "LSWI_MAX", lswimax1)
att.put.nc(nc.out_wrf, "LSWI_MAX", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "LSWI_MAX", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "LSWI_MAX", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "LSWI_MAX", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "LSWI_MAX", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "LSWI_MAX", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "LSWI_MIN", lswimin1)
att.put.nc(nc.out_wrf, "LSWI_MIN", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "LSWI_MIN", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "LSWI_MIN", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "LSWI_MIN", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "LSWI_MIN", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "LSWI_MIN", "units", "NC_CHAR", "mol km^-2 hr^-1")

var.put.nc(nc.out_wrf, "VEGFRA_VPRM", vegfra_vprm1)
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "FieldType", "NC_INT", 104)
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "MemoryOrder", "NC_CHAR", "XYZ")
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "coordinates", "NC_CHAR", "XLONG XLAT")
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "description", "NC_CHAR", "EMISSIONS")
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "stagger", "NC_CHAR", " ")
att.put.nc(nc.out_wrf, "VEGFRA_VPRM", "units", "NC_CHAR", "mol km^-2 hr^-1")

print("hello10")
var.put.nc(nc.out_wrf,"Times",times3)
var.put.nc(nc.out_wrf,"XLAT",xlat)
var.put.nc(nc.out_wrf,"XLONG",xlong)

att.put.nc(nc.out_wrf, "NC_GLOBAL", "TITLE", "NC_CHAR", "OUTPUT FROM *   PROGRAM:WRF/CHEM V3.4 MODEL")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "START_DATE", "NC_CHAR", "2006-02-01_00:00:00")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "SIMULATION_START_DATE", "NC_CHAR", "2006-02-01_00:00:00")
att.put.nc(nc.out_wrf, "NC_GLOBAL", "WEST-EAST_GRID_DIMENSION", "NC_FLOAT", 59)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "SOUTH-NORTH_GRID_DIMENSION", "NC_FLOAT", 59)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "BOTTOM-TOP_GRID_DIMENSION", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DX", "NC_FLOAT", 30000)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DY", "NC_FLOAT", 30000)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "GRID_ID", "NC_INT", 1)

att.put.nc(nc.out_wrf, "NC_GLOBAL", "PARENT_ID", "NC_INT", 0)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "I_PARENT_START", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "J_PARENT_START", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "PARENT_GRID_RATIO", "NC_INT", 1)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "DT", "NC_FLOAT", 120)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "CEN_LAT", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "CEN_LON", "NC_FLOAT", -123.64)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "TRUELAT1", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "TRUELAT2", "NC_FLOAT", 81.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "MOAD_CEN_LAT", "NC_FLOAT", 39.73)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "STAND_LON", "NC_FLOAT", -123.64)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "POLE_LAT", "NC_FLOAT", 90.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "POLE_LON", "NC_FLOAT", 0.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "GMT", "NC_FLOAT", 0.)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "JULYR", "NC_INT", 2006)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "JULDAY", "NC_INT", 32)
att.put.nc(nc.out_wrf, "NC_GLOBAL", "MAP_PROJ", "NC_INT", 1)



close.nc(nc.out_wrf)
}

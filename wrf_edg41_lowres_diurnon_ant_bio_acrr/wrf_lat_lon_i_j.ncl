    ;************************************************************************************
    ;* NCL script to generate netCDF file containing wind speed at specific location    *
    ;* Example for wrfout_d01_2000-01-24_12:00:00.nc                                    *
    ;*                                        *
    ;************************************************************************************


    ; --------------  LOAD FUNCTIONS AND PROCEDURES ----------------

    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
    load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

    ; --------------  BEGINING OF NCL SCRIPT ----------------


    begin
    ;********************************************************
    ; read in netCDF file and make a loop for all time steps
    ;********************************************************
      in     = addfile("./wrfout_d02_edg41_2006-feb1-5_1way_lowres_diurnon.nc","r")
      ;in     = addfile("./geo_em.d02.nc","r")



    ;************************************************
    ;  - Select lon & lat of the point of interest -
    ;************************************************


     res = True      
     res@returnInt = True               ; False : return real values, True: return interger values
     ; lat = 39.73                  ; Latitude of the point of interest
     ; lon = -123.64                  ; Longitude of the point of interest
    
      ; lat = 33.6
     ; lon = -118.4
    
      lat = 37.9   ; Richmond
     lon = -122.3
    
     ; lat = 40.7  ; SLC
     ; lon = -111.64
     point = wrf_user_ll_to_ij(in,lon,lat,res)       ; wrf_user_ll_to_ij(nc_file,lon,lat,opt)
     
     x = point(1)
     y = point(0)

     print("X location is: " + x)            ; print the value of X at the screen
     print("Y location is: " + y)            ; print the value of Y at the screen


   end

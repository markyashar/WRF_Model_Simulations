; load functions and procedures

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

f=addfile("./wrfout_d02_2006-feb1-27_co2_ant_bio.nc","r")
; Read in variables
ts_in     = f->CO2_ANT(:,:,:,:)
lat2d_in  = f->XLAT(:,:,:)
lon2d_in  = f->XLONG(:,:,:)
times_in  = wrf_user_getvar(f,"times",-1)  ; get all times in the file
times = wrf_user_getvar(f,"times",-1)  ; get all times in the file
; print(ts_in)
; print(lat2d_in)
; print(lon2d_in)
; print(times_in)
; print(times)

; Get dimensions from input variable
dimvar_in = dimsizes(ts_in)
ntim_in   = dimvar_in(0)
vert_in   = dimvar_in(1)

; print(vert_in) 
; print(dimvar_in)
; print(ntim_in)
; print(jlat_in)
; print(ilon_in)

;trim buffer (buffer size is set by the user above)
; ts = ts_in(:,:)
ts = ts_in(:)
ts1= ts_in(:)
lat2d = lat2d_in(:)
lon2d = lon2d_in(:)

; print(ts)
; print(lat2d)
; print(lon2d)

;get dimensions of trimmed variable
dimvar = dimsizes(ts)
ntim = dimvar(0)
; jlat = dimvar(1)
; ilon = dimvar(2)

; print(dimvar)
; print(ntim)
; print(jlat)
; print(ilon)

;use ncl functions to compute simple statistics on the data
;whole domain for all time planes
d_avg = avg(ts)                                         ;domain average over all time planes
d_max = max(ts)                                         ;domain max over all time planes
d_min = min(ts)                                         ;domain min over all time planes
d_stddev = stddev(ts)                                   ;domain standard deviation over all time planes
d_var = variance(ts)                                    ;domain variance over all time planes
d_median = dim_median(ts1)

print("ALL STATS LISTED HERE ARE FOR THE CENTER POINT OF THE NESTED DOMAIN ONLY, (I,J)=(32,44) in ACRR")
print("      ")
print ("min/max Anthropogenic CO2 mixing ratio (ppmv) at domain center (ACRR) over all k planes and over all time planes: " + d_min+ "/" + d_max)
print("Average CO2 anthropogenic mixing ratio (ppmv) at domain center (ACRR) over all k planes and over all time planes: " + " " + d_avg)
print("Standard Deviation of CO2 anthropogenic mixing (ppmv) ratio at domain center (ACRR) over all k planes and over all time planes:" + " " + d_stddev)
print("Variance of CO2 anthropogenic mixing (ppmv) ratio at domain center (ACRR) over all k planes and over all time planes:" + " " + d_var)
print("Median of CO2 anthropogenic mixing (ppmv) ratio at domain center (ACRR) for k=1 and over all time planes:" + " " + d_median)
print("      ")

;Whole domain at each time plane
;Allocate memory for variables. Memory allocation in NCL occurs during a 'new' call, assignment, subscripting, and expression evaluation.
; The 'typeof' function returns the string name of the input type. This is especially useful when you need to create a new variable 
; whose type is the same as an existing variable's type.
co2_avg = new((/ntim/),typeof(ts))
co2_max = new((/ntim/),typeof(ts))
co2_min = new((/ntim/),typeof(ts))
co2_stddev = new((/ntim/),typeof(ts))
co2_var = new((/ntim/),typeof(ts))
co2_value = new((/ntim/),typeof(ts))
co2_median = new((/ntim/),typeof(ts1))
 
; print(t_avg)
; print(t_max)
; print(t_min)
; print(t_stddev)
; print(t_var)
; print(ntim)
do i=0,ntim-1
        co2_avg(i) = avg(ts(i,:))                               ;domain average over k for each time plane
        co2_max(i) = max(ts(i,:))                               ;domain max for each time plane
        co2_min(i) = min(ts(i,:))                               ;domain min for each time plane
        co2_stddev(i) = stddev(ts(i,:))                 ;domain standard deviation for each time plane
        co2_var(i) = variance(ts(i,:))                  ;domain variance for each time plane
        co2_value(i) = ts(i,1)                          ;domain value for k=1 for each time plane 
        co2_median(i) = dim_median(ts(i,1))             ;domain median for k=1 for each time plane
print("Average CO2 anthro mixing ratio (ppmv) at domain center (ACRR) over all vertical planes for" + " " + times_in(i) + ":" + "  " + co2_avg(i))
print("Minimum CO2 anthro mixing ratio (ppmv) at domain center (ACRR) over all vertical planes for" + " " + times_in(i) + ":" + "  " + co2_min(i))
print("Maximum CO2 anthro mixing ratio (ppmv) at domain center (ACRR) over all vertical planes for " + " " + times_in(i) + ":" + "  " + co2_max(i))
print("Standard Deviation of CO2 anthropogenic mixing ratio (ppmv) at domain center (ACRR) over all vertical planes for" + " " + times_in(i) + ":" + "  " + co2_stddev(i))
print("Variance of CO2 anthro mixing ratio (ppmv) at domain center (ACRR) over all vertical planes for" + " " + times_in(i) + ":" + "  " + co2_var(i))
print("Value of CO2 anthro mixing ratio (ppmv) at domain center (ACRR) at k=1 for" + " " + times_in(i) + ":" + "  " + co2_value(i))


end do

end
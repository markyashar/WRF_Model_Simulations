load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/ccm_func.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/diagnostics_cam.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/popRemap.ncl"

begin

system("/bin/rm -f CO2_richmond_ant_diurnon_feb1-5_2006_change_d01.nc")  ; remove any pre-existing file

; open file to write variable to (file should not already exist)
out = addfile("CO2_richmond_ant_diurnon_feb1-5_2006_change_d01.nc","c")

;open output netcdf data file to get data from 
data = addfile("./CO2_edg41_ant_diurnon_feb1-5_2006_d01.nc","r")


; First get the variables we will need, and assume they have associated metadata:
; (a) coordinate variables time, bottom_top,XLAT(south_north) ,XLON (east_west) and 
; (b) attributes.



dims = dimsizes(data->E_CO2)
print(dims)

pme = new((/dims(0),dims(1),dims(2),dims(3)/),float)
dimNames=(/"Time","bottom_top","south_north","west_east"/) ; coordinate variables
dimSizes = (/ -1 , dims(1), dims(2), dims(3) /); coordinate dimensions
dimUnlim = (/ True, False, False, False /) ; define unlimited d
filedimdef(out, dimNames, dimSizes, dimUnlim)


times = wrf_user_getvar(data,"times",-1)  ; get all times in the file
ntimes = dimsizes(times)         ; number of times in the file
print(times) 
printVarSummary(times)
print(ntimes)
; do it=0,ntimes-1
; print("Working on time: " + times(it) )      
E_CO2= wrf_user_getvar(data,"E_CO2",-1)
XLAT=wrf_user_getvar(data,"XLAT",-1)
XLONG=wrf_user_getvar(data,"XLONG",-1)


nlat = dimsizes(XLAT)
nlon = dimsizes(XLONG)
printVarSummary(XLAT)
printVarSummary(XLONG)
Times = wrf_user_getvar(data,"Times",-1)

print("hello1")

var= data->E_CO2

var(:,:,:,:) = 0
var(:,:,23,34)=max(E_CO2(:,0,:,:)) ; parent domain
;; var(:,:,34,23)=max(E_CO2(:,0,:,:)) ; parent domain
; var(:,:,24,44)=max(E_CO2(:,0,:,:)) ; nested domain
;; var(:,:,44,24)=max(E_CO2(:,0,:,:)) ; nested domain
print(max(E_CO2(:,0,:,:)))

; output variables directly
out->E_CO2= var
out->XLAT = XLAT
out->XLONG =XLONG
out->Times = Times

; assign named dimensions
E_CO2!0="Time"
E_CO2!1="bottom_top"
E_CO2!2="south_north"
E_CO2!3="west_east"
E_CO2@FieldType = 104
E_CO2@MemoryOrder = "XYZ"
E_CO2@coordinates="XLAT XLONG"
E_CO2@description="EMISSIONS"
E_CO2@stagger=""
E_CO2@units = "mol km^-2 hr^-1"

XLAT!0="Time"
XLAT!1="south_north"
XLAT!2="west_east"

XLONG!0="Time"
XLONG!1="south_north"
XLONG!2="west_east"

times!0="Time"

fAtt = True ; assign file atributes
fAtt@TITLE = "OUTPUT FROM * PROGRAM:WRF/CHEM V3.4 MODEL"
fAtt@START_DATE = "2006-02-01_00:00:00"
fAtt@SIMULATION_START_DATE = "2006-02-01_00:00:00"
fAtt@WEST_EAST_GRID_DIMENSION = 59
fAtt@SOUTH_NORTH_GRID_DIMENSION = 59
fAtt@BOTTOM_TOP_GRID_DIMENSION = 1
fAtt@DX = 30000
fAtt@DY = 30000
fAtt@GRID_ID = 1
fAtt@PARENT_ID = 0
fAtt@I_PARENT_START = 1
fAtt@J_PARENT_START = 1
fAtt@PARENT_GRID_RATIO = 1
fAtt@DT = 120
fAtt@CEN_LAT = 39.73
fAtt@CEN_LON = -123.64
fAtt@TRUELAT1 = 39.73
fAtt@TRUELAT2 = 81.
fAtt@MOAD_CEN_LAT = 39.73
fAtt@STAND_LON = -123.64
fAtt@POLE_LAT = 90.
fAtt@POLE_LON = 0.
fAtt@GMT = 0.
fAtt@JULYR = 2006
fAtt@JULDAY = 32
fAtt@MAP_PROJ = 1
fileattdef(out, fAtt) ; copy file attributes



end 

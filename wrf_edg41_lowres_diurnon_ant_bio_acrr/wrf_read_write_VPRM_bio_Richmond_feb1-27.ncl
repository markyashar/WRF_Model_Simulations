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

system("/bin/rm -f CO2_richmond_only_vprm_bio_feb1-27_2006_d01.nc")  ; remove any pre-existing file

; open file to write variable to (file should not already exist)
out = addfile("CO2_richmond_only_vprm_bio_feb1-27_2006_d01.nc","c")

;open output netcdf data file to get data from 
data = addfile("./CO2_biogenic_02_01-27_2006_d01.nc","r")


; First get the variables we will need, and assume they have associated metadata:
; (a) coordinate variables time, bottom_top,XLAT(south_north) ,XLON (east_west) and 
; (b) attributes.



dims = dimsizes(data->EVI)
print(dims)

pme = new((/dims(0),dims(1),dims(2),dims(3)/),float)
dimNames=(/"Time","vprm_vgcls","south_north","west_east"/) ; coordinate variables
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
EVI= wrf_user_getvar(data,"EVI",-1)
EVI_MAX= wrf_user_getvar(data,"EVI_MAX",-1)
EVI_MIN= wrf_user_getvar(data,"EVI_MIN",-1)
LSWI= wrf_user_getvar(data,"LSWI",-1)
LSWI_MAX= wrf_user_getvar(data,"LSWI_MAX",-1)
LSWI_MIN= wrf_user_getvar(data,"LSWI_MIN",-1)
VEGFRA_VPRM= wrf_user_getvar(data,"VEGFRA_VPRM",-1)
XLAT=wrf_user_getvar(data,"XLAT",-1)
XLONG=wrf_user_getvar(data,"XLONG",-1)


nlat = dimsizes(XLAT)
nlon = dimsizes(XLONG)
printVarSummary(XLAT)
printVarSummary(XLONG)
Times = wrf_user_getvar(data,"Times",-1)

print("hello1")

var1= data->EVI
var2= data->EVI_MAX
var3= data->EVI_MIN
var4= data->LSWI
var5= data->LSWI_MAX
var6= data->LSWI_MIN
var7= data->VEGFRA_VPRM

var1(:,0,:,:) = 0
var1(:,0,23,34)=max(EVI(:,0,23,34)) ; parent domain
var1(:,1,:,:) = 0
var1(:,1,23,34)=max(EVI(:,1,23,34))
var1(:,2,:,:) = 0
var1(:,2,23,34)=max(EVI(:,2,23,34))
var1(:,3,:,:) = 0
var1(:,3,23,34)=max(EVI(:,3,23,34))
var1(:,4,:,:) = 0
var1(:,4,23,34)=max(EVI(:,4,23,34))
var1(:,5,:,:) = 0
var1(:,5,23,34)=max(EVI(:,5,23,34))
var1(:,6,:,:) = 0
var1(:,6,23,34)=max(EVI(:,6,23,34))
var1(:,7,:,:) = 0
var1(:,7,23,34)=max(EVI(:,7,23,34))

var2(:,0,:,:) = 0
var2(:,0,23,34)=max(EVI_MAX(:,0,:,:)) ; parent domain
var2(:,1,:,:) = 0
var2(:,1,23,34)=max(EVI_MAX(:,1,:,:))
var2(:,2,:,:) = 0
var2(:,2,23,34)=max(EVI_MAX(:,2,:,:))
var2(:,3,:,:) = 0
var2(:,3,23,34)=max(EVI_MAX(:,3,:,:))
var2(:,4,:,:) = 0
var2(:,4,23,34)=max(EVI_MAX(:,4,:,:))
var2(:,5,:,:) = 0
var2(:,5,23,34)=max(EVI_MAX(:,5,:,:))
var2(:,6,:,:) = 0
var2(:,6,23,34)=max(EVI_MAX(:,6,:,:))
var2(:,7,:,:) = 0
var2(:,7,23,34)=max(EVI_MAX(:,7,:,:))

var3(:,0,:,:) = 0
var3(:,0,23,34)=max(EVI_MIN(:,0,:,:)) ; parent domain
var3(:,1,:,:) = 0
var3(:,1,23,34)=max(EVI_MIN(:,1,:,:))
var3(:,2,:,:) = 0
var3(:,2,23,34)=max(EVI_MIN(:,2,:,:))
var3(:,3,:,:) = 0
var3(:,3,23,34)=max(EVI_MIN(:,3,:,:))
var3(:,4,:,:) = 0
var3(:,4,23,34)=max(EVI_MIN(:,4,:,:))
var3(:,5,:,:) = 0
var3(:,5,23,34)=max(EVI_MIN(:,5,:,:))
var3(:,6,:,:) = 0
var3(:,6,23,34)=max(EVI_MIN(:,6,:,:))
var3(:,7,:,:) = 0
var3(:,7,23,34)=max(EVI_MIN(:,7,:,:))

var4(:,0,:,:) = 0
var4(:,0,23,34)=max(LSWI(:,0,:,:)) ; parent domain
var4(:,1,:,:) = 0
var4(:,1,23,34)=max(LSWI(:,1,:,:))
var4(:,2,:,:) = 0
var4(:,2,23,34)=max(LSWI(:,2,:,:))
var4(:,3,:,:) = 0
var4(:,3,23,34)=max(LSWI(:,3,:,:))
var4(:,4,:,:) = 0
var4(:,4,23,34)=max(LSWI(:,4,:,:))
var4(:,5,:,:) = 0
var4(:,5,23,34)=max(LSWI(:,5,:,:))
var4(:,6,:,:) = 0
var4(:,6,23,34)=max(LSWI(:,6,:,:))
var4(:,7,:,:) = 0
var4(:,7,23,34)=max(LSWI(:,7,:,:))

var5(:,0,:,:) = 0
var5(:,0,23,34)=max(LSWI_MAX(:,0,:,:)) ; parent domain
var5(:,1,:,:) = 0
var5(:,1,23,34)=max(LSWI_MAX(:,1,:,:))
var5(:,2,:,:) = 0
var5(:,2,23,34)=max(LSWI_MAX(:,2,:,:))
var5(:,3,:,:) = 0
var5(:,3,23,34)=max(LSWI_MAX(:,3,:,:))
var5(:,4,:,:) = 0
var5(:,4,23,34)=max(LSWI_MAX(:,4,:,:))
var5(:,5,:,:) = 0
var5(:,5,23,34)=max(LSWI_MAX(:,5,:,:))
var5(:,6,:,:) = 0
var5(:,6,23,34)=max(LSWI_MAX(:,6,:,:))
var5(:,7,:,:) = 0
var5(:,7,23,34)=max(LSWI_MAX(:,7,:,:))

var6(:,0,:,:) = 0
var6(:,0,23,34)=max(LSWI_MIN(:,0,:,:)) ; parent domain
var6(:,1,:,:) = 0
var6(:,1,23,34)=max(LSWI_MIN(:,1,:,:))
var6(:,2,:,:) = 0
var6(:,2,23,34)=max(LSWI_MIN(:,2,:,:))
var6(:,3,:,:) = 0
var6(:,3,23,34)=max(LSWI_MIN(:,3,:,:))
var6(:,4,:,:) = 0
var6(:,4,23,34)=max(LSWI_MIN(:,4,:,:))
var6(:,5,:,:) = 0
var6(:,5,23,34)=max(LSWI_MIN(:,5,:,:))
var6(:,6,:,:) = 0
var6(:,6,23,34)=max(LSWI_MIN(:,6,:,:))
var6(:,7,:,:) = 0
var6(:,7,23,34)=max(LSWI_MIN(:,7,:,:))

var7(:,0,:,:) = 0
var7(:,0,23,34)=max(VEGFRA_VPRM(:,0,:,:)) ; parent domain
var7(:,1,:,:) = 0
var7(:,1,23,34)=max(VEGFRA_VPRM(:,1,:,:))
var7(:,2,:,:) = 0
var7(:,2,23,34)=max(VEGFRA_VPRM(:,2,:,:))
var7(:,3,:,:) = 0
var7(:,3,23,34)=max(VEGFRA_VPRM(:,3,:,:))
var7(:,4,:,:) = 0
var7(:,4,23,34)=max(VEGFRA_VPRM(:,4,:,:))
var7(:,5,:,:) = 0
var7(:,5,23,34)=max(VEGFRA_VPRM(:,5,:,:))
var7(:,6,:,:) = 0
var7(:,6,23,34)=max(VEGFRA_VPRM(:,6,:,:))
var7(:,7,:,:) = 0
var7(:,7,23,34)=max(VEGFRA_VPRM(:,7,:,:))



;var1(:,:,:,:) = 0
;var1(:,0,24,44)=max(EVI(:,0,:,:)) ; nested domain
;var1(:,1,24,44)=max(EVI(:,1,:,:))
;var1(:,2,24,44)=max(EVI(:,2,:,:))
;var1(:,3,24,44)=max(EVI(:,3,:,:))
;var1(:,4,24,44)=max(EVI(:,4,:,:))
;var1(:,5,24,44)=max(EVI(:,5,:,:))
;var1(:,6,24,44)=max(EVI(:,6,:,:))
;var1(:,7,24,44)=max(EVI(:,7,:,:))

;var2(:,:,:,:) = 0
;var2(:,0,24,44)=max(EVI_MAX(:,0,:,:)) ; nested domain
;var2(:,1,24,44)=max(EVI_MAX(:,1,:,:))
;var2(:,2,24,44)=max(EVI_MAX(:,2,:,:))
;var2(:,3,24,44)=max(EVI_MAX(:,3,:,:))
;var2(:,4,24,44)=max(EVI_MAX(:,4,:,:))
;var2(:,5,24,44)=max(EVI_MAX(:,5,:,:))
;var2(:,6,24,44)=max(EVI_MAX(:,6,:,:))
;var2(:,7,24,44)=max(EVI_MAX(:,7,:,:))

;var3(:,:,:,:) = 0
;var3(:,0,24,44)=max(EVI_MIN(:,0,:,:)) ; nested domain
;var3(:,1,24,44)=max(EVI_MIN(:,1,:,:))
;var3(:,2,24,44)=max(EVI_MIN(:,2,:,:))
;var3(:,3,24,44)=max(EVI_MIN(:,3,:,:))
;var3(:,4,24,44)=max(EVI_MIN(:,4,:,:))
;var3(:,5,24,44)=max(EVI_MIN(:,5,:,:))
;var3(:,6,24,44)=max(EVI_MIN(:,6,:,:))
;var3(:,7,24,44)=max(EVI_MIN(:,7,:,:))

;var4(:,:,:,:) = 0
;var4(:,0,24,44)=max(LSWI(:,0,:,:))
;var4(:,1,24,44)=max(LSWI(:,1,:,:)) ; nested domain
;var4(:,2,24,44)=max(LSWI(:,2,:,:))
;var4(:,3,24,44)=max(LSWI(:,3,:,:))
;var4(:,4,24,44)=max(LSWI(:,4,:,:))
;var4(:,5,24,44)=max(LSWI(:,5,:,:))
;var4(:,6,24,44)=max(LSWI(:,6,:,:))
;var4(:,7,24,44)=max(LSWI(:,7,:,:))

;var5(:,:,:,:) = 0
;var5(:,0,24,44)=max(LSWI_MAX(:,0,:,:))
;var5(:,1,24,44)=max(LSWI_MAX(:,1,:,:)) ; nested domain
;var5(:,2,24,44)=max(LSWI_MAX(:,2,:,:))
;var5(:,3,24,44)=max(LSWI_MAX(:,3,:,:))
;var5(:,4,24,44)=max(LSWI_MAX(:,4,:,:))
;var5(:,5,24,44)=max(LSWI_MAX(:,5,:,:))
;var5(:,6,24,44)=max(LSWI_MAX(:,6,:,:))
;var5(:,7,24,44)=max(LSWI_MAX(:,7,:,:))

;var6(:,:,:,:) = 0
;var6(:,0,24,44)=max(LSWI_MIN(:,0,:,:))
;var6(:,1,24,44)=max(LSWI_MIN(:,1,:,:)) ; nested domain
;var6(:,2,24,44)=max(LSWI_MIN(:,2,:,:))
;var6(:,3,24,44)=max(LSWI_MIN(:,3,:,:))
;var6(:,4,24,44)=max(LSWI_MIN(:,4,:,:))
;var6(:,5,24,44)=max(LSWI_MIN(:,5,:,:))
;var6(:,6,24,44)=max(LSWI_MIN(:,6,:,:))
;var6(:,7,24,44)=max(LSWI_MIN(:,7,:,:))

;var7(:,:,:,:) = 0
;var7(:,0,24,44)=max(VEGFRA_VPRM(:,0,:,:))
;var7(:,1,24,44)=max(VEGFRA_VPRM(:,1,:,:)) ; nested domain
;var7(:,2,24,44)=max(VEGFRA_VPRM(:,2,:,:))
;var7(:,3,24,44)=max(VEGFRA_VPRM(:,3,:,:))
;var7(:,4,24,44)=max(VEGFRA_VPRM(:,4,:,:))
;var7(:,5,24,44)=max(VEGFRA_VPRM(:,5,:,:))
;var7(:,6,24,44)=max(VEGFRA_VPRM(:,6,:,:))
;var7(:,7,24,44)=max(VEGFRA_VPRM(:,7,:,:))


; output variables directly
out->EVI= var1
out->EVI_MAX= var2
out->EVI_MIN= var3
out->LSWI = var4
out->LSWI_MAX = var5
out->LSWI_MIN = var6
out->VEGFRA_VPRM = var7


out->XLAT = XLAT
out->XLONG =XLONG
out->Times = Times

; assign named dimensions
EVI!0="Time"
EVI!1="vprm_vgcls"
EVI!2="south_north"
EVI!3="west_east"
EVI@FieldType = 104
EVI@MemoryOrder = "XYZ"
EVI@coordinates="XLAT XLONG"
EVI@description="EMISSIONS"
EVI@stagger=""
EVI@units = "mol km^-2 hr^-1"

EVI_MAX!0="Time"
EVI_MAX!1="vprm_vgcls"
EVI_MAX!2="south_north"
EVI_MAX!3="west_east"
EVI_MAX@FieldType = 104
EVI_MAX@MemoryOrder = "XYZ"
EVI_MAX@coordinates="XLAT XLONG"
EVI_MAX@description="EMISSIONS"
EVI_MAX@stagger=""
EVI_MAX@units = "mol km^-2 hr^-1"

EVI_MIN!0="Time"
EVI_MIN!1="vprm_vgcls"
EVI_MIN!2="south_north"
EVI_MIN!3="west_east"
EVI_MIN@FieldType = 104
EVI_MIN@MemoryOrder = "XYZ"
EVI_MIN@coordinates="XLAT XLONG"
EVI_MIN@description="EMISSIONS"
EVI_MIN@stagger=""
EVI_MIN@units = "mol km^-2 hr^-1"

LSWI!0="Time"
LSWI!1="vprm_vgcls"
LSWI!2="south_north"
LSWI!3="west_east"
LSWI@FieldType = 104
LSWI@MemoryOrder = "XYZ"
LSWI@coordinates="XLAT XLONG"
LSWI@description="EMISSIONS"
LSWI@stagger=""
LSWI@units = "mol km^-2 hr^-1"

LSWI_MAX!0="Time"
LSWI_MAX!1="vprm_vgcls"
LSWI_MAX!2="south_north"
LSWI_MAX!3="west_east"
LSWI_MAX@FieldType = 104
LSWI_MAX@MemoryOrder = "XYZ"
LSWI_MAX@coordinates="XLAT XLONG"
LSWI_MAX@description="EMISSIONS"
LSWI_MAX@stagger=""
LSWI_MAX@units = "mol km^-2 hr^-1"

LSWI_MIN!0="Time"
LSWI_MIN!1="vprm_vgcls"
LSWI_MIN!2="south_north"
LSWI_MIN!3="west_east"
LSWI_MIN@FieldType = 104
LSWI_MIN@MemoryOrder = "XYZ"
LSWI_MIN@coordinates="XLAT XLONG"
LSWI_MIN@description="EMISSIONS"
LSWI_MIN@stagger=""
LSWI_MIN@units = "mol km^-2 hr^-1"

VEGFRA_VPRM!0="Time"
VEGFRA_VPRM!1="vprm_vgcls"
VEGFRA_VPRM!2="south_north"
VEGFRA_VPRM!3="west_east"
VEGFRA_VPRM@FieldType = 104
VEGFRA_VPRM@MemoryOrder = "XYZ"
VEGFRA_VPRM@coordinates="XLAT XLONG"
VEGFRA_VPRM@description="EMISSIONS"
VEGFRA_VPRM@stagger=""
VEGFRA_VPRM@units = "mol km^-2 hr^-1"



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

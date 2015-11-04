; load functions and procedures

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

 
begin

 

; WRF ARW input file (NOTE, your wrfout file does not need
; the .nc, but NCL needs it so make sure to add it in the
; line below)

; a = addfile("./wrfout_d02_2006-feb1-5_1way_lowres_diurnoff.nc","r")
; a = addfile("./wrfout_edg41_diurnoff_1way_d01_2006-feb1-5.nc","r")
; a =  addfile("./wrfout_d01_edg3_2006-feb1-5_1way_lowres_diurn_diff.nc","r")
; a =  addfile("./wrfout_d02_edg3_2006-feb1-5_1way_lowres_diurn_diff.nc","r")
; a = addfile("./wrfout_d02_edg41_2006-feb1-5_1way_lowres_diurnon.nc","r")
; a = addfile("./wrfout_d01_edg41_2006-feb1-5_1way_lowres_diurnoff.nc","r")
; a = addfile("./wrfout_CO2_richmond_src_d01_2006-feb1-5.nc","r")
a = addfile("./wrfout_d01_2006-feb1-27_co2_ant_bio_run4.nc","r")
; Output on screen. Output will be called "plt_Surface1"

; type = "x11"
type = "pdf"
; type = "ps"
; type = "ncgm"

; wks = gsn_open_wks(type,"wrfout_CO2_ANT_Wind_feb1-27_2006_d01")
wks = gsn_open_wks(type,"wrfout_EBIO_RES_Wind_feb1-27_2006_d01")
; wks = gsn_open_wks(type,"wrfout_d02_2006-feb1-5_SST_test")
; wks = gsn_open_wks(type,"wrfout_d01_edg3_k2_2006-feb1-5_1way_lowres_diurn_diff") 
; wks = gsn_open_wks(type,"wrfout_d02_edg41_k1_2006-feb1-5_1way_lowres_diurnon")
; wks = gsn_open_wks(type,"wrfout_d01_edg41_k0_2006-feb1-5_1way_lowres_diurnoff")
; wks = gsn_open_wks(type,"wrfout_d02_edg3_k2_2006-feb1-5_1way_lowres_diurn_diff") 
; Set basic resources

res = True

res@MainTitle = "Feb 1-27,2006 EBIO_RES,u,v"          ; Give plot a main title

res@Footer = False                              ; Set Footers off
pltres = True                                ; Plotting resources
mpres = True                                      ; Map resources

 

;---------------------------------------------------------------

times = wrf_user_getvar(a,"times",-1)    ; get all times in the file
; times = wrf_user_getvar(a,"times",(/5,121,1/))
ntimes = dimsizes(times)                  ; number of times in the file
; it = (/5,121,1/)                            ; only interested in first time
; do it=8,ntimes-1
do it=1,40
;print(it)
;print(times)
;print(times(it))
print("Working on time: " + times(it) )
res@TimeLabel = times(it)            ; keep some time information

   

;---------------------------------------------------------------

; Get variables

 
; CO2_ANT = wrf_user_getvar(a,"CO2_ANT",-1)                ; Get CO2_ANT
; co2_ant = wrf_user_getvar(a,"CO2_ANT",-1)
; co2_bio = wrf_user_getvar(a,"CO2_BIO",-1)
ebio_res =  wrf_user_getvar(a,"EBIO_RES",-1)


u10 = wrf_user_getvar(a,"U10",it)                       ; Get U10

v10 = wrf_user_getvar(a,"V10",it)                       ; Get V10

   u10 = u10*1.94386                           ; Convert to knots

   v10 = v10*1.94386

   u10@units = "kts"

   v10@units = "kts"

 
;---------------------------------------------------------------

;  ; Plotting options for EBIO_RES
 opts=res
 opts@cnFillOn = True
 opts@cnLinesOn = True
 opts@lbLabelAutoStride    = True    ; Nice spacing for labelbar labels.
 opts@lbBoxLinesOn         = True   ; Turn off labelbar box lines
 opts@gsnSpreadColorEnd = -3  ; End third from the last color in color map
 opts@lbLabelFontHeightF = 0.02
 opts@InitTime = True                     ; Do not plot time or footers
 opts@ValidTime = True
 opts@Footer = True
 opts@gsnSpreadColors      = True               ; span full colormap

; opts@cnLevelSelectionMode = "ManualLevels"
; opts@gsnSpreadColors      = True               ; span full colormap
; opts@cnMinLevelValF = 500.0
; opts@cnMaxLevelValF = 25000.0
; opts@cnLevelSpacingF = 3


 mpres@mpGeophysicalLineColor = "Black"
 mpres@mpNationalLineColor    = "Black"
 mpres@mpUSStateLineColor     = "Black"
 mpres@mpGridLineColor = "Black"
 mpres@mpLimbLineColor = "Black"
 mpres@mpPerimLineColor = "Black"
 
 mpres@mpOutlineOn = True
 mpres@mpGeophysicalLineThicknessF = 2.0
 mpres@mpUSStateLineThicknessF = 2.0 

 ; contour_co2ant = wrf_contour(a,wks,CO2_ANT(it,1,:,:),opts)
 ; contour_co2ant = wrf_contour(a,wks,co2_ant(it,0,:,:),opts)
 contour_ebiores = wrf_contour(a,wks,ebio_res(it,0,:,:),opts)
 opts@gsnStringFontHeightF = 0.030
 res@tmXBLabelFontHeightF = 0.037                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.037
 delete(opts)

 opts = res
 ; contour_co2ant = wrf_contour(a,wks,CO2_ANT(it,1,:,:),opts)
 ; contour_co2bio = wrf_contour(a,wks,co2_bio(it,0,:,:),opts)
 contour_ebiores  = wrf_contour(a,wks,ebio_res(it,0,:,:),opts) 
 opts@gsnStringFontHeightF = 0.030
 res@tmXBLabelFontHeightF = 0.037                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.037
 delete(opts)


; Plotting options for Wind Vectors

opts = res                                  ; Add basic resources

opts@FieldTitle = "Winds"             ; Overwrite the field title

opts@NumVectors = 47                      ; Density of wind barbs

vector = wrf_vector(a,wks,u10,v10,opts)             ; Create plot

delete(opts)

 

; MAKE PLOTS 

; plot = wrf_map_overlays(a,wks,  \
      ; (/contour_co2ant,contour_tc,contour_psl,vector/),pltres,mpres)
      ; (/contour_co2ant,vector/),pltres,mpres)

 plot = wrf_map_overlays(a,wks, (/contour_ebiores,vector/),pltres,mpres)
;---------------------------------------------------------------

 
end do
 

end
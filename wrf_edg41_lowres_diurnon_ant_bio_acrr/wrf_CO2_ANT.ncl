; load functions and procedures

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

 

; WRF ARW input file (NOTE, your wrfout file does not need
; the .nc, but NCL needs it so make sure to add it in the
; line below)

a = addfile("./wrfout_CO2_richmond_src_d01_2006-feb1-5.nc","r")

; Output on screen. Output will be called "plt_Surface1"

; type = "x11"
type = "pdf"
; type = "ps"
; type = "ncgm"


wks = gsn_open_wks(type,"CO2_Richmond_ant_concent_feb1-5_2006_d01")
 

; Set basic resources

res = True

res@MainTitle = "Feb 1-5,2006 Richmod Test CO2 Conc."          ; Give plot a main title

res@Footer = False                              ; Set Footers off
pltres = True                                ; Plotting resources
mpres = True                                      ; Map resources

 

;---------------------------------------------------------------

times = wrf_user_getvar(a,"times",-1)    ; get all times in the file
ntimes = dimsizes(times)                  ; number of times in the file
; it = (/5,121,1/)                            ; only interested in first time
do it=1,ntimes-1
; do it=1,2
;print(it)
;print(times)
;print(times(it))
print("Working on time: " + times(it) )
res@TimeLabel = times(it)            ; keep some time information

   

;---------------------------------------------------------------

; Get variables

 
CO2_ANT = wrf_user_getvar(a,"CO2_ANT",-1)                ; Get CO2_ANT
; co2_ant = wrf_user_getvar(a,"E_CO2",-1)
; co2_ant = wrf_user_getvar(a,"E_CO2",(/2/))

;---------------------------------------------------------------

; Plotting options for E_CO2
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
 opts@cnLevelSelectionMode = "ManualLevels"
 opts@gsnSpreadColors      = True               ; span full colormap
 opts@cnMinLevelValF = 377.0
 opts@cnMaxLevelValF = 440.0
 opts@cnLevelSpacingF = 10
 mpres@mpGeophysicalLineColor = "Black"
 mpres@mpNationalLineColor    = "Black"
 mpres@mpUSStateLineColor     = "Black"
 mpres@mpGridLineColor = "Black"
 mpres@mpLimbLineColor = "Black"
 mpres@mpPerimLineColor = "Black"
 
 mpres@mpOutlineOn = True
 mpres@mpGeophysicalLineThicknessF = 2.0
 mpres@mpUSStateLineThicknessF = 2.0 

 contour_co2ant = wrf_contour(a,wks,CO2_ANT(it,0,:,:),opts)
 opts@gsnStringFontHeightF = 0.030
 res@tmXBLabelFontHeightF = 0.037                 ; resize tick labels
 res@tmYLLabelFontHeightF = 0.037
 delete(opts)
 

; MAKE PLOTS 

 plot = wrf_map_overlays(a,wks, (/contour_co2ant/),pltres,mpres)
;---------------------------------------------------------------

 
end do
 

end
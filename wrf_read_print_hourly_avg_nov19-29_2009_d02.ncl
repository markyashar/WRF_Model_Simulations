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

;open netcdf data file to get data from 
data = addfile("./wrfout_d02_jm_slc_ant_nov1-29_2009_ic_bc.nc","r")

CO2_ANT_00h = wrf_user_getvar(data,"CO2_ANT",(/456,480,504,528,552,576,600,624,648,672/))
CO2_ANT_01h = wrf_user_getvar(data,"CO2_ANT",(/457,481,505,529,553,577,601,625,649,673/))
CO2_ANT_02h = wrf_user_getvar(data,"CO2_ANT",(/458,482,506,530,554,578,602,626,650,674/))
CO2_ANT_03h = wrf_user_getvar(data,"CO2_ANT",(/459,483,507,531,555,579,603,627,651,675/))
CO2_ANT_04h = wrf_user_getvar(data,"CO2_ANT",(/460,484,508,532,556,580,604,628,652,676/))
CO2_ANT_05h = wrf_user_getvar(data,"CO2_ANT",(/461,485,509,533,557,581,605,629,653,677/))
CO2_ANT_06h = wrf_user_getvar(data,"CO2_ANT",(/462,486,510,534,558,582,606,630,654,678/))
CO2_ANT_07h = wrf_user_getvar(data,"CO2_ANT",(/463,487,511,535,559,583,607,631,655,679/))
CO2_ANT_08h = wrf_user_getvar(data,"CO2_ANT",(/464,488,512,536,560,584,608,632,656,680/))
CO2_ANT_09h = wrf_user_getvar(data,"CO2_ANT",(/465,489,513,537,561,585,609,633,657,681/))
CO2_ANT_10h = wrf_user_getvar(data,"CO2_ANT",(/466,490,514,538,562,586,610,634,658,682/))
CO2_ANT_11h = wrf_user_getvar(data,"CO2_ANT",(/467,491,515,539,563,587,611,635,659,683/))
CO2_ANT_12h = wrf_user_getvar(data,"CO2_ANT",(/468,492,516,540,564,588,612,636,660,684/))
CO2_ANT_13h = wrf_user_getvar(data,"CO2_ANT",(/469,493,517,541,565,589,613,637,661,685/))
CO2_ANT_14h = wrf_user_getvar(data,"CO2_ANT",(/470,494,518,542,566,590,614,638,662,686/))
CO2_ANT_15h = wrf_user_getvar(data,"CO2_ANT",(/471,495,519,543,567,591,615,639,663,687/))
CO2_ANT_16h = wrf_user_getvar(data,"CO2_ANT",(/472,496,520,544,568,592,616,640,664,688/))
CO2_ANT_17h = wrf_user_getvar(data,"CO2_ANT",(/473,497,521,545,569,593,617,641,665,689/))
CO2_ANT_18h = wrf_user_getvar(data,"CO2_ANT",(/474,498,522,546,570,594,618,642,666,690/))
CO2_ANT_19h = wrf_user_getvar(data,"CO2_ANT",(/475,499,523,547,571,595,619,643,667,691/))  
CO2_ANT_20h = wrf_user_getvar(data,"CO2_ANT",(/476,500,524,548,572,596,620,644,668,692/))
CO2_ANT_21h = wrf_user_getvar(data,"CO2_ANT",(/477,501,525,549,573,597,621,645,669,693/))
CO2_ANT_22h = wrf_user_getvar(data,"CO2_ANT",(/478,502,526,550,574,598,622,646,670,694/))
CO2_ANT_23h = wrf_user_getvar(data,"CO2_ANT",(/479,503,527,551,575,599,623,647,671,695/))
CO2_ANT_24h = wrf_user_getvar(data,"CO2_ANT",(/480,504,528,552,576,600,624,648,672,696/))


CO2_ANT_avg_hourly_nov1_14_00h = dim_avg_n_Wrap(CO2_ANT_00h,0)+380
CO2_ANT_avg_hourly_nov1_14_01h = dim_avg_n_Wrap(CO2_ANT_01h,0)+380
CO2_ANT_avg_hourly_nov1_14_02h = dim_avg_n_Wrap(CO2_ANT_02h,0)+380
CO2_ANT_avg_hourly_nov1_14_03h = dim_avg_n_Wrap(CO2_ANT_03h,0)+380
CO2_ANT_avg_hourly_nov1_14_04h = dim_avg_n_Wrap(CO2_ANT_04h,0)+380
CO2_ANT_avg_hourly_nov1_14_05h = dim_avg_n_Wrap(CO2_ANT_05h,0)+380
CO2_ANT_avg_hourly_nov1_14_06h = dim_avg_n_Wrap(CO2_ANT_06h,0)+380
CO2_ANT_avg_hourly_nov1_14_07h = dim_avg_n_Wrap(CO2_ANT_07h,0)+380
CO2_ANT_avg_hourly_nov1_14_08h = dim_avg_n_Wrap(CO2_ANT_08h,0)+380
CO2_ANT_avg_hourly_nov1_14_09h = dim_avg_n_Wrap(CO2_ANT_09h,0)+380
CO2_ANT_avg_hourly_nov1_14_10h = dim_avg_n_Wrap(CO2_ANT_10h,0)+380
CO2_ANT_avg_hourly_nov1_14_11h = dim_avg_n_Wrap(CO2_ANT_11h,0)+380
CO2_ANT_avg_hourly_nov1_14_12h = dim_avg_n_Wrap(CO2_ANT_12h,0)+380
CO2_ANT_avg_hourly_nov1_14_13h = dim_avg_n_Wrap(CO2_ANT_13h,0)+380
CO2_ANT_avg_hourly_nov1_14_14h = dim_avg_n_Wrap(CO2_ANT_14h,0)+380
CO2_ANT_avg_hourly_nov1_14_15h = dim_avg_n_Wrap(CO2_ANT_15h,0)+380
CO2_ANT_avg_hourly_nov1_14_16h = dim_avg_n_Wrap(CO2_ANT_16h,0)+380
CO2_ANT_avg_hourly_nov1_14_17h = dim_avg_n_Wrap(CO2_ANT_17h,0)+380
CO2_ANT_avg_hourly_nov1_14_18h = dim_avg_n_Wrap(CO2_ANT_18h,0)+380
CO2_ANT_avg_hourly_nov1_14_19h = dim_avg_n_Wrap(CO2_ANT_19h,0)+380
CO2_ANT_avg_hourly_nov1_14_20h = dim_avg_n_Wrap(CO2_ANT_20h,0)+380
CO2_ANT_avg_hourly_nov1_14_21h = dim_avg_n_Wrap(CO2_ANT_21h,0)+380
CO2_ANT_avg_hourly_nov1_14_22h = dim_avg_n_Wrap(CO2_ANT_22h,0)+380
CO2_ANT_avg_hourly_nov1_14_23h = dim_avg_n_Wrap(CO2_ANT_23h,0)+380
CO2_ANT_avg_hourly_nov1_14_24h = dim_avg_n_Wrap(CO2_ANT_24h,0)+380 

print ("00h" + "	"  + CO2_ANT_avg_hourly_nov1_14_00h(0,43,31))
print ("01h" + "	"  + CO2_ANT_avg_hourly_nov1_14_01h(0,43,31))
print ("02h" + "	"  + CO2_ANT_avg_hourly_nov1_14_02h(0,43,31))
print ("03h" + "	"  + CO2_ANT_avg_hourly_nov1_14_03h(0,43,31))
print ("04h" + "	"  + CO2_ANT_avg_hourly_nov1_14_04h(0,43,31))
print ("05h" + "	"  + CO2_ANT_avg_hourly_nov1_14_05h(0,43,31))
print ("06h" + "	"  + CO2_ANT_avg_hourly_nov1_14_06h(0,43,31))
print ("07h" + "	"  + CO2_ANT_avg_hourly_nov1_14_07h(0,43,31))
print ("08h" + "	"  + CO2_ANT_avg_hourly_nov1_14_08h(0,43,31))
print ("09h" + "	"  + CO2_ANT_avg_hourly_nov1_14_09h(0,43,31))
print ("10h" + "	"  + CO2_ANT_avg_hourly_nov1_14_10h(0,43,31))
print ("11h" + "	"  + CO2_ANT_avg_hourly_nov1_14_11h(0,43,31))
print ("12h" + "	"  + CO2_ANT_avg_hourly_nov1_14_12h(0,43,31))
print ("13h" + "	"  + CO2_ANT_avg_hourly_nov1_14_13h(0,43,31))
print ("14h" + "	"  + CO2_ANT_avg_hourly_nov1_14_14h(0,43,31))
print ("15h" + "	"  + CO2_ANT_avg_hourly_nov1_14_15h(0,43,31))
print ("16h" + "	"  + CO2_ANT_avg_hourly_nov1_14_16h(0,43,31))
print ("17h" + "	"  + CO2_ANT_avg_hourly_nov1_14_17h(0,43,31))
print ("18h" + "	"  + CO2_ANT_avg_hourly_nov1_14_18h(0,43,31))
print ("19h" + "	"  + CO2_ANT_avg_hourly_nov1_14_19h(0,43,31))
print ("20h" + "	"  + CO2_ANT_avg_hourly_nov1_14_20h(0,43,31))
print ("21h" + "	"  + CO2_ANT_avg_hourly_nov1_14_21h(0,43,31))
print ("22h" + "	"  + CO2_ANT_avg_hourly_nov1_14_22h(0,43,31))
print ("23h" + "	"  + CO2_ANT_avg_hourly_nov1_14_23h(0,43,31))
print ("24h" + "	"  + CO2_ANT_avg_hourly_nov1_14_24h(0,43,31))

end 

import numpy as np
import pylab as pl
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.dates import date2num
import datetime
data = np.loadtxt('wrf_read_print_hourly_avg_nov19-29_2009_d01.txt',usecols = (0,1))
print data[:,0]
print data[:,1]
# plot the first column as x, and second column as y
fig = pl.figure()
ax = fig.add_subplot(111)
for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontsize(8.5)
my_xticks = ['00h','01h','02h','03h','04h','05h','06h','07h','08h','09h','10h','11h','12h','13h','14h','15h','16h','17h','18h','19h','20h','21h','22h','23h','24h']
pl.xticks(data[:,0],my_xticks)
pl.plot(data[:,0], data[:,1], 'ro')
pl.xlabel('hour of day')
pl.ylabel('Total hourly mean Anthropogenic CO2 concentration (ppm)')
pl.title('Total hourly mean CO2 Concentration, SLC, November 19-29, 2009 (all times UTC), Parent Domain', fontsize = 8.0)
pl.grid(True)
plt.savefig('CO2_ANTHRO_TOTAL_nov19-29_2009_SLC_d01_hourly_mean.png', format = 'png',  bbox_inches='tight', pad_inches=0, dpi=300)

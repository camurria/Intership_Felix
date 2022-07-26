from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import pyscamp as MP
import mpl_toolkits.mplot3d

t0 = 1126259462.4
size = 4.
cut = 0.5

hdata = TimeSeries.fetch_open_data('H1', t0-size,t0+size)
ldata = TimeSeries.fetch_open_data('L1', t0-size,t0+size)
#fig0 = plt.figure(1)
#fig0 = plt.plot(hdata)
#plot.show() #questo per plottare la time series

#per la parte del filtro
from gwpy.signal import filter_design
notches = [filter_design.notch(line, hdata.sample_rate) for line in (60, 120, 180)]

bp_h = filter_design.bandpass(50, 250, hdata.sample_rate)
zpk_h = filter_design.concatenate_zpks(bp_h, *notches)
hfilt = hdata.filter(zpk_h, filtfilt=True)
hfilt = hfilt.crop(t0-size+cut, t0+size-cut)

bp_l = filter_design.bandpass(50,250,ldata.sample_rate)
zpk_l = filter_design.concatenate_zpks(bp_l, *notches)
lfilt = ldata.filter(zpk_l, filtfilt=True)
lfilt = lfilt.crop(t0-size+cut, t0+size-cut) #currently lenght = 1638 data points

Han = hfilt*1e20
Liv = lfilt*1e20
index = np.arange(len(Han))

plt.figure(0)
plt.subplot(1,1,1)
plt.plot(index,hfilt)
plt.title("Hanford data for GW150914")
plt.xlabel('Time (s)')
plt.ylabel('Strain (dimensionless)')

#plt.subplot(1,2,2)
#plt.plot(lfilt)
#plt.title("Livingstone data for GW150914")
#plt.xlabel('Time (s)')
#plt.ylabel('Strain (dimensionless)')
print(type(Han))
i=1
Sub=[]
output=[]
Tim =[]
M=1000.
matr = MP.selfjoin(Han,M)
matrix = matr[0]

#while i < 8:
#	M = 200*i
#	Sub.append(M)
#	print('M = {} s'.format(M*hdata.dt.value))
#	matrix_and_index = MP.selfjoin(Han,M)
#	matrix = matrix_and_index[0]
#	while len(matrix) < len(Han):
#		matrix = np.append(matrix,0)
#	output.append(matrix[25000])
#	i = i+1

plt.figure(1)
plt.plot(matrix)
#ax.title('Matrix profile Han Selfjoin M = {}'.format(M))
#plt.scatter(Sub,output)
plt.show()


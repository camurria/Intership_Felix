from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import pyscamp as MP

t0 = 1126259462.4
size = 3.5
cut = 0.5

sdata = TimeSeries.fetch_open_data('H1', t0-size,t0+size) #signal data
bdata = TimeSeries.fetch_open_data('H1', t0+size,t0+3*size) #background data
#fig0 = plt.figure(1)
#fig0 = plt.plot(hdata)
#plot.show() #questo per plottare la time series

#per la parte del filtro
from gwpy.signal import filter_design
notches = [filter_design.notch(line, sdata.sample_rate) for line in (60, 120, 180)]

bp_s = filter_design.bandpass(50, 250, sdata.sample_rate)
zpk_s = filter_design.concatenate_zpks(bp_s, *notches)
sfilt = sdata.filter(zpk_s, filtfilt=True)
sfilt = sfilt.crop(t0-size+cut, t0+size-cut)

bp_b = filter_design.bandpass(50,250,bdata.sample_rate)
zpk_b = filter_design.concatenate_zpks(bp_b, *notches)
bfilt = bdata.filter(zpk_b, filtfilt=True)
bfilt = bfilt.crop(t0+size+cut, t0+3*size-cut)
Sign = sfilt*1e20
Back = bfilt*1e20

index = np.arange(len(Sign))
print(index)

plt.figure(0)
plt.subplot(1,2,1)
plt.plot(index,sfilt)
plt.title("Signal data for GW150914")
plt.xlabel('Time (s)')
plt.ylabel('Strain (dimensionless)')
plt.subplot(1,2,2)
plt.plot(index,bfilt)
plt.title("Background data")
plt.xlabel('Time (s)')
plt.ylabel('Strain (dimensionless)')

i=1

for M in [1000]:
	print('Size of M: {}={} s'.format(M,M*sdata.dt.value))
	Join = MP.abjoin(Sign,Back,M)
	matrix_join = Join[0]
	Self = MP.selfjoin(Back,M)
	matrix_self = Self[0]
	times = t0-size+cut +(np.arange(len(Sign)-M+1))*Sign.dt.value
	contrast = matrix_join - matrix_self
	#print('Matrix join:',matrix_join)
	#print('Matrix self:',matrix_self)
	
	plt.figure(i)
	plt.subplot(1,2,1)
	plt.plot(matrix_join)
	plt.title('Joined MP w/ subl {}'.format(M))
	plt.subplot(1,2,2)
	plt.plot(matrix_self)
	plt.title('Self MP w/ subl {}'.format(M))
	
	plt.figure(i+1)
	plt.plot(contrast)
	plt.title('Contrast Profile w/ subl {}'.format(M))
	i = i+2
	
plt.show()


from gwpy.timeseries import TimeSeries
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import pyscamp as MP

t0 = 1126259462.4
hdata = TimeSeries.fetch_open_data('H1', t0-0.4,t0+0.4)
ldata = TimeSeries.fetch_open_data('L1', t0-0.4,t0+0.4)
#fig0 = plt.figure(1)
#fig0 = plt.plot(hdata)
#plot.show() #questo per plottare la time series

#per la parte del filtro

from gwpy.signal import filter_design
notches = [filter_design.notch(line, hdata.sample_rate) for line in (60, 120, 180)]

bp_h = filter_design.bandpass(50, 250, hdata.sample_rate)
zpk_h = filter_design.concatenate_zpks(bp_h, *notches)
hfilt = hdata.filter(zpk_h, filtfilt=True)
hfilt = hfilt.crop(t0-0.2,t0+0.2)

bp_l = filter_design.bandpass(50,250,ldata.sample_rate)
zpk_l = filter_design.concatenate_zpks(bp_l, *notches)
lfilt = ldata.filter(zpk_l, filtfilt=True)
lfilt = lfilt.crop(t0-0.2,t0+0.2) #currently lenght = 1638 data points

M = 600
Han = hfilt*1e20
Liv = lfilt*1e20
matrix_and_index = MP.abjoin(Han,Liv,M)
matrix = matrix_and_index[0]
times = t0-0.2+np.arange(len(Han)-M+1)*hdata.dt.value
#print(times)

plt.figure(2)
plt.subplot(1,2,1)
plt.plot(hfilt)
plt.title("Hanford data for GW150914")
plt.xlabel('Time (s)')
plt.ylabel('Strain (dimensionless)')

plt.subplot(1,2,2)
plt.plot(lfilt)
plt.title("Livingstone data for GW150914")
plt.xlabel('Time (s)')
plt.ylabel('Strain (dimensionless)')

plt.figure(3)
plt.plot(times,matrix)
plt.title('Matrix profile Han vs Liv')

plt.show()


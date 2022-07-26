import pyscamp as MP
import matplotlib.pyplot as plt
import gwpy
import pylab
import numpy as np
from gwpy.timeseries import TimeSeries
from gwosc.datasets import event_gps
from pycbc.waveform import get_td_waveform
#-------------END: import packages

#------------START: Create Time Series
m = 36 # Solar masses
sample_rate = 4096 #Hz
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                     mass1=m,
                     mass2=m,
                     delta_t=1./sample_rate,
                     f_lower=30)

hp =(hp*5e18).crop(0,0.1) #length is now 906 data points
gauss= np.random.normal(loc=0.0,scale=2.,size=len(hp)*20) #length is appr. 18 000 data points >> len(hp)

only_signal = np.random.normal(loc=0.0,scale=2.,size=len(hp)*20)
times = np.arange(len(gauss))/float(sample_rate)
wf_start = np.random.randint(0, len(gauss) - len(hp))
gauss[wf_start:wf_start+len(hp)] += hp.numpy()
print(only_signal[wf_start:wf_start+len(hp)])
print(gauss[wf_start:wf_start+len(hp)])
#-------------END: Create Time Series

#------------START: Matrix Profile
mat_and_ind = MP.selfjoin(gauss,int(len(hp)))
matrix = mat_and_ind[0]
print(matrix, len(matrix))
#--------------END: Matrix Profile

#------------START: Contrast Profile
ABjoin = MP.abjoin(gauss,only_signal,int(len(hp)))
contrast = ABjoin[0] - matrix
print(ABjoin[0])
#--------------END: Contrast Profile

#----------START: plot setup
fig1=plt.figure(1)
fig2=plt.figure(1)
fig1=plt.subplot(1,2,1)
fig2=plt.subplot(1,2,2)

fig1.plot(gauss)
fig2.plot(hp.sample_times, gauss[wf_start:wf_start+len(hp)])
fig2.plot(hp.sample_times, hp)

fig3 = plt.figure(2)
fig3 = plt.plot(matrix)
plt.title('Matrix Profile (SelfJoin)')

plt.figure(4)
plt.plot(ABjoin[0])
plt.title('Joined Matrix w/ and w/o signal')

plt.figure(3)
plt.plot(contrast)
plt.title('Contrast Profile')
#----------END: plot setup
plt.show()

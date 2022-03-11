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
gauss= np.random.normal(loc=0.0,scale=2.,size=len(hp)*60) #length is appr. 18 000 data points >> len(hp)

times = np.arange(len(gauss))/float(sample_rate)
wf_start = np.random.randint(0, len(gauss) - len(hp))
gauss[wf_start:wf_start+len(hp)] += hp.numpy()
#-------------END: Create Time Series

#------------START: Matrix Profile
mat_and_ind = MP.selfjoin(gauss,len(hp))
matrix = mat_and_ind[0]
print(matrix, len(matrix))
#--------------END: Matrix Profile

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
#----------END: plot setup
plt.show()

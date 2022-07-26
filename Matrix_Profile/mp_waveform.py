import pyscamp as MP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gwpy
import pylab
from gwpy.timeseries import TimeSeries
from gwosc.datasets import event_gps
from pycbc.waveform import get_td_waveform

gps = event_gps('GW150914')
start=gps-10
end=gps+10

data = TimeSeries.fetch_open_data('H1',start=start,end=end,verbose=True,cache=True)

A = data*1e20 #A is the signal 

#print(A.value,'1/dt: ',data.sample_rate.value,'dt: ',data.dt.value)

m = 36 # Solar masses
hp, hc = get_td_waveform(approximant="SEOBNRv4_opt",
                     mass1=m,
                     mass2=m,
                     delta_t=data.dt.value,
                     f_lower=20)
                     
B = hp*1e20 #Time Series B is the expected waveform = query
B = B.crop(0.,0.1)
#hp = hp / max(np.correlate(hp, hp, mode='full'))**0.5

subl=[] #different sublenghts m to Matrix Profile
for x in [int(len(B)/4)]:
	subl.append(x)

print('Length of A (signal) = {} & length of B (query) = {}'.format(len(A),len(B)))
#pylab.plot(B.sample_times,B)
#pylab.show()

#print(A.sample_times)

for j in range(len(subl)):
	join = MP.abjoin_sum(A,B,subl[j])
	#joinp =np.pad(join,(len(A)-len(join)-m))
	fig = plt.subplot(1,len(subl),j+1)
	print('sublength = {}, length matrix = {}'.format(subl[j],len(join)))
	#print(join)
	plt.plot(join)
	
	minima=0
	i=0
	while i < len(join):
		if join[i]==join.min():
			print('minimum of {} found at {}'.format(join.min(),i))
			minima=i
		i = i+1

	t_min = start + minima*data.dt.value
	idx_gps = (gps-start)*data.sample_rate.value
	#idx_merg = 
	print('true merger time: {} and index {}, and found minimum: {}'.format(gps,idx_gps,t_min))

join_matrix = MP.abjoin_matrix(A,B,subl[0])
print(join_matrix)
#plt.show()

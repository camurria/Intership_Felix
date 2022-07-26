import pyscamp as MP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gwpy
from gwpy.timeseries import TimeSeries
from gwosc.datasets import event_gps

subl=[]
for m in [100,200,500]:
	subl.append(m)

gps = event_gps('GW190412')
segment = (int(gps)-1,int(gps)+1)
print('segment of time: ',segment,'centered on gps time: ',gps)
#sample_rate = 4096

#test=TimeSeries(np.random.random(1000),sample_rate=100)
data = TimeSeries.fetch_open_data('L1',*segment,verbose=True,cache=True)
data = data*1e20
print(data.value)
data.plot()

for j in range(len(subl)):
	profile, index = MP.selfjoin(data.value,subl[j])
	fig = plt.subplot(1,len(subl),j+1)
	print(profile)
	fig.plot(profile)


#print('profile min',profile.min())
#i=1
#mini=[]
#while i < len(data)-subl[0]+1:
#	if profile[i]==profile.min():
#		print('minimums: ',i)
#		mini.append(i)
#	i = i+1

#print(mini)
#fig1 = plt.subplot(1,3,1)
#fig2 = plt.subplot(1,3,2)
#fig3 = plt.subplot(1,3,3)


#fig1.plot(data)
#fig2.plot(profile)
#fig3.plot(profile)
#plt.xlim(mini[0]-subl[0]/2,mini[0]+subl[0]/2)

plt.show()

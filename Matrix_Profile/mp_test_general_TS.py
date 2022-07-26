import pyscamp as MP
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

has = MP.gpu_supported()

m = 500 #sublength

periods = 2000
np.random.seed([3,1415])
tidx = pd.date_range('2016-07-01', periods=periods, freq='T')
data = np.random.randn(periods)
ts = pd.Series(data=data, index=tidx, name='TimeSeriesTest')

fig1 = plt.subplot(1,3,1)
fig2 = plt.subplot(1,3,2)
fig3 = plt.subplot(1,3,3)

profile, index = MP.selfjoin(ts,m)
fig1.plot(ts)
fig2.plot(profile)

print(profile.min())
i=1
while i < periods-m:
	if profile[i]==profile.min():
		minimum = i
		print(minimum)
	i = i+1

fig3.plot(profile)
plt.xlim(minimum-m/2,minimum+m/2)
plt.show()

import numpy as np
from scipy import stats
import pandas as pd
import math
import matplotlib.pyplot as plt

def movstd(a,window):
    left= window[0]
    right = window[1]
    result = []
    for i in range(len(a)):
        r=0
        if i >= left and i +right<len(a):
            r=np.std(a[i - left:i + right+1],ddof=1)
        elif i < left:
            r= np.std(a[:i + right+1],ddof=1)
        elif i + right >= len(a):
            r= np.std(a[i - left:],ddof=1)
        else:
            r=np.std(a[:],ddof=1)
        result.append(0 if math.isnan(r)  else r)
    return result


def findInT(query,target):
    m = len(query);
    n = len(target);
    Q = stats.zscore(query,ddof=1)#zNorm ???
    stdv = movstd(target,(0,m-1))
    Q= np.append(Q[::-1],np.zeros(n-m))
    dots =np.convolve(target,Q)
    dist =2 * (m - (dots[m-1:n])/ stdv[0:n - m +1])
    return np.sqrt(dist)


#----------Part GW
from gwpy.timeseries import TimeSeries
t0 = 1126259462.4
size = 1.5
cut = 0.5

sdata = TimeSeries.fetch_open_data('H1', t0-size,t0+size) #signal data
bdata = TimeSeries.fetch_open_data('H1', t0+size,t0+3*size) #background data
print(len(sdata),len(bdata))
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

qry = sfilt*1e20
tgt = bfilt*1e20
#-----------------------

#qry = np.array([1, 10, 5]);
#tgt = [4,8,6,-1,-2,-3,-1,3,4,5];
output=findInT(qry,tgt)
print(output)
plt.figure(1)
plt.plot(output)
plt.show()
#disp(Stdv) 2.0000    2.0000    1.0000   18.7705   21.0000    4.2426         0
#Q = Q(end:-1:1); % Reverse
#the query
#disp(Q);
#Q(m + 1:n) = 0; % pad
#zeros
#disp(Q);
#dots = conv(T, Q);
#disp(X(m:n))
#dist = 2 * (m - (dots(m:n))./ Stdv(1:n - m + 1));
#dist = sqrt(dist);


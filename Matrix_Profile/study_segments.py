import h5py
import gwpy
import numpy as np
from pycbc.frame import read_frame
import matplotlib.pyplot as plt
import sys
import pyscamp as MP

type = str(sys.argv[1])
number = int(sys.argv[2])

if type == 'signal':
	f = h5py.File('/home/felixb/gravitational_waves/Intership_Felix/create_segments/signal_50_55.hdf5','r')
elif type == 'noise':
	f = h5py.File('/home/felixb/gravitational_waves/Intership_Felix/create_segments/noise_tests.hdf5','r')
elif type == 'glitch':
	f = h5py.File('/home/felixb/gravitational_waves/Intership_Felix/create_segments/glitches_tests.hdf5','r')
else:
	print('wrong argument : please choose between <signal>, <noise> or <glitch>')
	sys.exit()

groups = f.keys()
if number <= 0 or number > len(groups):
	print('Wrong argument : number of plot asked (int) should be between <1> and <{}>'.format(len(groups)))
	sys.exit()

name_metadata=['gps_start_segment','gps_begin_of_ts','duration',
'mass1','mass2','spin1z','spin2z','inclination','coa_phase',
'distance','right_ascension','declination','polarization','SNR']

data = [0]*number
times = [0]*number

for key,i in zip(groups,range(number)):
	grp=f[key]
	timeseries = grp['time series'][:]
	metadata = grp['metadata']
	gpsstart = metadata[1]
	duration = metadata[2]
	SNR = metadata[13]
	print('SNR',SNR)
	fs = 2048.
	t = np.arange(gpsstart, gpsstart + duration, 1./fs)
	data[i] = timeseries
	times[i] = t
f.close()

datan = [0]*number
n = h5py.File('/home/felixb/gravitational_waves/Intership_Felix/create_segments/noise_30_32.hdf5','r')
groupsn = n.keys()
for keyn,i_n in zip(groupsn,range(number)):
	grpn = n[keyn]
	timeseriesn = grpn['time series'][:]
	datan[i_n] = timeseriesn
n.close()
sig_count = 0
sig_tab = []
M = 100 #Sublength for Matrix Profile
for j in range(number):
	matrix_and_index_self = MP.selfjoin(data[j],M)
	Self = matrix_and_index_self[0]
	#matrix_and_index_join = MP.abjoin(data[j],datan[np.random.randint(1,number)],M)
	matrix_and_index_join = MP.abjoin(data[j],datan[j],M)
	Join = matrix_and_index_join[0]
	Contrast = Join - Self
	#plt.figure(j+1)
	#plt.suptitle('Segment #{} of length {}. Type: {}'.format(j,len(data[j]),type))
	#plt.subplot(1,2,1)
#	plt.plot(times[j],data[j])
#	plt.title('Time Series')
	sig_tab = np.append(sig_tab,np.std(Contrast))
#	plt.subplot(1,2,2)
#	plt.plot(times[j][:len(Contrast)],Contrast)
#	plt.title('Self joined Matrix Profile')
	#if np.std(Contrast) < 0.50:
	plt.figure(j)
	#plt.subplot(1,2,1)
	plt.plot(times[j][:len(Contrast)],Contrast)
	plt.xlabel('GPS Time [s]')
	plt.ylabel('Amplitude of the Matrix Profile')
	plt.title('Matrix Profile: Joined profile of noise and signal samples')
	#plt.subplot(1,2,2)
	#plt.hist(Contrast/np.std(Contrast)-np.mean(Contrast),bins=30)
	#plt.title('Contrast MP cumul')
	#print(j,'\t std:',np.std(Contrast),'\t mean:',np.mean(Contrast))
	#sig_count +=1

plt.figure()
plt.hist(sig_tab,bins = 50)
print('bad behaviors :',sig_count)
plt.show()

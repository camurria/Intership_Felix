import h5py
import gwpy
import numpy as np
import sys
from pycbc.frame import read_frame
import matplotlib.pyplot as plt

type = str(sys.argv[1])
number = int(sys.argv[2])

if type == 'signal':
	f = h5py.File('/home/felixb/matrix_profile/Intership_Felix/create_segments/signal_40_42.hdf5','r')
elif type == 'noise':
	f = h5py.File('/home/felixb/matrix_profile/Intership_Felix/create_segments/noise_tests.hdf5','r')
elif type == 'glitch':
	f = h5py.File('/home/felixb/matrix_profile/Intership_Felix/create_segments/glitches_tests.hdf5','r')
else: 
	print('wrong argument : please choose between <signal>, <noise> or <glitch>')
	sys.exit()
	
name_metadata=['gps_start_segment','gps_begin_of_ts','duration',
'mass1','mass2','spin1z','spin2z','inclination','coa_phase',
'distance','right_ascension','declination','polarization','SNR']
groups = f.keys()
i=0

if number <= 0 or number > len(groups):
	print('Wrong argument : number of plot asked (int) should be between <1> and <{}>'.format(len(groups)))
	sys.exit()


for key,i in zip(groups,range(0,number)):
	grp=f[key]
	data = grp['time series'][:]
	metadata = grp['metadata']
	gpsstart = metadata[1]
	duration = metadata[2]
	for j in range(len(metadata)):
		print('[{}] {} = {}'.format(j,name_metadata[j],metadata[j]))
	fs = 2048.
	t = np.arange(len(data))
	#t = np.arange(gpsstart, gpsstart + duration, 1./fs)
	print(len(t))
	print(metadata)
	plt.figure(i+1)
	plt.plot(t,data,label=key)
	print(key)
	i=i+1
f.close()
plt.show()

# METADATA
#1:start
#2:duration
#13:snr

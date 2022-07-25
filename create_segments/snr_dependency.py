import h5py
import gwpy
import numpy as np
import matplotlib.pyplot as plt
import sys
import pyscamp as MP
import os

sigfile = str(sys.argv[1])	#name of signal file
noifile = str(sys.argv[2])	#name of noise file
snr_min = int(sys.argv[3])	#value of minimum SNR
snr_max = int(sys.argv[4])	#value of maximum SNR
number = int(sys.argv[5])	#number of events computed
M = int(sys.argv[6]) 			#Sublength for Matrix Profile
data = [0]*number 		#Init empty table for signal timeseries
times = [0]*number		#Init empty table for times corr to timeseries
datan = [0]*number		#Init empty table for noise timeseries
name_metadata=['gps_start_segment','gps_begin_of_ts','duration',
'mass1','mass2','spin1z','spin2z','inclination','coa_phase',
'distance','right_ascension','declination','polarization','SNR']	#names of parameters stored in "metadata"
sig_count = 0			#A counter for values of sigma above certain treshold
thresh = 0.40			#Said threshold
sig_tab = []			#Init empty table for sigma values of Contrast Matrix Profile
spathfile = os.path.join('/home/felixb/matrix_profile/Intership_Felix/create_segments',sigfile)
npathfile = os.path.join('/home/felixb/matrix_profile/Intership_Felix/create_segments',noifile)

#________________Noise File: retrieve Time Series (stored in datan),
# to compute "Joined Matrix Profile", with signal segment + noisy segment
n = h5py.File(npathfile,'r')
groupsn = n.keys()
for keyn,i_n in zip(groupsn,range(number)):
	grpn = n[keyn]
	timeseriesn = grpn['time series'][:]
	datan[i_n] = timeseriesn
n.close()

#_______________Signal File
f = h5py.File(spathfile,'r')
groups = f.keys()

for key,i in zip(groups,range(number)):
	grp=f[key]
	timeseries = grp['time series'][:]
	metadata = grp['metadata']
	gpsstart = metadata[1]
	duration = metadata[2]
	fs = 2048.
	t = np.arange(gpsstart, gpsstart + duration, 1./fs)
	data[i] = timeseries #get time series for signal
	times[i] = t
	
	matrix_and_index_self = MP.selfjoin(data[i],M)	# Compute self joined MP
	#matrix_and_index_join = MP.abjoin(data[i],datan[np.random.randint(1,number)],M) #in case the JMP is done with random noise segment
	matrix_and_index_join = MP.abjoin(data[i],datan[i],M)	#Compute Joined MP
	Contrast = matrix_and_index_self[0] - matrix_and_index_join[0]				#Compute JMP
	
	sig_tab = np.append(sig_tab,np.std(Contrast))	#Store sigmas in a table
	print('step #{} of {}'.format(i+1,number))
	if sig_tab[i] < thresh:
		sig_count += 1
	
	plt.figure(i)
	plt.subplot(2,2,3)
	plt.hist(Contrast,bins=40)
	plt.subplot(2,2,2)
	plt.plot(timeseries)
	plt.subplot(2,2,1)
	plt.plot(matrix_and_index_join[0])
	plt.subplot(2,2,4)
	plt.plot(Contrast)
	print(snr_min,snr_max, sig_count)
	print(np.std(Contrast))
f.close()
#print(sig_tab)
print(snr_min,snr_max,sig_count)
#with open('out_m{}.txt'.format(M), 'a') as o:
#    print(snr_min,snr_max,sig_count, file=o)
#o.close()
    
plt.show()



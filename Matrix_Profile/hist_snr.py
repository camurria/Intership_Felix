import matplotlib.pyplot as plt
import numpy as np
import os,sys

qry = str(sys.argv[1])

if qry == 'std':
	var = ['50']
	#var = ['35','40','45','50']
elif qry == 'm':
	var = ['50','100','200','300']
else:
	print('arg = <std> or <m>')
	sys.exit()

fig = plt.figure()
markers = ['x','o','D','v','*','s','<','>','p','D']
colors = ['red','blue','green','orange','black','purple','brown','yellow','pink','grey']
i=0 #just to change scatter style between datasets
for t in var:
	snr_min=[]
	snr_max=[]
	count=[]
	with open('/home/felixb/gravitational_waves/Intership_Felix/create_segments/out_'+qry+t+'.txt','r') as f:
		lines = f.readlines()
		for x in lines:
			snr_min = np.append(snr_min,float(x.split(' ')[0]))
			snr_max = np.append(snr_max,float(x.split(' ')[1]))
			count = np.append(count,int(x.split(' ')[2]))
	f.close()
	width = snr_max[0]-snr_min[0]
	snr_cent = snr_min+width/2
	label = '{}: '.format(qry)+str(int(t))#/100)
	plt.scatter(snr_min,1-count/1000,s=20,marker=markers[i],color=colors[i],label=label)
	i += 1
	plt.legend(loc='upper left')
	plt.xlabel('SNR')
	plt.ylabel('metric: std of cumulated MP above threhsold of 0.40')
plt.show()

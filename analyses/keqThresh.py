## Check half functional rule for binding parameters
## Written by Madeline Galbraith
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import matplotlib.font_manager as font_manager
fontpath = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf'

prop = font_manager.FontProperties(fname=fontpath)
plt.rcParams['font.family'] = prop.get_name()



#############
## Check keq > 1/Sox2 and keq < 1/Sox2
## 9


files =["../sf1/stem_parameters_sf1.dat","../sf2/stem_parameters_sf2.dat","../sf3/stem_parameters_sf3.dat"]#"../sf4/stem_parameters_sf4.dat"]
dirs = ['sf1','sf2','sf3']
#files =["/home/madeline/Research/RACIPE/test_og/stem_parameters_og.dat","../sf1/stem_parameters_sf1.dat","../sf2/stem_parameters_sf2.dat","../sf3/stem_parameters_sf3.dat"]#,"../sf4/stem_parameters_sf4.dat"]
#dirs = ['og','sf1','sf2','sf3']

def getD(file,dir):

	par = pd.read_csv(filepath_or_buffer=file,header=None,sep='\t')
	sim = pd.read_csv(filepath_or_buffer='simple_'+dir+'.dat',header=0,sep=',')
	## sim has model, Gcnf, cdx2, gata6, pbx1, klf4, nanog, oct4, oct4sox2, sox2

	par = par.values
	sim = sim.values

	##  get the raw results for 
	## 	kp== 
	## 	si==complex/(sox2*oct4)
	## 	oc==oct4
	## 	os==complex
	## 	sx==sox2

	[kp,ind] = getP(par,'P')
	#[si,ind2] = getP(sim,'S')
	[oc,ind3] = getP(sim,'1')
	[os,ind4] = getP(sim,'2')
	[sx,ind5] = getP(sim,'3')
	## now we only want to keep the data that can be used for all sets
	## so take the union of all the sets to get the indices for the proper set
	ind = list(set(ind) & set(ind3) & set(ind4) & set(ind5)) # & set(ind2))
	kp = kp[np.sort(ind)]
	oc = 2**oc[np.sort(ind)]
	os = 2**os[np.sort(ind)]
	sx = 2**sx[np.sort(ind)]

	return [oc,os,sx,kp]

def getP(dat,type=''):
	length = len(dat[:,0])
	kp = np.zeros(10000)
	inds = []
	for i in range(length):
		index = (dat[i,0]-1)
		if type=='P' and index>=0 and index<10000:
			index = int(index)
			kp[index]=(dat[i,10]/dat[i,11])
			inds +=[index]
		elif type=='S' and index>=0 and index<10000:
			index = int(index)
			kp[index]=(dat[i,8]/(dat[i,7]*dat[i,9]))
			inds +=[index]

		elif type=='1' and index>=0 and index<10000:
			index = int(index)
			kp[index]=dat[i,7]
			inds +=[index]
		elif type=='2' and index>=0 and index<10000:
			index = int(index)
			kp[index]=dat[i,8]
			inds +=[index]
		elif type=='3' and index>=0 and index<10000:
			index = int(index)
			kp[index]=dat[i,9]
			inds +=[index]
	return [kp,inds]


def plotEverything(data,labels):
	plt.rc('font', size=20)     # fontsize of the axes title
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=50)    # fontsize of the x and y labels
	#plt.rc('legend', fontsize=20)    # legend fontsize
	plt.rc('figure', titlesize=20) 

	fig = plt.figure(figsize=(9,7))
	count=0
	for i in data:
		dd = data[i][1]/(data[i][0]*data[i][2])
		kp = data[i][3]
		if count>0:
			plt.bar(labels[i],dd*0.+1.,color='orange')#,label="r$_{eq}$>[Oct4Sox2]/[Oct4][Sox2]")
			plt.bar(labels[i],np.mean(kp<dd),color='g')#,label="r$_{eq}$<[Oct4Sox2]/[Oct4][Sox2]")
		else:
			plt.bar(labels[i],dd*0.+1.,color='orange',label="r$_{eq}$>[Oct4Sox2]/[Oct4][Sox2]")
			plt.bar(labels[i],np.mean(kp<dd),color='g',label="r$_{eq}$<[Oct4Sox2]/[Oct4][Sox2]")
		count+=1
	plt.axhline(y=0.5, color='r',linestyle='--')
	plt.ylim([0,1.5])
	plt.legend()
	plt.xticks(rotation=30)
	plt.ylabel("Probability",fontsize=18)
	plt.title('Percentage of r$_{eq}$ < or > [Oct4Sox2]/[Oct4][Sox2]',fontsize=18)
        plt.tight_layout()
	#plt.show()
	fig.savefig('percent_keq_sim.jpg')
	plt.close()

	'''
	fig = plt.figure(figsize=(9,7))
	count=0
	for i in data:
		dd = data[i][2]
		kp = data[i][3]
		if count>0:
			plt.bar(labels[i],dd*0.+1.,color='orange')
			plt.bar(labels[i],np.mean(kp<1./dd),color='g')
		else:
			plt.bar(labels[i],dd*0.+1.,color='orange',label="keq>1/[Sox2]")
			plt.bar(labels[i],np.mean(kp<1./dd),color='g',label="keq<1/[Sox2]")
		count+=1
	plt.axhline(y=0.5, color='r',linestyle='--')
	plt.ylim([0,1.5])
        plt.legend()
	plt.xticks(rotation=30)
	plt.title('Percentage of keq < or > 1/[Sox2]',fontsize=20)
	plt.tight_layout()
	#plt.show()
	fig.savefig('percent_keq_sox2.jpg')
	plt.close()
	
	fig = plt.figure(figsize=(9,7))
	count=0
	for i in data:
		dd = data[i][0]
		kp = data[i][3]
		if count>0:
		        plt.bar(labels[i],dd*0.+1.,color='orange')
			plt.bar(labels[i],np.mean(kp<1./dd),color='g')
		else:
		        plt.bar(labels[i],dd*0.+1.,color='orange',label="keq>1/[Oct4]")
			plt.bar(labels[i],np.mean(kp<1./dd),color='g',label="keq<1/[Oct4]")
		count+=1
	plt.axhline(y=0.5, color='r',linestyle='--')
	plt.ylim([0,1.5])
	plt.legend()
	plt.xticks(rotation=30)
	plt.title('Percentage of keq < or > 1/[Oct4]',fontsize=20)
	#plt.show()
	plt.tight_layout()
        fig.savefig('percent_keq_oct4.jpg')
	plt.close()
	'''

def plotE2(data,labels):
	plt.rc('font', size=20)     # fontsize of the axes title
	plt.rc('axes', titlesize=20)     # fontsize of the axes title
	plt.rc('axes', labelsize=50)    # fontsize of the x and y labels
	#plt.rc('legend', fontsize=20)    # legend fontsize
	plt.rc('figure', titlesize=20) 

	fig = plt.figure(figsize=(9,7))
	for i in data:
		dd = data[i][1]/(data[i][0]*data[i][2])
		kp = data[i][3]
		xx = np.abs(kp-dd)/(kp+dd)*200.
		plt.hist(xx,label=labels[i],alpha=0.4,bins=50)#,density=True)#,label="keq<[Oct4Sox2]/[Oct4][Sox2]")
	plt.legend()
	plt.xlim(-2,30)
	plt.ylabel("Frequency",fontsize=20)
	plt.xlabel("Percentage",fontsize=20)
	plt.title('Percent Difference between r$_{eq}$  and [Oct4Sox2]/[Oct4][Sox2]')
        plt.tight_layout()
	fig.savefig('hist_keq_sim.jpg')
	plt.close()

	'''
	fig = plt.figure(figsize=(9,7))
	for i in data:
		dd = data[i][2]
		kp = data[i][3]
		xx = np.abs(kp-dd)/(kp+dd)*200.
		plt.hist(xx,label=labels[i],alpha=0.4)#,label="keq<[Oct4Sox2]/[Oct4][Sox2]")
        plt.legend()
	plt.ylabel("Frequency",fontsize=20)
	plt.xlabel("Percentage",fontsize=20)
	plt.title('Percent Difference between keq  and 1/[Sox2]')
	plt.tight_layout()
	fig.savefig('hist_keq_sox2.jpg')
	plt.close()
	
	fig = plt.figure(figsize=(9,7))
	for i in data:
		dd = data[i][0]
		kp = data[i][3]
		xx = np.abs(kp-dd)/(kp+dd)*200.
		plt.hist(xx,label=labels[i],alpha=0.4)#,label="keq<[Oct4Sox2]/[Oct4][Sox2]")
	plt.legend()
	plt.xticks(rotation=30)
	plt.ylabel("Frequency",fontsize=20)
	plt.xlabel("Percentage",fontsize=20)
	plt.title('Percent Difference between keq  and 1/[Oct4]')
	plt.tight_layout()
        fig.savefig('hist_keq_oct4.jpg')
	plt.close()
	'''

	
files =["/home/madeline/Research/RACIPE/test_og/stem_parameters_og.dat","../sf1/stem_parameters_sf1.dat","../sf2/stem_parameters_sf2.dat","../sf3/stem_parameters_sf3.dat"]#,"../sf4/stem_parameters_sf4.dat"]
dirs = ['og','sf1','sf2','sf3']#,'sf4']
#files =["../sf1/stem_parameters_sf1.dat","../sf2/stem_parameters_sf2.dat","../sf3/stem_parameters_sf3.dat","../sf4/stem_parameters_sf4.dat"]
#dirs = ['sf1','sf2','sf3','sf4']


data ={}
for i in range(len(files)):
	data[i] = getD(files[i],dirs[i])
	print files[i]
print data.keys()

labs = ["RACIPE","RACIPE-wb (1)","RACIPE-wb (2)","RACIPE-wb (3)"]#,"RACIPE-wb (4)"]
plotEverything(data,labs)
plotE2(data,labs)

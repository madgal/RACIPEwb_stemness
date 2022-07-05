## Plot the entropy
## Written by Madeline Galbraith
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,copy

import matplotlib.font_manager as font_manager


fontpath = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf'

prop = font_manager.FontProperties(fname=fontpath)
plt.rcParams['font.family'] = prop.get_name()

params = {'legend.handlelength': 2}
plt.rcParams.update(params)



###

def plotEnt(labels,E,avgS,stitle):
	## get scaling for x value
	scaling = []
	for el in labels:
		if 'central' not in el:
			#print el.split('_')
			if 'over' in el:
				scaling +=[int(el.split('_')[-1][1:])*10.]
			else:
				scaling +=[1./(int(el.split('_')[-1][1:])*10.)]
				
		else:
			scaling +=[1.0]

	### reformat to be more useful

	set1 = {'Oct4':[[],[]],'Sox2':[[],[]],'Cdx2':[[],[]],'Gata6':[[],[]],'Nanog':[[],[]],'Pbx1':[[],[]],'Klf4':[[],[]],'Gcnf':[[],[]]}
	set2 = {'Oct4':[[],[]],'Sox2':[[],[]],'Cdx2':[[],[]],'Gata6':[[],[]],'Nanog':[[],[]],'Pbx1':[[],[]],'Klf4':[[],[]],'Gcnf':[[],[]]}
	set3 = {'Oct4':[[],[]],'Sox2':[[],[]],'Cdx2':[[],[]],'Gata6':[[],[]],'Nanog':[[],[]],'Pbx1':[[],[]],'Klf4':[[],[]],'Gcnf':[[],[]]}


	avgCentral=0
	for key in set1:
		for i in range(len(labels)):
			if key.lower() == labels[i].split('_')[0]:
				set1[key][0] += [scaling[i]]
				set1[key][1] += [E[i]]
				if scaling[i]<=1:
					set2[key][0] += [avgS[i]]
					set2[key][1] += [E[i]]
				else:
					set3[key][0] += [avgS[i]]
					set3[key][1] += [E[i]]
	
			else:
        			print key,labels[i]
	
        	for j in range(len(labels)):
			if 'central' in labels[j]:
				break
			
        	set1[key][0] += [scaling[j]]
        	set1[key][1] += [E[j]]
		set2[key][0] += [avgS[j]]
		set2[key][1] += [E[j]]
		avgCentral=avgS[j]

	set1 = sorting(set1)
	set2 = sorting(set2)
	set3 = sorting(set3)

	set4 =combine(set2,set3)
	fig, (ax1,ax2) = plt.subplots(1,2,figsize=(10,4))
	for k in set1:
		ax1.plot(set1[k][0],set1[k][1],'--.',label=k)
		ax2.plot(set4[k][0],set4[k][1],'--.',label=k)
        
	
        ax1.axvline(x=1,color='k',linestyle='--')
	ax2.axvline(x=avgCentral,color='k',linestyle='--')
	
	ax1.set_xscale('log')
        ax1.set_xlabel('Scaling of gene expression',fontsize=15)
	ax1.set_ylabel('Information entropy',fontsize=15)
        ax2.legend(loc='upper left',frameon=False)
	
	ax2.set_xlabel('Average number of stable states',fontsize=15)
        ax2.set_ylabel('Information entropy',fontsize=15)
	
	ax1.tick_params(labelsize=15)
	ax2.tick_params(labelsize=15)
        plt.show()
	fig.savefig('fig_entropy_final.png',bbox_inches='tight',dpi=1200)

def sorting(dict):

	for key in dict:
		a = dict[key][0]
		b = dict[key][1]

		a = np.array(a)
		b = np.array(b)

		ind = np.argsort(a)

		dict[key][0] = list(a[ind])
		dict[key][1] = list(b[ind])

	return dict

def combine(d1,d2):

	final={}
	for k in d1:
		final[k] = [[],[]]
		final[k][0] += d1[k][0]
		final[k][1] += d1[k][1]
		d2[k][0].reverse()
		d2[k][1].reverse()
		final[k][0] += d2[k][0]
		final[k][1] += d2[k][1]

	return final

#print pattern,',', -np.sum(probs*np.log(probs)),',',avgst

df = pd.read_csv('probs_6-18.dat',sep=',',header=0)

labels = df['Pattern'].values
avgS = df['avgS'].values

E=df['he'].values
plotEnt(labels,E,avgS,'e')

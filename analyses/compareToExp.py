import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt


def getAvg(filename):
	dfr = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t',skiprows=[1])
	data = dfr.values[:,4:]

	return np.mean(data,axis=0)

def getOG_wtc():

	avg={}
	for i in range(1,10):
		data = pd.read_csv("/home/madeline/Research/RACIPE/everythingFromDongya_EMT/c_code/euler/CIE/log/wt/Guo_exp_data/cluster_"+str(i)+".dat",header=None,sep='\t').values

		avg[i] = np.mean(data,axis=0)

        return avg

def getN_wtc():

	avg={}
	for i in range(1,14):
		data = pd.read_csv("Wt_sf3_cluster_"+str(i)+".dat",header=None,sep='\t').values
		tmp = np.array([data[:,1],data[:,2],data[:,4],data[:,5],data[:,6],data[:,8]])
		tmp = np.transpose(tmp)
		avg[i] = np.mean(tmp,axis=0)

        return avg


def generateRandom(filename):

	## get data
	#data = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t',skiprows=[1])
	data = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t')
	alldata={}
	#alldata['Gcnf'] = data['Gcnf'].values
	alldata['Cdx2'] = data['Cdx2'].values
	alldata['Gata6'] = data['Gata6'].values
	#alldata['Pbx1'] = data['Pbx1'].values
	alldata['Klf4'] = data['Klf4'].values
	alldata['Nanog'] = data['Nanog'].values
	alldata['Oct4'] = data['Oct4'].values
	#alldata['Oct4-Sox2'] = data['Oct4-Sox2'].values
	alldata['Sox2'] = data['Sox2'].values


	histD = {}
	## sample from histogram 
	## based on Daniel's answer at https://stackoverflow.com/questions/17821458/random-number-from-histogram
	for k in alldata:
		values,indices=np.histogram(alldata[k],bins=20)
		values=values.astype(np.float32)
		weights=values/np.sum(values)

		#Below, 5 is the dimension of the returned array.
		histD[k]=np.random.choice(indices[1:],10000,p=weights)

	randomV = np.empty((10000,6))
	randomV[:,0] = histD['Cdx2']
	randomV[:,1] = histD['Gata6']
	randomV[:,2] = histD['Klf4']
	randomV[:,3] = histD['Nanog']
	randomV[:,4] = histD['Oct4']
	randomV[:,5] = histD['Sox2']

	return randomV

def getMSE_XZ(avgClust,allDat):
	
	mse_xz={}
	ng = 9 ## genes

	for k in avgClust:
		mse_xz[k] =[]
		for i in range(len(allDat[:,0])):
			mse_xz[k] += [(1./ng)*np.mean((avgClust[k] - allDat[i])**2)]
	
	return mse_xz

def getSimilarity(avgC_og,avgC_wb):

	comp ={}
	ng=9

	for kog in avgC_og:
		comp[kog]={}
		for kwb in avgC_wb:
			
			#print kog,kwb,avgC_og[kog],avgC_wb[kwb]
			comp[kog][kwb] = (1./ng)*np.mean((avgC_og[kog] - avgC_wb[kwb])**2)

	return comp

def getData(filename):

	avgClust={}
	for file in os.listdir("."):
		if 'Wt' in file:
			clust = file.split('_')[2]
			clust = int(clust.split('.')[0])
			data = pd.read_csv(file,header=None,sep='\t').values
			avgClust[clust] = np.mean(data,axis=0)

	return avgClust
	
def compClust(mse_xz,mse_xy):

	simV={}
	for i in mse_xy:
		simV[i]=[]
		for j in mse_xy[i]:

			mse_xz[i] = np.array(mse_xz[i])
			temp= np.sum(mse_xy[i][j]<mse_xz[i])*1./len(mse_xz[i])
			#print temp,i,j
			if temp >=0.99:
				simV[i] +=[j]

	return simV



fout = open("ClusterCompResults.dat","w")
file_og='clusterSC_og.dat'
randomVector = generateRandom(file_og)
avgClust_og = getOG_wtc()
mse_xz=getMSE_XZ(avgClust_og,randomVector)

fout.write("Original:WithBinding\n")
avgClust_wb = getN_wtc()
mse_xy = getSimilarity(avgClust_og,avgClust_wb)

print mse_xy
simV = compClust(mse_xz,mse_xy)
for k in simV:
	print k, simV[k],len(simV[k])
		
for i in simV:
	fout.write(str(i)+": ")
	maxN = len(simV[i])-1
	if maxN>=0:
		for j in range(maxN):
			fout.write('%s,' %simV[i][j])
	fout.write('%s\n' %simV[i][maxN])
fout.close()

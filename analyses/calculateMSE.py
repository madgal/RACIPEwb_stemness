import pandas as pd
import numpy as np
import os
import sys


def getAvg(filename):
	dfr = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t',skiprows=[1])
	data = dfr.values[:,4:]

	return np.mean(data,axis=0)
def generateRandom(filename):

	## get data
	data = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t',skiprows=[1])
	alldata={}
	alldata['Gcnf'] = data['Gcnf'].values
	alldata['Cdx2'] = data['Cdx2'].values
	alldata['Gata6'] = data['Gata6'].values
	alldata['Pbx1'] = data['Pbx1'].values
	alldata['Klf4'] = data['Klf4'].values
	alldata['Nanog'] = data['Nanog'].values
	alldata['Oct4'] = data['Oct4'].values
	alldata['Oct4-Sox2'] = data['Oct4-Sox2'].values
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

	randomV = np.empty((10000,9))
	randomV[:,0] = histD['Gcnf']
	randomV[:,1] = histD['Cdx2']
	randomV[:,2] = histD['Gata6']
	randomV[:,3] = histD['Pbx1']
	randomV[:,4] = histD['Klf4']
	randomV[:,5] = histD['Nanog']
	randomV[:,6] = histD['Oct4']
	randomV[:,7] = histD['Oct4-Sox2']
	randomV[:,8] = histD['Sox2']

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
			
			comp[kog][kwb] = (1./ng)*np.mean((avgC_og[kog] - avgC_wb[kwb])**2)

	return comp

def getData(filename):


	fst = filename.split(".")[0]

	avgClust={}
	for file in os.listdir("."):
		if (fst in file) and ("_C" in file) and ("cdt" in file):
			clust = file.replace(fst,"")
			clust = clust.replace("_C","")
			clust = clust.replace(".cdt","")
			avgClust[clust] = getAvg(file)

	return avgClust
	
def compClust(mse_xz,mse_xy):

	simV={}
	for i in mse_xy:
		simV[i]=[]
		for j in mse_xy[i]:

			mse_xz[i] = np.array(mse_xz[i])
			temp= np.sum(mse_xy[i][j]<mse_xz[i])*1./len(mse_xz[i])
			if temp >=0.99:
				simV[i] +=[j]

	return simV



fout = open("ClusterCompResults.dat","w")
for file in os.listdir("."):
	file_og = "resultSC_og.cdt"
	randomVector = generateRandom(file_og)
	avgClust_og = getData(file_og)	
	mse_xz=getMSE_XZ(avgClust_og,randomVector)

	if ("result" in file) and ("cdt" in file) and ("_C" not in file):
		print file
		fout.write(file+"\n")
		fout.write("Original:WithBinding\n")
		avgClust_wb = getData(file)	
		mse_xy = getSimilarity(avgClust_og,avgClust_wb)

		simV = compClust(mse_xz,mse_xy)

		for i in simV:
			fout.write(str(i)+": ")
			maxN = len(simV[i])-1
			if maxN>=0:
				for j in range(maxN):
					fout.write(simV[i][j])
					fout.write(",")
				fout.write(simV[i][maxN])
			fout.write("\n")
	fout.flush()
fout.close()

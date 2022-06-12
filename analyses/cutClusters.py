import pandas as pd
import numpy as np
import os
import sys
import time
import matplotlib.pyplot as plt


def getData(filename):

        df_cdt = pd.read_csv(filepath_or_buffer=filename,header=0,sep='\t',skiprows=[1])
        filename = filename.replace('cdt','gtr')
	df_gtr = pd.read_csv(filepath_or_buffer=filename,header=None,sep='\t')

	if df_gtr.loc[0][0]=='NODEID':
		df_gtr = df_gtr.iloc[1:]
	if len(df_gtr.loc[1])==4:
		df_gtr.columns = ["NODEID","LEFT","RIGHT","CORRELATION"]
	else:
		df_gtr.columns = ["NODEID","LEFT","RIGHT","CORRELATION","NODECOLOR"]
	return [df_gtr,df_cdt]

def getEx(g1,g2,df_cdt):

	res=[]
	for gene in g1:
		if 'GENE' in gene:
			ind = np.argwhere(df_cdt.values[:,0]==gene)[0][0]	
			res +=[ list(df_cdt.values[ind,4:])]
	for gene in g2:
		if 'GENE' in gene:
			ind = np.argwhere(df_cdt.values[:,0]==gene)[0][0]	
			res +=[ list(df_cdt.values[ind,4:])]
	return res
def getG(df,color,df_cdt):
	
	data=df.values
	ind = np.argwhere(data[:,4]==color)[:,0]
	res= getEx(data[ind,1],data[ind,2],df_cdt)
	return res

def printFile(cluN,res):

	res = np.array(res)
	res[:,0] = (res[:,0])#-np.mean(res[:,0]))/np.std(res[:,0])
	res[:,1] = (res[:,1])#-np.mean(res[:,1]))/np.std(res[:,1])
	res[:,2] = (res[:,2])#-np.mean(res[:,2]))/np.std(res[:,2])
	res[:,3] = (res[:,3])#-np.mean(res[:,3]))/np.std(res[:,3])
	res[:,4] = (res[:,4])#-np.mean(res[:,4]))/np.std(res[:,4])
	res[:,5] = (res[:,5])#-np.mean(res[:,5]))/np.std(res[:,5])
	res[:,6] = (res[:,6])#-np.mean(res[:,6]))/np.std(res[:,6])
	res[:,7] = (res[:,7])#-np.mean(res[:,7]))/np.std(res[:,7])
	res[:,8] = (res[:,8])#-np.mean(res[:,8]))/np.std(res[:,8])


	fileout = open('Wt_sf3_cluster_'+str(cluN)+'.dat','w')
	for line in res:
		fileout.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[7],line[8]))

	fileout.close()

for file in os.listdir("."):
	if ("resultSC_sf3" in file) and ("cdt" in file) and ("_C" not in file):
		print file
		[df,df_cdt]=getData(file)	

labs=np.unique(df.values[:,4])
ctypes=[]
for el in labs:
	if type(el)=='float':
		print 'b'
	else:
		ctypes+=[el]
	
num = len(ctypes)
cnt,colors=[],[]
for i in range(0,num):

	id = np.argwhere(df.values[:,4]==ctypes[i])[:,0]
	tmpcnt=0
	for j in range(0,len(id)):
		for h in range(1,3):
			if 'GENE' in df.values[id[j],h]:
				tmpcnt+=1
	cnt+=[tmpcnt]
	colors+=[ctypes[i]]

res = np.array(cnt)*1./np.sum(cnt)
print np.sum(cnt)

wtC = {}
count=1
print np.sort(res)
for i in range(len(res)):
	if res[i]>=0.005:
		#wtC=getG(df,colors[i],df_cdt)
		print colors[i],count
		#printFile(count,wtC)
		count+=1
		#print count

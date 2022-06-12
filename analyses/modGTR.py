import pandas as pd
import numpy as np
import os
import sys
import time


def  getClust(dfg,dfc,comps):

	r1 = dfg.values[:,1]
	r2 = dfg.values[:,2]
	cols=dfg.values[:,4]

	clust = map2clust(r1,r2,cols,comps)

	#print comps
	#for i in range(len(r1)):
	#	print clust[i],cols[i],comps[cols[i]]
	return clust

def map2clust(gene1,gene2,colors,comps):

	clusts = []
	print 'Mapping'
	for i in range(len(gene1)):
		nt = colors[i]
		if nt =='#000000':
			clusts+=[nt]
		elif comps[nt]>=0.003:
			clusts+=[nt]
		else:
			clusts+=['#ffffff']
		
	return clusts


def modGTRFiles(file,clust,df):

	colors=["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#A4E804", "#004D43", "#8FB0FF", "#3B5DFF", "#B903AA", "#BA0900", "#FFB500", "#6B7900", "#FF90C9", "#BC23FF", "#99ADC0"]

	if 'NODECOLOR' not in df.columns:
		df['NODECOLOR'] = np.empty(len(df['NODEID']),dtype='object')

	#for i in clust:
	for i in range(len(clust)):
		df['NODECOLOR'][i] = clust[i]
		#df['NODECOLOR'][i] = colors[(clust[i])%len(colors)]
		#print clust[i],colors[(clust[i])%len(colors)]

	df.to_csv(file,sep='\t',index=False)

def getComp(df):
	labs = np.unique(df.values[:,4])
	ctypes = []
	for el in labs:
		if type(el)=='float':
			print 'b'
		else:
			ctypes+=[el]
	num=len(ctypes)
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

	cnt = np.array(cnt)
	res={}
	for i in range(len(cnt)):
		res[colors[i]] = cnt[i]*1./np.sum(cnt)

	return res
		
		
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

	df_gtr=df_gtr.reset_index()
	df_gtr=df_gtr.drop(columns=['index'])
	comps = getComp(df_gtr)
	clust = getClust(df_gtr,df_cdt,comps)

	modGTRFiles(filename,clust,df_gtr)


for file in os.listdir("."):
	if ("result" in file) and ("cdt" in file) and ('SC_sf3' in file):
		print file
		getData(file)	

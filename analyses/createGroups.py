import pandas as pd
import numpy as np
import os
import sys


def add(res,n,l,r):

	id=-1
	if not res:
		res[1]=[n,l,r]
		id=1

	else:
		skN,skL,skR=False,False,False
		for k in res.keys():
			if (n in res[k]):
				if (id!=-1) and (id!=k):
					print 'Eror in indexing???'
					exit()
				else:
					id=k
					skN=True
			if (l in res[k]):
				if (id!=-1) and (id!=k):
					print 'Eror in indexing???'
					exit()
				else:
					id=k
					skL=True
			if (r in res[k]):
				if (id!=-1) and (id!=k):
					print 'Eror in indexing???'
					exit()
				else:
					id=k
					skR=True

			if skN & skL & skR:
				break

		if id==-1:
			id = np.max(res.keys())+1
			res[id]=[]

		if not skN:
			res[id]+=[n]
		if not skL:
			res[id]+=[l]
		if not skR:
			res[id]+=[r]


	return [res,id]
					


def getClust(data):

	results ={}
	data = data.values
	
	## select out the data
	node = data[:,0]
	left = data[:,1]
	right = data[:,2]
	corr = data[:,3]
	corr = corr.astype(float)
	cl = {}

	for j in range(len(node)):
		i = len(node)-j-1
		
		if corr[i] > 0.8:
			[results,clust] = add(results,node[i],left[i],right[i])
			cl[i] = clust
	return [results,cl]

def getRes(data,df):

	labs = df.values[:,0]
	indice = []
	for el in data:
		if "GENE" in el:
			ind = np.argwhere(el==labs)[0]
			indice+=[ind[0]]


	return indice
	
def writeCDTFiles(filename,res,df_cdt):

	ff = filename.replace('gtr','cdt')

	## save the header
	df = pd.read_csv(filepath_or_buffer=ff,sep='\t')
	df = df.iloc[:1]

	for k in res.keys():
		end = "_C"+str(k)+".cdt"
		filex = filename.replace('.gtr',end)

		fi=open(filex,'w') 
		fi.write(df.to_csv(index=False,sep='\t'))
		dd = getRes(res[k],df_cdt)
		for ind in dd:
			sf = pd.DataFrame(df_cdt, index=[ind])
			fi.write(sf.to_csv(header=False,index=False,sep='\t'))

		fi.close()
def modGTRFiles(file,clust,df):

	colors=["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#A4E804", "#004D43", "#8FB0FF", "#3B5DFF", "#B903AA", "#BA0900", "#FFB500", "#6B7900", "#FF90C9", "#BC23FF", "#99ADC0"]

	if 'NODECOLOR' not in df.columns:
		df['NODECOLOR'] = np.empty(len(df['NODEID']),dtype='object')

	for i in clust:
		df['NODECOLOR'][i] = colors[(clust[i])%len(colors)]
		#print clust[i],colors[(clust[i])%len(colors)]

	df.to_csv(file,sep='\t',index=False)

		
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
	[res,clust] = getClust(df_gtr)

	print np.max(res.keys())
	writeCDTFiles(filename,res,df_cdt)
	modGTRFiles(filename,clust,df_gtr)

for file in os.listdir("."):
	## generate the groups
	if ("result" in file) and ("cdt" in file) and ("_C" not in file):
		print file
		getData(file)	

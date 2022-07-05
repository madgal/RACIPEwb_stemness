## Find and label clusters from hierarchical clustering
## Written by Madeline Galbraith
import pandas as pd
import numpy as np
import os
import sys

def main(filename,df,groups):
        fp = open(filename,'w')
        for index,row in df.iterrows():
                id = row['NAME']
                fp.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %(row['GWEIGHT'],id,row['Gcnf'],row['Cdx2'],row['Gata6'],row['Pbx1'],row['Klf4'],row['Nanog'],row['Oct4'],row['Oct4-Sox2'],row['Sox2']))
        fp.close()

def getData(filelist,countlist):
	print countlist,filelist
        count=int(countlist[0])
        df = pd.read_csv(filepath_or_buffer=filelist[0],header=0,sep='\t',skiprows=[1])
	df = df.astype({'GWEIGHT': 'int32'})
        df['GWEIGHT']=df['GWEIGHT']*0.+count
        for i in range(1,len(filelist)):
                count=int(countlist[i])
                dia = pd.read_csv(filepath_or_buffer=filelist[i],header=0,sep='\t',skiprows=[1])
		dia = dia.astype({'GWEIGHT': 'int32'})
                dia['GWEIGHT']=dia['GWEIGHT']*0.+count
                df = pd.concat([df,dia],axis=0,ignore_index=True)


        df=df.replace([np.inf,-np.inf],np.nan)
        df=df.dropna()
        return df
def getModelInfo(filename):
        data = pd.read_csv(filepath_or_buffer=filename,header=None,sep='\t').values
        g={}
        for i in range(len(data[:,1])):
                g[data[i,0]] = data[i,1]
        return g



groups = getModelInfo("/home/madeline/Research/RACIPE/test_og/stem_parameters_og.dat")
filelistA,filelistB=[],[]
countA,countB=[],[]
for file in os.listdir("."):
	if "result" in file and "og" in file and "cdt" in file and "_C" in file and "SC" in file:
		filelistA+=[file]
		countA+=[file.replace(".cdt","").split("_")[2][1:]]
	elif "result" in file and "og" in file and "cdt" in file and "_C" in file:
		filelistB+=[file]
		countB+=[file.replace(".cdt","").split("_")[2][1:]]

outputA = "clustPCA_SC_og.dat"
print filelistA,countA
df = getData(filelistA,countA)
main(outputA,df,groups)

'''
print 'Makde it past A'
outputB = "clustPCA_og.dat"
print filelistB,countB
df = getData(filelistB,countB)
main(outputB,df,groups)
'''


print 'Made it past OG'

numF = 3
for i in range(1,numF+1):

	groups = getModelInfo("../sf"+str(i)+"/stem_parameters_sf"+str(i)+".dat")
	
	countA,countB=[],[]
	filelistA,filelistB=[],[]
	for file in os.listdir("."):
		if ("result" in file) and ("sf"+str(i)+"_" in file) and ("cdt" in file) and ("_C" in file) and ("N" in file):
			filelistA+=[file]
			countA+=[file.replace(".cdt","").split("_")[2][1:]]
		elif ("result" in file) and ("sf"+str(i)+"_" in file) and ("cdt" in file) and ("_C" in file) and ("SC" in file):
			filelistB+=[file]
			countB+=[file.replace(".cdt","").split("_")[2][1:]]

	outputA = "clustPCA_N_sf"+str(i)+".dat"
	print i, filelistA,countA
	df = getData(filelistA,countA)
	main(outputA,df,groups)

	outputB = "clustPCA_SC_sf"+str(i)+".dat"
	print i, filelistB,countB
	df = getData(filelistB,countB)
	main(outputB,df,groups)

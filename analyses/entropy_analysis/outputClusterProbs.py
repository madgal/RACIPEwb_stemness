import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,sys
from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()


###

meanC,stdC=0,0

def main(data,avgst,pattern,clustersWT,clustersWTO):

	title = 'cl_'+pattern+'.dat'
	fileo = open(title,'w')
	## first determine cluster
	fileo.write('Model,prob,1/ss,randomB,clusterP,clusterSS')
	for i in range(len(data)):
		ci = map2Clust(clustersWT,data[i][4:])
		ciO = map2Clust(clustersWTO,data[i][4:])

		#print ci,ciO
		#print data[i],'\n'
		fileo.write("%s,%s,%s,%s,%s,%s\n' %(data[i][0],data[i][1],data[i][2],data[i][3],ci,ciO))
		#probs[ci]+=data[i][1]
		#probsO[ci]+=data[i][2]


	fileo.close()


def map2Clust(clustersWT,datapt):

	tmpD=[]
	for i in range(len(clustersWT)):
			tmpD += [np.min(np.sum((datapt-clustersWT[i])**2,axis=1))]

	ind = np.argwhere(tmpD==np.min(tmpD))[0][0]
	return  (ind+1)

def getClusterRes_OLD():
	c1 = pd.read_csv("../wt_prev/Wt_cluster_1.dat",header=None,sep='\t').values
	c2 = pd.read_csv("../wt_prev/Wt_cluster_2.dat",header=None,sep='\t').values
	c3 = pd.read_csv("../wt_prev/Wt_cluster_3.dat",header=None,sep='\t').values
	c4 = pd.read_csv("../wt_prev/Wt_cluster_4.dat",header=None,sep='\t').values
	c5 = pd.read_csv("../wt_prev/Wt_cluster_5.dat",header=None,sep='\t').values
	c6 = pd.read_csv("../wt_prev/Wt_cluster_6.dat",header=None,sep='\t').values
	c7 = pd.read_csv("../wt_prev/Wt_cluster_7.dat",header=None,sep='\t').values
	c8 = pd.read_csv("../wt_prev/Wt_cluster_8.dat",header=None,sep='\t').values
	c9 = pd.read_csv("../wt_prev/Wt_cluster_9.dat",header=None,sep='\t').values
	c10 = pd.read_csv("../wt_prev/Wt_cluster_10.dat",header=None,sep='\t').values
	c11 = pd.read_csv("../wt_prev/Wt_cluster_11.dat",header=None,sep='\t').values
	c12 = pd.read_csv("../wt_prev/Wt_cluster_12.dat",header=None,sep='\t').values
	c13 = pd.read_csv("../wt_prev/Wt_cluster_13.dat",header=None,sep='\t').values
	c14 = pd.read_csv("../wt_prev/Wt_cluster_14.dat",header=None,sep='\t').values
	c15 = pd.read_csv("../wt_prev/Wt_cluster_15.dat",header=None,sep='\t').values

	return [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15]


def getClusterRes():
	c1 = pd.read_csv("../Wt_cluster_1.dat",header=None,sep='\t').values
	c2 = pd.read_csv("../Wt_cluster_2.dat",header=None,sep='\t').values
	c3 = pd.read_csv("../Wt_cluster_3.dat",header=None,sep='\t').values
	c4 = pd.read_csv("../Wt_cluster_4.dat",header=None,sep='\t').values
	c5 = pd.read_csv("../t_cluster_5.dat",header=None,sep='\t').values
	c6 = pd.read_csv("../Wt_cluster_6.dat",header=None,sep='\t').values
	c7 = pd.read_csv("../Wt_cluster_7.dat",header=None,sep='\t').values
	c8 = pd.read_csv("../Wt_cluster_8.dat",header=None,sep='\t').values
	c9 = pd.read_csv("../Wt_cluster_9.dat",header=None,sep='\t').values
	c10 = pd.read_csv("../Wt_cluster_10.dat",header=None,sep='\t').values
	c11 = pd.read_csv("../Wt_cluster_11.dat",header=None,sep='\t').values
	c12 = pd.read_csv("../Wt_cluster_12.dat",header=None,sep='\t').values
	c13 = pd.read_csv("../Wt_cluster_13.dat",header=None,sep='\t').values
	c14 = pd.read_csv("../Wt_cluster_14.dat",header=None,sep='\t').values
	c15 = pd.read_csv("../Wt_cluster_15.dat",header=None,sep='\t').values
	c16 = pd.read_csv("../Wt_cluster_16.dat",header=None,sep='\t').values

	return [c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16]

def getGroups(file):
	parms = pd.read_csv(file,sep='\t',header=None)
	models = parms.values[:,0]
	weights = parms.values[:,-10:]

	paramW = {}
	paramO = {}
	for i in range(len(models)):
		paramW[models[i]] = weights[i]/100
		paramO[models[i]] = 1./parms.values[i,1]

	avgStates = np.mean((parms.values[:,1]<=6)*parms.values[:,1])
	return [paramW,paramO,avgStates]

def normalize(df,mean,std):
	d2 = df.values

	d2[:,3]  =  (d2[:,3] - meanC[3] )/stdC[3]
	d2[:,4]  =  (d2[:,4] - meanC[4] )/stdC[4]
	d2[:,5]  =  (d2[:,5] - meanC[5] )/stdC[5]
	d2[:,6]  =  (d2[:,6] - meanC[6] )/stdC[6]
	d2[:,7]  =  (d2[:,7] - meanC[7] )/stdC[7]
	d2[:,8]  =  (d2[:,8] - meanC[8] )/stdC[8]
	d2[:,9]  =  (d2[:,9] - meanC[9] )/stdC[9]
	d2[:,10] = (d2[:,10] - meanC[10])/stdC[10]
	d2[:,11] = (d2[:,11] - meanC[11])/stdC[11]


	rows= len(d2)
	arr = np.zeros(shape=(rows,10))

	arr[:,0] = d2[:,2]
	arr[:,1] = d2[:,4]
	arr[:,2] = d2[:,5]
	arr[:,3] = d2[:,3]
	arr[:,4] = d2[:,8]
	arr[:,5] = d2[:,6]
	arr[:,6] = d2[:,7]
	arr[:,7] = d2[:,9]
	arr[:,8] = d2[:,10]
	arr[:,9] = d2[:,11]


	return arr


def normalize_og(df):
	d2 = df.values

	meanC=[np.mean(d2[:,0]),0,np.mean(d2[:,2]),np.mean(d2[:,3]),np.mean(d2[:,4]),np.mean(d2[:,5]),np.mean(d2[:,6]),np.mean(d2[:,7]),np.mean(d2[:,8]),np.mean(d2[:,9]),np.mean(d2[:,10]),np.mean(d2[:,11])]
	stdC =[np.std(d2[:,0]),1, np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9]), np.std(d2[:,10]), np.std(d2[:,11])]
	d2[:,3]  =  (d2[:,3] - meanC[3] )/stdC[3]
	d2[:,4]  =  (d2[:,4] - meanC[4] )/stdC[4]
	d2[:,5]  =  (d2[:,5] - meanC[5] )/stdC[5]
	d2[:,6]  =  (d2[:,6] - meanC[6] )/stdC[6]
	d2[:,7]  =  (d2[:,7] - meanC[7] )/stdC[7]
	d2[:,8]  =  (d2[:,8] - meanC[8] )/stdC[8]
	d2[:,9]  =  (d2[:,9] - meanC[9] )/stdC[9]
	d2[:,10] = (d2[:,10] - meanC[10])/stdC[10]
	d2[:,11] = (d2[:,11] - meanC[11])/stdC[11]


	rows= len(d2)
	arr = np.zeros(shape=(rows,10))

	arr[:,0] = d2[:,2]
	arr[:,1] = d2[:,4]
	arr[:,2] = d2[:,5]
	arr[:,3] = d2[:,3]
	arr[:,4] = d2[:,8]
	arr[:,5] = d2[:,6]
	arr[:,6] = d2[:,7]
	arr[:,7] = d2[:,9]
	arr[:,8] = d2[:,10]
	arr[:,9] = d2[:,11]


	return [arr,meanC,stdC]

def readData(filestart):
	df1 = pd.read_csv(filepath_or_buffer=filestart+'_1.dat',header=None,sep='\t')#,names=['slen','swid','plen','pwid','target'])
	temp1= df1[0]*0.+1.
	df = pd.concat([temp1,temp1,df1[0],df1[1],df1[2],df1[3],df1[4],df1[5],df1[6],df1[7],df1[8],df1[9]],axis=1,ignore_index=True)

	###  2
	di = pd	.read_csv(filepath_or_buffer=filestart+'_2.dat',header=None,sep='\t')#,names=['slen','swid','plen','pwid','target'])
	temp2= di[0]*0.+2.
	dia = pd.concat([temp2,temp2*0.+1,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	dia = pd.concat([temp2,temp2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	
	## solutions with 3 steady states
	di = pd.read_csv(filepath_or_buffer=filestart+'_3.dat',header=None,sep='\t')
	temp = di[0]*0.+3.
	dia = pd.concat([temp,temp*0+1.,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	dia = pd.concat([temp,temp*0+2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	dia = pd.concat([temp,temp*0+3.,di[0],di[19],di[20],di[21],di[22],di[23],di[24],di[25],di[26],di[27]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	
	## solutions with 4 steady states
	di = pd.read_csv(filepath_or_buffer=filestart+'_4.dat',header=None,sep='\t')
	temp = di[0]*0.+4.
	dia = pd.concat([temp,temp*0+1,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
	dia = pd.concat([temp,temp*0+2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+3.,di[0],di[19],di[20],di[21],di[22],di[23],di[24],di[25],di[26],di[27]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
	dia = pd.concat([temp,temp*0+4.,di[0],di[28],di[29],di[30],di[31],di[32],di[33],di[34],di[35],di[36]],axis=1,ignore_index=True)
    	df = pd.concat([df,dia],axis=0,ignore_index=True)

	## solutions with 5 steady states
	di = pd.read_csv(filepath_or_buffer=filestart+'_5.dat',header=None,sep='\t')
	temp = di[0]*0.+5.
	dia = pd.concat([temp,temp*0+1,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
 	dia = pd.concat([temp,temp*0+2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+3.,di[0],di[19],di[20],di[21],di[22],di[23],di[24],di[25],di[26],di[27]],axis=1,ignore_index=True)
 	df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+4.,di[0],di[28],di[29],di[30],di[31],di[32],di[33],di[34],di[35],di[36]],axis=1,ignore_index=True)
       	df = pd.concat([df,dia],axis=0,ignore_index=True)
 	dia = pd.concat([temp,temp*0+5.,di[0],di[37],di[38],di[39],di[40],di[41],di[42],di[43],di[44],di[45]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
	
	## solutions with 6 steady states
	di = pd.read_csv(filepath_or_buffer=filestart+'_6.dat',header=None,sep='\t')
	temp = di[0]*0.+6.
	dia = pd.concat([temp,temp*0+1,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+3.,di[0],di[19],di[20],di[21],di[22],di[23],di[24],di[25],di[26],di[27]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+4.,di[0],di[28],di[29],di[30],di[31],di[32],di[33],di[34],di[35],di[36]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+5.,di[0],di[37],di[38],di[39],di[40],di[41],di[42],di[43],di[44],di[45]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)
        dia = pd.concat([temp,temp*0+6.,di[0],di[46],di[47],di[48],di[49],di[50],di[51],di[52],di[53],di[54]],axis=1,ignore_index=True)
        df = pd.concat([df,dia],axis=0,ignore_index=True)

	return df
	
def getData(filestart,groups,og_groups,selected,mean,std):

	data = readData(filestart)
	data = normalize(data,mean,std)
	allG = []
	allO = []

	final = []
	## paramW = groups

	modelN = data[:,0]

	for k in selected:

		ind = np.argwhere(modelN==k)[:,0]
		
		if len(ind)>0:
			for i in range(len(ind)):
				inc = ind[i]
			
				if i==selected[k]:
					final +=[np.append([data[inc,0],groups[k][selected[k]],og_groups[k],1],data[inc,1:] )]
				else:
					final +=[np.append([data[inc,0],groups[k][i],og_groups[k],0],data[inc,1:] )]

	return np.array(final)

	
def getDataCen(filestart,groups,og_groups,selected):

	data = readData(filestart)
	[data,mean,std] = normalize_og(data)
	allG = []
	allO = []

	final = []
	## paramW = groups

	modelN = data[:,0]

	for k in selected:

		ind = np.argwhere(modelN==k)[:,0]
		
		if len(ind)>0:
			for i in range(len(ind)):
				inc = ind[i]
			
				if i==selected[k]:
					final +=[np.append([data[inc,0],groups[k][selected[k]],og_groups[k],1],data[inc,1:] )]
				else:
					final +=[np.append([data[inc,0],groups[k][i],og_groups[k],0],data[inc,1:] )]

	return [np.array(final),mean,std]

def getRandomSample(prob,og_p):

	selected={}
	for k in prob:
		for i in range(int(1./og_p[k])):
			ss = int(prob[k][i]*100)
			if i==0:
                                l = np.ones(shape=ss)*i
                        elif ss>0:
                                l = np.insert(l,0,np.ones(shape=ss)*i)


		try:
                	temp = l.copy()
 	                for i in range(9):
        	                l = np.insert(l,0,temp)
                	np.random.shuffle(l)

 	                i = np.random.randint(0,len(l))
			selected[k] = int(l[i])
		except:
			selected[k]=0

        return  selected

##################
clusters    = getClusterRes()
clustersOld = getClusterRes_OLD()

res=[]
	
[groups,og_groups,avgState] = getGroups('../stem_parameters_0.dat')
selected = getRandomSample(groups,og_groups)
[data,meanC,stdC] = getDataCen('../stem_solution_0',groups,og_groups,selected)

if my_rank==0:
	main(data,avgState,'central',clusters,clustersOld)

exit()
	
dirs = ['gata6','gcnf','cdx2','klf4','nanog','pbx1','oct4','sox2']
scales = ['s1','s2','s3','s4','s5','s6','s7','s8','s9','s10']

myListB={0:['down'],1:['down'],2:['over'],3:['over']}
myListC={0:['s1','s2','s3','s4','s5'],1:['s6','s7','s8','s9','s10'],2:['s1','s2','s3','s4','s5'],3:['s6','s7','s8','s9','s10']}

for a in dirs:
        for b in myListB[my_rank]:#['down','over']:
                for c in myListC[my_rank]:#scales:
			pattern = a+'_'+b+'_'+c
			fst = '../'+a+'/'+b+'/stem_solution_'+c
			[groups,og_groups,avgState] = getGroups('../'+a+'/'+b+'/stem_parameters_'+c+'.dat')
			selected = getRandomSample(groups,og_groups)
			data = getData(fst,groups,og_groups,selected,meanC,stdC)
			main(data,avgState,pattern,clusters,clustersOld)

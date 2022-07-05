## Parse and collate data files for analysis
## Written by Madeline Galbraith
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os,sys

def main(file0,file1,file2,file3,file4,df):

	d2 = df.values
        mean=[np.mean(d2[:,0]),0,np.mean(d2[:,2]),np.mean(d2[:,3]),np.mean(d2[:,4]),np.mean(d2[:,5]),np.mean(d2[:,6]),np.mean(d2[:,7]),np.mean(d2[:,8]),np.mean(d2[:,9]),np.mean(d2[:,10]),np.mean(d2[:,11])]
        std =[np.std(d2[:,0]),1, np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9]), np.std(d2[:,10]), np.std(d2[:,11])]

	fp  = open(file0,'w')
	fp1 = open(file1,'w')
	fp2 = open(file2,'w')
	fp3 = open(file3,'w')
	fp4 = open(file4,'w')

        fp.write( "id\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("Gcnf","Cdx2","Gata6","Pbx1","Klf4","Nanog","Oct4","Oct4-Sox2","Sox2"))
        fp1.write("id\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("Gcnf","Cdx2","Gata6","Pbx1","Klf4","Nanog","Oct4","Oct4-Sox2","Sox2"))
        fp4.write("id,%s,%s,%s,%s,%s,%s,%s,%s,%s\n" %("Gcnf","Cdx2","Gata6","Pbx1","Klf4","Nanog","Oct4","Oct4-Sox2","Sox2"))
	for index,row in df.iterrows():
		temp = row[2]
		w = row[0]
		fp2.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %( temp,row[4],row[5],row[3],row[8],row[6],row[7],row[9],row[10],row[11]))

                row3  = (row[3]  -mean[3])  /std[3]
                row4  = (row[4]  -mean[4])  /std[4]
                row5  = (row[5]  -mean[5])  /std[5]
                row6  = (row[6]  -mean[6])  /std[6]
                row7  = (row[7]  -mean[7])  /std[7]
                row8  = (row[8]  -mean[8])  /std[8]
                row9  = (row[9]  -mean[9])  /std[9]
                row10 = (row[10] -mean[10]) /std[10]
                row11 = (row[11] -mean[11]) /std[11]
                #fp1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( temp,row4/w,row5/w,row3/w,row8/w,row6/w,row7/w,row9/w,row10/w,row11/w))
                fp1.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( temp,row4,row5,row3,row8,row6,row7,row9,row10,row11))
                fp4.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %( temp,row4,row5,row3,row8,row6,row7,row9,row10,row11))

                row3  = (row[3]  -mean_og[3])  /std_og[3]
                row4  = (row[4]  -mean_og[4])  /std_og[4]
                row5  = (row[5]  -mean_og[5])  /std_og[5]
                row6  = (row[6]  -mean_og[6])  /std_og[6]
                row7  = (row[7]  -mean_og[7])  /std_og[7]
                row8  = (row[8]  -mean_og[8])  /std_og[8]
                row9  = (row[9]  -mean_og[9])  /std_og[9]
                row10 = (row[10] -mean_og[10]) /std_og[10]
                row11 = (row[11] -mean_og[11]) /std_og[11]
                #fp.write( '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( temp,row4/w,row5/w,row3/w,row8/w,row6/w,row7/w,row9/w,row10/w,row11/w))
                fp.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %( temp,row4,row5,row3,row8,row6,row7,row9,row10,row11))
                fp3.write('%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n' %( temp,row4,row5,row3,row8,row6,row7,row9,row10,row11))

	fp.close()
	fp1.close()
	fp2.close()
	fp3.close()
	fp4.close()

	return [mean,std]

## read in the data  #################################
def getData(filestart,skip=False):
	print "Start"
	df1 = pd.read_csv(filepath_or_buffer=filestart+'_1.dat',header=None,sep='\t')#,names=['slen','swid','plen','pwid','target'])
	temp1= df1[0]*0.+1.
	df = pd.concat([temp1,temp1,df1[0],df1[1],df1[2],df1[3],df1[4],df1[5],df1[6],df1[7],df1[8],df1[9]],axis=1,ignore_index=True)




	try:
		di = pd	.read_csv(filepath_or_buffer=filestart+'_2.dat',header=None,sep='\t')#,names=['slen','swid','plen','pwid','target'])
		temp2= di[0]*0.+2.
		dia = pd.concat([temp2,temp2*0.+1,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
		df = pd.concat([df,dia],axis=0,ignore_index=True)
        	dia = pd.concat([temp2,temp2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
	        df = pd.concat([df,dia],axis=0,ignore_index=True)
	except:
		print "No rows in 1 and/or 2" 
	
	try:
		## solutions with 3 steady states
		di = pd.read_csv(filepath_or_buffer=filestart+'_3.dat',header=None,sep='\t')
		temp = di[0]*0.+3.
		dia = pd.concat([temp,temp*0+1.,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
		df = pd.concat([df,dia],axis=0,ignore_index=True)
                dia = pd.concat([temp,temp*0+2,di[0],di[10],di[11],di[12],di[13],di[14],di[15],di[16],di[17],di[18]],axis=1,ignore_index=True)
                df = pd.concat([df,dia],axis=0,ignore_index=True)
                dia = pd.concat([temp,temp*0+3.,di[0],di[19],di[20],di[21],di[22],di[23],di[24],di[25],di[26],di[27]],axis=1,ignore_index=True)
                df = pd.concat([df,dia],axis=0,ignore_index=True)

	except:
		print "No rows in " + filestart+"_3.dat"
	
	try:
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

	except:
		print "No rows in " + filestart+"_4.dat"


	try:
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


	except:
		print "No rows in " + filestart+"_5.dat"
	
	
	try:
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

	except:
		print "No rows in " + filestart+"_6.dat"
	
	
	if False:
		try:
			## solutions with 7 steady states
			di = pd.read_csv(filepath_or_buffer=filestart+'_7.dat',header=None,sep='\t')
			temp = di[0]*0.+7.
			dia = pd.concat([temp,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
			df = pd.concat([df,dia],axis=0,ignore_index=True)
		except:
			print "No rows in " + filestart+"_7.dat"
		
		try:
			## solutions with 8 steady states
			di = pd.read_csv(filepath_or_buffer=filestart+'_8.dat',header=None,sep='\t')
			temp = di[0]*0.+8.
			dia = pd.concat([temp,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
	                df = pd.concat([df,dia],axis=0,ignore_index=True)
		except:
			print "No rows in " + filestart+"_7.dat"

		try:
			## solutions with 9 steady states
			di = pd.read_csv(filepath_or_buffer=filestart+'_9.dat',header=None,sep='\t')
			temp = di[0]*0.+9.
			dia = pd.concat([temp,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
			df = pd.concat([df,dia],axis=0,ignore_index=True)
		except:
			print "No rows in " + filestart+"_7.dat"

		try:
			## solutions with 10 steady states
			di = pd.read_csv(filepath_or_buffer=filestart+'_10.dat',header=None,sep='\t')
			temp = di[0]*0.+10.
			dia = pd.concat([temp,di[0],di[1],di[2],di[3],di[4],di[5],di[6],di[7],di[8],di[9]],axis=1,ignore_index=True)
			df = pd.concat([df,dia],axis=0,ignore_index=True)
		except:
			print "No rows in " + filestart+"_7.dat"
	
	df=df.replace([np.inf,-np.inf],np.nan)
	df=df.dropna()
	return df

def setog(filename,df):
	d2 = df.values
        mean=[0,np.mean(d2[:,1]),np.mean(d2[:,2]),np.mean(d2[:,3]),np.mean(d2[:,4]),np.mean(d2[:,5]),np.mean(d2[:,6]),np.mean(d2[:,7]),np.mean(d2[:,8]),np.mean(d2[:,9])]
        std =[1, np.std(d2[:,1]), np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9])]
        fp2 = open(filename,'w')
        fp2.write( "id\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %("Gcnf","Cdx2","Gata6","Pbx1","Klf4","Nanog","Oct4","Oct4-Sox2","Sox2"))
        for index,row in df.iterrows():
                row1  = (row[1]  -mean[1])  /std[1]
                row2  = (row[2]  -mean[2])  /std[2]
                row3  = (row[3]  -mean[3])  /std[3]
                row4  = (row[4]  -mean[4])  /std[4]
                row5  = (row[5]  -mean[5])  /std[5]
                row6  = (row[6]  -mean[6])  /std[6]
                row7  = (row[7]  -mean[7])  /std[7]
                row8  = (row[8]  -mean[8])  /std[8]
                row9  = (row[9]  -mean[9])  /std[9]
                fp2.write('%s\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\t%3.6f\n' %( index,row1,row2,row3,row4,row5,row6,row7,row8,row9))

        fp2.close()

def setog2(filename,df):
	d2 = df.values
        mean=[0,np.mean(d2[:,1]),np.mean(d2[:,2]),np.mean(d2[:,3]),np.mean(d2[:,4]),np.mean(d2[:,5]),np.mean(d2[:,6]),np.mean(d2[:,7]),np.mean(d2[:,8]),np.mean(d2[:,9])]
        std =[1, np.std(d2[:,1]), np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9])]
        fp2 = open(filename,'w')
        for index,row in df.iterrows():
                row1  = (row[1]  -mean[1])  /std[1]
                row2  = (row[2]  -mean[2])  /std[2]
                row3  = (row[3]  -mean[3])  /std[3]
                row4  = (row[4]  -mean[4])  /std[4]
                row5  = (row[5]  -mean[5])  /std[5]
                row6  = (row[6]  -mean[6])  /std[6]
                row7  = (row[7]  -mean[7])  /std[7]
                row8  = (row[8]  -mean[8])  /std[8]
                row9  = (row[9]  -mean[9])  /std[9]
                fp2.write('%s,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f\n' %( index,row1,row2,row3,row4,row5,row6,row7,row8,row9))

        fp2.close()

def setog3(filename,df):
        fp2 = open(filename,'w')
        for index,row in df.iterrows():
                fp2.write('%s,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f\n' %( index,row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9]))

        fp2.close()



def getThresh1(filename,key):
        if key!="og":
                tdata = pd.read_csv(filepath_or_buffer=filename,header=None,sep=' ').values[0]
                return   {'Gata6':np.log2(tdata[5]),  'Gcnf':np.log2(tdata[6]),   'Cdx2':np.log2(tdata[7]), 'Klf4':np.log2(tdata[8]), 'Nanog':np.log2(tdata[9]),  'Pbx1':np.log2(tdata[10]), 'Oct4':np.log2(tdata[11]), 'Oct4-Sox2':np.log2(tdata[12]),   'Sox2':np.log2(tdata[13])}
        else:
                return   {'Gata6':np.log2(0.83),  'Gcnf':np.log2(12.4),   'Cdx2':np.log2(2.39), 'Klf4':np.log2(3.95), 'Nanog':np.log2(0.08),  'Pbx1':np.log2(38.13), 'Oct4':np.log2(2.39), 'Oct4-Sox2':np.log2(12.4),   'Sox2':np.log2(38.13)}



#########################################################
#########################################################
#########################################################
#########################################################


'''
df = getData('/home/madeline/Research/RACIPE/test_og/stem_solution_og')
d2 = df.values
'''
df = pd.read_csv(filepath_or_buffer="/home/madeline/Research/RACIPE/file1/stem_solution_0_center.dat",header=0,sep="\t")
d2 = df.values
mean_og=[0,0,0,np.mean(d2[:,1]),np.mean(d2[:,2]),np.mean(d2[:,3]),np.mean(d2[:,4]),np.mean(d2[:,5]),np.mean(d2[:,6]),np.mean(d2[:,7]),np.mean(d2[:,8]),np.mean(d2[:,9])]
std_og =[1,1,1, np.std(d2[:,1]), np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9])]
setog("clusterSC_og.dat",df)
setog("cluster_og.dat",df)
setog3("simple_og.dat",df)
setog2("simpleSC_og.dat",df)
setog2("testSC_og.dat",df)

dirs=['og','sf1','sf2','sf3']
threshs={}
for key in dirs:
	threshs[key]=getThresh1('../'+key+'/input_params_'+key+'.txt',key)

threshs['og']['Gcnf']=(threshs['og']['Gcnf']-mean_og[3])/std_og[3]
threshs['og']['Cdx2']=(threshs['og']['Cdx2']-mean_og[4])/std_og[4]
threshs['og']['Gata6']=(threshs['og']['Gata6']-mean_og[5])/std_og[5]
threshs['og']['Pbx1']=(threshs['og']['Pbx1']-mean_og[6])/std_og[6]
threshs['og']['Klf4']=(threshs['og']['Klf4']-mean_og[7])/std_og[7]
threshs['og']['Nanog']=(threshs['og']['Nanog']-mean_og[8])/std_og[8]
threshs['og']['Oct4']=(threshs['og']['Oct4']-mean_og[9])/std_og[9]
threshs['og']['Oct4-Sox2']=(threshs['og']['Oct4-Sox2']-mean_og[10])/std_og[10]
threshs['og']['Sox2']=(threshs['og']['Sox2']-mean_og[11])/std_og[11]

numDir = 3
for fname in range(1,numDir+1):
	df = getData('../sf'+str(fname)+'/stem_solution_sf'+str(fname))
	[mean,std]=main("clusterN_sf"+str(fname)+".dat","clusterSC_sf"+str(fname)+".dat",'simple_sf'+str(fname)+'.dat','simpleN_sf'+str(fname)+'.dat','simpleSC_sf'+str(fname)+'.dat',df)
	fkey = 'sf'+str(fname)
	'''
	threshs[fkey]['Gcnf']=(threshs[fkey]['Gcnf']-mean[3])/std[3]
	threshs[fkey]['Cdx2']=(threshs[fkey]['Cdx2']-mean[4])/std[4]
	threshs[fkey]['Gata6']=(threshs[fkey]['Gata6']-mean[5])/std[5]
	threshs[fkey]['Pbx1']=(threshs[fkey]['Pbx1']-mean[6])/std[6]
	threshs[fkey]['Klf4']=(threshs[fkey]['Klf4']-mean[7])/std[7]
	threshs[fkey]['Nanog']=(threshs[fkey]['Nanog']-mean[8])/std[8]
	threshs[fkey]['Oct4']=(threshs[fkey]['Oct4']-mean[9])/std[9]
	threshs[fkey]['Oct4-Sox2']=(threshs[fkey]['Oct4-Sox2']-mean[10])/std[10]
	threshs[fkey]['Sox2']=(threshs[fkey]['Sox2']-mean[11])/std[11]
	'''

df = pd.DataFrame.from_dict(threshs)
df.to_csv('modified_thresholds.dat')

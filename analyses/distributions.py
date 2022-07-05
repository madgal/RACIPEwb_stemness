import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

import matplotlib.font_manager as font_manager


fontpath = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf'

prop = font_manager.FontProperties(fname=fontpath)
plt.rcParams['font.family'] = prop.get_name()



################3
def getInputDat(filename):
	tdata = pd.read_csv(filepath_or_buffer=filename,header=None,sep=' ').values[0]
	return   [tdata[1], tdata[2],tdata[3], tdata[4]]
def getThresh1(fkey):
	df = pd.read_csv(filepath_or_buffer='modified_thresholds.dat',header=0,sep=',')
	df = df.to_dict('list')

	return {df['Unnamed: 0'][0]:df[fkey][0],df['Unnamed: 0'][1]:df[fkey][1],df['Unnamed: 0'][2]:df[fkey][2],df['Unnamed: 0'][3]:df[fkey][3],df['Unnamed: 0'][4]:df[fkey][4],df['Unnamed: 0'][5]:df[fkey][5],df['Unnamed: 0'][6]:df[fkey][6],df['Unnamed: 0'][7]:df[fkey][7],df['Unnamed: 0'][8]:df[fkey][8]}

def getThresh1a(filename,key):
	if key!="og":
		tdata = pd.read_csv(filepath_or_buffer=filename,header=None,sep=' ').values[0]
		return   {'Gata6':np.log2(tdata[5]),  'Gcnf':np.log2(tdata[6]),   'Cdx2':np.log2(tdata[7]), 'Klf4':np.log2(tdata[8]), 'Nanog':np.log2(tdata[9]),  'Pbx1':np.log2(tdata[10]), 'Oct4':np.log2(tdata[11]), 'Oct4-Sox2':np.log2(tdata[12]),   'Sox2':np.log2(tdata[13])}
	else:
		return   {'Gata6':np.log2(0.83),  'Gcnf':np.log2(12.4),   'Cdx2':np.log2(2.39), 'Klf4':np.log2(3.95), 'Nanog':np.log2(0.08),  'Pbx1':np.log2(38.13), 'Oct4':np.log2(2.39), 'Oct4-Sox2':np.log2(12.4),   'Sox2':np.log2(38.13)}

def getThresh2(filename,index):
	tdata = pd.read_csv(filepath_or_buffer=filename,header=None,sep='\t')
	tdata = tdata.sort_values(by=0)
	index = index.astype(int) -1
	tdata = tdata.values
	threshA =(tdata[index,44]+tdata[index,48]+tdata[index,60])/3.
	threshB =tdata[index,63]
	threshC =(tdata[index,49]+tdata[index,50]+tdata[index,61]+tdata[index,64])/4.
	threshD =tdata[index,59]
	threshE =(tdata[index,46]+tdata[index,51]+tdata[index,54]+tdata[index,56]+tdata[index,62])/5.
	threshF =tdata[index,57]
	threshG =(tdata[index,47]+tdata[index,52]+tdata[index,55])/3.
	threshH =(tdata[index,45]+tdata[index,58]+tdata[index,65]+tdata[index,66])/4.
	threshI =tdata[index,53]
	
	## find threshs
	return  {'Gata6':np.log2(threshA),  'Gcnf':np.log2(threshB),   'Cdx2':np.log2(threshC), 'Klf4':np.log2(threshD), 'Nanog':np.log2(threshE),  'Pbx1':np.log2(threshF), 'Oct4':np.log2(threshG), 'Oct4-Sox2':np.log2(threshH),   'Sox2':np.log2(threshI)}


#############
names = ["Gcnf","Cdx2","Gata6","Pbx1","Klf4","Nanog","Oct4","Oct4-Sox2","Sox2"]

dirs=["og","sf1","sf2","sf3"]

########################## GENE EXPRESSION - ONE FILE #####################
fig = plt.figure(figsize=(8.27, 11.69), dpi=100)
gridX,gridY=3,3
grid_size=(gridX,gridY)
rowC,colC=0,0

colors=['r','b','g','k','y','m','c','orange']
datas,threshs = [],[]
for key in dirs:
	datas+=[pd.read_csv(filepath_or_buffer='simple_'+key+'.dat',header=None,sep=',')]
	#threshs+=[getThresh1('../'+key+'/input_params_'+key+'.txt',key)]
	threshs+=[getThresh1(key)]


colC=0
for k in range(1,len(datas[0].columns)):

	plt.subplot2grid(grid_size,(rowC,colC))
	for rd in range(len(datas)):
		#print rd,colors[rd]
		#print datas[rd][k]
		plt.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
		mxx = 0.2
		#print rd,names[k-1]
		plt.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd])
		#plt.vlines(np.median(rdata[k]),0,mxx,linestyles='solid')
	plt.xlabel('Gene expression (log2)')#,fontsize=20)
	plt.ylabel('Frequency', ha='center')#,fontsize=20)
	plt.title(names[k-1])#,fontsize=30)
	plt.legend()

	rowC+=1
	if rowC==gridX:
		rowC=0
		colC+=1

#inpD = getInputDat('../s'+str(i)+'/input_params_'+str(i)+'.txt')
#title = "S"+str(i)+" B=["+str(inpD[0])+"-"+str(inpD[1])+"] U=["+str(inpD[2])+"-"+str(inpD[3])+"]"
#fig.suptitle(title,va='top')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.savefig("GeneExpression.png")
plt.close()

########################## GENE EXPRESSION - MANY FILES #####################
for k in range(1,len(datas[0].columns)):
	fig,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,figsize=(12, 8), dpi=100)
	mxx = 0.2
	for rd in range(len(datas)):
		ax1.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
		ax1.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd])
		ax1.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd])
		if rd==0:
			ax2.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
			ax2.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd],label="Threshold_estimated")
			ax2.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd],label="Median")
			ax2.legend()
		if rd==1:
			ax3.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
			ax3.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd],label="Threshold_estimated")
			ax3.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd],label="Median")
			ax3.legend()
		if rd==2:
			ax4.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
			ax4.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd],label="Threshold_estimated")
			ax4.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd],label="Median")
			ax4.legend()
		if rd==3:
			ax5.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
			ax5.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd],label="Threshold_estimated")
			ax5.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd],label="Median")
			ax5.legend()
		if rd==4:
			ax6.hist(datas[rd][k],label=dirs[rd],bins=50,alpha=0.5,density=True,color=colors[rd])
			ax6.vlines(threshs[rd][names[k-1]],0,mxx,linestyles='dashed',color=colors[rd],label="Threshold_estimated")
			ax6.vlines(np.median(datas[rd][k]),0,mxx,linestyles='solid',color = colors[rd],label="Median")
			ax6.legend()
			

	#plt.xlabel('Gene expression (log2)')#,fontsize=20)
	#plt.ylabel('Frequency', ha='center')#,fontsize=20)
	fig.text(0.5,0.04, 'Gene expression (log2)',ha='center',va='center',fontsize=15)
	fig.text(0.01,0.5, 'Frequency',ha='center',va='center',rotation=90,fontsize=15)
	plt.suptitle(names[k-1],fontsize=20)
	plt.legend()
	plt.tight_layout(rect=[0, 0.03, 1, 0.95])
	fig.savefig("pdf_s_"+names[k-1]+".png")
	plt.close()


########################## THRESHOLDS ####################################
fig = plt.figure(figsize=(5, 5), dpi=100)
above,below=[],[]
for rd in range(len(datas)):
        t1,t2=[],[] 
	for k in range(1,len(datas[rd].columns)):
		#print rd,names[k-1],np.mean(datas[rd][k])
		temp= np.mean(datas[rd][k] > threshs[rd][names[k-1]])
		t1+=[temp]
		temp= np.mean(datas[rd][k] < threshs[rd][names[k-1]])
		t2+=[temp]
	above  += [t1]
	below  += [t2]


#print above
#print below
titles=["RACIPE","RACIPE-wb (par1)","RACIPE-wb (par2)","RACIPE-wb (par3)"]
#print len(datas)
for rd in range(len(datas)):
	fig = plt.figure()
	plt.bar(names,np.array(above[rd])*0.+1.,color='orange',label="x/x0>1")
        plt.bar(names,np.array(below[rd]),color='g',label="x/x0<1")
	plt.axhline(y=0.5, color='r',linestyle='--')
	plt.ylim([0,1])
	plt.legend(loc='upper right',fontsize=15)
	plt.xticks(rotation=60)
	ax = plt.gca()
	ax.tick_params(axis='both',labelsize=15)
	titl = titles[rd]
	#print titl
	plt.title(titl)
	plt.ylabel('Probability',fontsize=15)
	plt.xlabel('Gene',fontsize=15)
	plt.tight_layout()
	fig.savefig('threshold_'+dirs[rd]+'.png')
	plt.close()

#PCAresults_s1.dat

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
from matplotlib.backends.backend_pdf import PdfPages

import matplotlib.font_manager as font_manager


fontpath = '/usr/share/fonts/truetype/freefont/FreeSerif.ttf'

prop = font_manager.FontProperties(fname=fontpath)
plt.rcParams['font.family'] = prop.get_name()

import sys,os


def plot(data,xmin1,xmax1,ymin1,ymax1,cmax1,sett=False):
        [HH,xh,yh] = np.histogram2d(data[:,2],data[:,3],bins=(100,100),density=False)
        HH=HH.T

        if sett:
                xmin1, xmax1 = xh[0], xh[-1]
                ymin1, ymax1 = yh[0], yh[-1]
                cmax1 = np.max(HH)

	xmin1,xmax1=-4,4
	ymin1,ymax1=-3,3
	print HH
        plt.imshow(HH,interpolation='nearest',origin='low',extent=[xmin1,xmax1,ymin1,ymax1],cmap='jet',vmin=0., vmax=7.5,aspect='auto')#cmax1)
	ax = plt.gca()
	ax.xaxis.set_ticks([-4,-2,0,2,4])
	ax.tick_params(axis='both', labelsize=25)

        #plt.colorbar()

        return[xmin1,xmax1,ymin1,ymax1,cmax1]

def plot2(data,xmin1,xmax1,ymin1,ymax1,name):
        #colors=['k','b','g','r','c','m','y','plum','navy','line']
	colors=["#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#A4E804", "#004D43", "#8FB0FF", "#3B5DFF", "#B903AA", "#BA0900", "#FFB500", "#6B7900", "#FF90C9", "#BC23FF", "#99ADC0"]

        maxi=0
        cols=[]
        for i in data[:,0]:
                i = int(i)
                cols+=[colors[i%len(colors)]]
                if i> maxi:
                        maxi=i

        for i in range(maxi):
                plt.scatter(-5,-5,color=colors[(i+1)%len(colors)],label='C'+str(i+1),alpha=0.6)
	ax = plt.gca()
	ax.xaxis.set_ticks([-4,-2,0,2,4])
	ax.yaxis.set_ticks([-3,-2,-1,0,1,2,3])
	ax.tick_params(axis='both', labelsize=25)
        #plt.scatter(data[1,2],data[1,3],color='y',label='C4')
        plt.scatter(data[:,2],data[:,3],color=cols,alpha=0.2)
        plt.xlim(xmin1,xmax1)
        plt.ylim(ymin1,ymax1)
	plt.gca().set_aspect('auto', adjustable='box')
	if 'SC' in name:
	        plt.legend(bbox_to_anchor=(1.04,0.8),fontsize=10)
	else:
	        plt.legend(bbox_to_anchor=(1.06,0.45),fontsize=10)

xmin1,xmax1,ymin1,ymax1,cmax1=0.,0.,0.,0.,0.
getSt=True
################################ GET PCA RESULTS - SC ################################
for file in os.listdir("."):
	if "PCAresults" in file:
                fig = plt.figure(figsize=(7.5, 6), dpi=400)
		data = pd.read_csv(filepath_or_buffer=file,header=None,sep=",").values
		
		p1,p2 = data[0][0],data[0][1]
		data = data[1:]

		if getSt:
			[xmin1,xmax1,ymin1,ymax1,cmax1]=plot(data,xmin1,xmax1,ymin1,ymax1,cmax1,True)
			getSt=False
		else:
			plot(data,xmin1,xmax1,ymin1,ymax1,cmax1)

		cbar = plt.colorbar()#fraction = 0.030,pad=0.04)
		cbar.set_label(label='Probability density', size=25,rotation=270,labelpad=30)
		ticks=[0,3.75,7.50]
		cbar.set_ticks(ticks)
		cbar.ax.set_yticklabels(ticks)	
		cbar.ax.tick_params(labelsize=25)


		xlab = "PC1 ({:2.2f}%)".format(p1)
		ylab = "PC2 ({:2.2f}%)".format(p2)
		plt.xlabel(xlab,size=25)
		plt.ylabel(ylab,size=25)
		#plt.title("Probability density of PCA")

		title = file.replace("PCAresults","pca")
		title = title.replace("dat","png")
		plt.tight_layout()
		fig.savefig(title)	
		plt.close()

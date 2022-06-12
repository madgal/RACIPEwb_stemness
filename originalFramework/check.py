import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

d1 = pd.read_csv(filepath_or_buffer = "/home/madeline/Research/RACIPE/file1/stem_solution_0_zscore.dat",header=0,sep="\t")
d2 = pd.read_csv(filepath_or_buffer = "/home/madeline/Research/RACIPE/latest_model/full_sims/analyses/clusterSC_og.dat",header=0,sep="\t")
d3 = pd.read_csv(filepath_or_buffer = "/home/madeline/Research/RACIPE/latest_model/full_sims/analyses/simpleSC_og.dat",header=None,sep=",")

d1 = d1.values[:,1:]
d2 = d2.values[:,1:]
d3 = d3.values[:,1:]


h1 = np.histogram(d1[:,0])
h2 = np.histogram(d2[:,0])
#h3 = np.histogram(d3[:,0])

plt.bar(h1[1][:-1],h1[0],label='1',alpha=0.7)
plt.bar(h2[1][:-1],h2[0],label='2',alpha=0.7)
#plt.bar(h3[1][:-1],h3[0],label='3',alpha=0.7)
plt.show()

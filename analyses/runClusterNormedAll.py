## The parallel part start
from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()
stat = MPI.Status()
import os,sys


from Bio import Cluster

def main(filein,fileout):
	with open(filein) as handle:
		record = Cluster.read(handle)

	genetree = record.treecluster(method='a',dist='e')
	genetree.scale()
	record.save(fileout,genetree)

allfiles=[]
for file in os.listdir("."):
	if ("cluster" in file) and ("dat" in file):
		allfiles+=[file]


## Determine which files it will process
numF =len(allfiles)
if numF >= nprocs:
	startV = my_rank*(numF/nprocs)
	stopV = (my_rank+1)*(numF/nprocs)
	if(my_rank==nprocs-1):
		stopV = numF
else:
	if my_rank<numF:
		startV = my_rank
		stopV = my_rank+1
	else:
		exit()

for i in range(startV,stopV):
	filein = allfiles[i]
	fileout = filein.replace("cluster","result")
	fileout = fileout.replace(".dat","")

	mystr = filein +" " + fileout
	print mystr
	main(filein,fileout)

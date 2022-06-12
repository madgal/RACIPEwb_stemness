import numpy as np
import pandas as pd



def main(file2,df):
	#df.replace(0, np.nan, inplace=True)
	#df.replace(np.inf, np.nan, inplace=True)
	#df.replace(-np.inf, np.nan, inplace=True)
	#df = df.dropna(axis=0)

        #d2 = df.values

	#mean =[(d2[1,0]), np.mean(d2[:,1]), np.mean(d2[:,2]), np.mean(d2[:,3]), np.mean(d2[:,4]), np.mean(d2[:,5]), np.mean(d2[:,6]), np.mean(d2[:,7]), np.mean(d2[:,8]), np.mean(d2[:,9])]
	#std =[(d2[1,0]), np.std(d2[:,1]), np.std(d2[:,2]), np.std(d2[:,3]), np.std(d2[:,4]), np.std(d2[:,5]), np.std(d2[:,6]), np.std(d2[:,7]), np.std(d2[:,8]), np.std(d2[:,9])]


        fp2 = open(file2,'w')

        for index,row in df.iterrows():
                temp = row[0]
		'''
		row[1] = ((row[1]) -mean[1])/std[1]
		row[2] = ((row[2]) -mean[2])/std[2]
		row[3] = ((row[3]) -mean[3])/std[3]
		row[4] = ((row[4]) -mean[4])/std[4]
		row[5] = ((row[5]) -mean[5])/std[5]
		row[6] = ((row[6]) -mean[6])/std[6]
		row[7] = ((row[7]) -mean[7])/std[7]
		row[8] = ((row[8]) -mean[8])/std[8]
		row[9] = ((row[9]) -mean[9])/std[9]
		'''
                fp2.write('%s,%s,%s,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f,%3.6f\n' %( index,index,index,row[1],row[2],row[3],row[4],row[5],row[6],row[7],row[8],row[9]))

        fp2.close()



#df = pd.read_csv(filepath_or_buffer="/home/madeline/Research/RACIPE/file1/stem_solution_0_raw.dat",header=0,sep="\t")
df = pd.read_csv(filepath_or_buffer="/home/madeline/Research/RACIPE/file1/stem_solution_0_zscore.dat",header=0,sep="\t")

main("testcPCA.dat",df)

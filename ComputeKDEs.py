##################################
##Compute kernel density estimates
##################################

##################################
##load functions

from scipy import stats
import pickle
import random
import numpy as np

##################################
##Set values to cycle through

num=10000
pnum=99
snum=59
ps=np.linspace(0.01,0.99,pnum)
ss=np.linspace(-12.0,17.0,snum)

##################################
##Populate qs

qs=np.empty((snum,pnum,num))
samp=np.loadtxt('qSamples/final_samples_cuda.txt')
for i in range(snum):
    for j in range(pnum):
        for k in range(num):
            qs[i,j,k]=samp[k+j*num+i*pnum*num]

##################################
##Save organized qs

#file1 = open('qSamples/qarray_cuda','w')
#pickle.dump(qs,file1)
#file1.close()   

##################################
##Compute kernel density estimates

kdes=np.empty((snum,pnum),dtype=object)
#[stats.gaussian_kde(qs[range(np.shape(qs)[0]),0])]
for i in range(snum):
    for j in range(pnum):
        qij=qs[i,j,range(num)]
        kdes[i,j]=stats.gaussian_kde(qij[qij<=1])

##################################
##Store 1001 values for each KDE

kdegrid=np.empty((np.shape(qs)[0],np.shape(qs)[1],1001))
for i in range(np.shape(qs)[0]):
    for j in range(np.shape(qs)[1]):
        kdegrid[i,j,0:1001]=kdes[i,j].evaluate(np.linspace(0,1,1001))

##################################
##Save kernel density estimates

file2 = open('DensityEstimates/kdes_fly1217','w')
pickle.dump(kdegrid,file2)
file2.close()
###########################
##MCMC for hypoxic fly data
###########################

###load functions###
execfile("WFwithSelection-ApproxAll.py")
from scipy import stats
import pickle
import random

###load kernel density estimates###
file1 = open('DensityEstimates/kdes_fly1217','r')
kdeg = pickle.load(file1)
file1.close()

###load neutral data###
flyn=np.loadtxt('Data/fly-n.txt')

###load data with selection###
flies=np.loadtxt('Data/fly-s.txt')
n=200 
Y=map(round,flies*n)

###dimensions###
snum=np.shape(kdeg)[0]
pnum=np.shape(kdeg)[1]
Ngrid=np.shape(kdeg)[2]

###rescale kernel density estimates###
for i in range(snum):
    for j in range(pnum):
        kdeg[i,j,range(Ngrid)]=kdeg[i,j,range(Ngrid)]/sum(kdeg[i,j,range(Ngrid)])

###empirical density for p###
pdensity=stats.gaussian_kde(flyn)
discretep=np.empty(99)
for i in range(99):
    discretep[i]=pdensity.evaluate(i/100.)
discretep=discretep/sum(discretep)

###integrate out p###
pout=np.empty((snum,Ngrid))
for j in range(snum):
    for k in range(Ngrid):
        pout[j,k]=sum(discretep*kdeg[j,range(pnum),k])

###functions for MCMC###

##interpolation of kde
def KernDenEst(q,sind):
    qi=math.modf(1000*q)
    below=pout[sind,qi[1]]
    above=pout[sind,qi[1]+1]
    return below+qi[0]*(above-below)

##likelihood functions
def logqgivesK(qs,sind):
    return sum(np.log(map(KernDenEst,qs,np.array(sind).repeat(len(qs)))))

def logqgives1(q,sind):
    return np.log(KernDenEst(q,sind))

##prior distributions
###prior on s is discrete uniform

##acceptance ratios
def ARs(ys,n,sicur,sinew,qs):
    return min(1,math.exp(logqgivesK(qs,sinew)-logqgivesK(qs,sicur)))

def ARq(y,n,qcur,qnew,sind):
    return min(1,math.exp(stats.binom.logpmf(y,n,qnew)+logqgives1(qnew,sind)-stats.binom.logpmf(y,n,qcur)-logqgives1(qcur,sind)))

###set up for MCMC###

##model parameters
steps=10000				#steps per chain
burn=1000		        	#burn-in steps per chain	
K=np.shape(flies)[0]			#number of SNPs
theta=.00014				#first mutation parameter
T=.1					#length of time
svalid=np.linspace(-12.0,17.0,snum)		#grid of possible values of s
sdq=0.05				#std. deviation of proposal density for each q_k

##create arrays to store parameter values
sinds=np.zeros(steps)	#indices for selection parameter
sinds=sinds.astype(int)	#treat indices as integers
ss=np.ndarray((steps,1))	#selection parameter
qs=np.empty((steps,K))	#q for each SNP

##initial values for parameters
sinds[0]=24		#index for selection parameter
ss[0]=svalid[sinds[0]]	#selection parameter
qs[0]=0.5		#q for each SNP
rejectS=0		#count of rejections of proposed s
rejectQ=np.zeros(K)	#count of rejections of proposed q_k

##run the MCMC
for j in range(1,steps):
    
    #propose new s given the q_k
    sinds[j]=np.random.random_integers(sinds[j-1]-5,sinds[j-1]+5)   
    if sinds[j]>snum-1 or sinds[j]<0:
        ss[j]=ss[j-1]
        sinds[j]=sinds[j-1]
        rejectS+=1
    else:
        ss[j]=svalid[sinds[j]]
        sU=np.random.uniform()
        salpha=ARs(Y,n,sinds[j-1],sinds[j],qs[j-1])
        if sU>salpha:
            ss[j]=ss[j-1]
            sinds[j]=sinds[j-1]
            rejectS+=1
        
    #propose new q_k given s
    for k in range(K):   
        qs[j,k]=np.random.normal(qs[j-1,k],sdq,1)
        if qs[j,k]>1.0 or qs[j,k]<0.0:
            qs[j,k]=qs[j-1,k]
            rejectQ[k]+=1
        else:
            qU=np.random.uniform()
            qalpha=ARq(Y[k],n,qs[j-1,k],qs[j,k],sinds[j])
            if qU>qalpha:
                qs[j,k]=qs[j-1,k]
                rejectQ[k]+=1            

##save the parameter values in the chain to a file
np.savetxt('Chains/MCMC_fly1s_output.txt',np.hstack((ss,qs)))

##save rejection counts
np.savetxt('Output/MCMC_fly1s_reject.txt',np.hstack((rejectS,rejectQ)))
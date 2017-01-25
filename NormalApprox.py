import numpy as np
import math

def Beta(t,theta):
    return 0.5*(theta-1)*t

def Eta(beta):
    return beta/(math.exp(beta)-1)

def VarM(t,beta,eta):
    return abs(2.0*eta/t*(eta+beta)**2*(1+eta/(eta+beta)-2*eta))
    
def NormApprox(t,theta):
    beta=Beta(t,theta)
    if beta==0:
        return np.random.normal(2.0/t,math.sqrt(2.0/3/t))
    else:
        eta=Eta(beta)
	if VarM(t,beta,eta)<0.0000000001:
            return 2.0*eta/t
        return np.random.normal(2.0*eta/t,math.sqrt(VarM(t,beta,eta)))

import numpy as np
import math

execfile("NormalApprox.py")

##Coefficients in expression for q_m
def a(k,m,theta):
    return (theta+2*k-1)*math.gamma(theta+m+k-1)/math.gamma(theta+m)/math.factorial(m)/math.factorial(k-m)

##Coefficients in expression for q_m using log gamma
def loga(k,m,theta):
    return math.log(theta+2*k-1)+math.lgamma(theta+m+k-1)-math.lgamma(theta+m)-math.lgamma(m+1)-math.lgamma(k-m+1)

##Terms in the alternating series
def b(k,m,t,theta):
    if (m > 0 or k > 0):
        return math.exp(loga(k,m,theta)-k*(k+theta-1)*t/2)
    else:
        return a(k,m,theta)*math.exp(-k*(k+theta-1)*t/2)

##Step at which terms in series start decreasing
def C(m,t,theta):
    i=0.
    while b(i+m+1,m,t,theta) >= b(i+m,m,t,theta):
        i+=1
    return i

##Terms in series
def Sterm(m,i,t,theta):
    return (-1)**i*b(m+i,m,t,theta)

##Algorithm to sample from transition density
def Algorithm2(t,theta1,theta2,x):
    m=round(NormApprox(t,theta1+theta2))
    L=np.random.binomial(m,x)
    Y=np.random.beta(theta1+L,theta2+m-L)
    return Y

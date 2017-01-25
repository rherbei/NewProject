execfile("NeutralWF-NormApproxAll.py")

##Simulate a Poisson Point Process
def PoissonPtProc(T,phi,rate=1):
    N=np.random.poisson(rate*T*phi)
    ts=np.random.uniform(0,T,N)
    psis=np.random.uniform(0,phi,N)
    return ts, psis, N

##These functions are specific to a Wright-Fisher diffusion with
##drift sigma*x*(1-x)+0.5*(theta1*(1-x)-theta2*x))
def phitilde(x,theta1,theta2,sigma):
    return sigma/2*(-sigma*x**2+x*(sigma-theta1-theta2)+theta1)

def Atilde(x,sigma):
    return sigma*x

def phiplus(sigma,theta1,theta2):
    if 0<=0.5/sigma*(sigma-theta1-theta2)<=1:
        phimax=1./8*(sigma-theta1-theta2)**2+theta1*sigma/2
    elif 0.5/sigma*(sigma-theta1-theta2)<0:
        phimax=theta1*sigma/2
    else:
        phimax=-theta2*sigma/2
    return phimax

def Aplus(sigma):
    if sigma >= 0:
        return sigma
    else:
        return 0

##Exact simulation of Wright-Fisher diffusion with selection and mutation
def Algorithm6(T,theta1,theta2,sigma,init):
    while True:
        reject=0
        ppp=PoissonPtProc(T,phiplus(sigma,theta1,theta2),rate=1)
        U=np.random.uniform()
        J=ppp[2]
        X=np.zeros(J+2)
        X[0]=init
        sort=np.argsort(ppp[0])
        t=np.hstack((0,ppp[0][sort]))
        psi=ppp[1][sort]
        for j in range(J):
            X[j+1]=Algorithm2(t[j+1]-t[j],theta1,theta2,X[j])
            if phitilde(X[j+1],theta1,theta2,sigma)>psi[j]:
                reject=1
                break
        if reject==1:
            continue
        X[J+1]=Algorithm2(T-t[-1],theta1,theta2,X[J])
        if U<=math.exp(Atilde(X[-1],sigma)-Aplus(sigma)):
            return X[J+1]

import numpy as np
from scipy import special
N,p=150,.5
def f(n):
    return special.binom(N,n)*p**n*(1-p)**(N-n)
s = zip(range(1,N),map(f,range(1,N)))
np.savetxt("bernoulli.dat",s,fmt='%0.5f')


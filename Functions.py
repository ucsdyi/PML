import numpy as np
import scipy.stats
log = np.log
exp = np.exp

def entropy(distribution):
    return scipy.stats.entropy(distribution)

def renyi_entropy(distribution,alpha):
    return log(sum(distribution**alpha))/(1-alpha)

def support_size(distribution,minimum_probability):
    return sum((distribution>=minimum_probability)*1)

def support_coverage(distribution,m):
    return sum(1-exp(-m*distribution))

def distance_to_uniformity(distribution,k):
    return sum(abs(distribution-1.0/k))

def sorted_L1_distance(distribution,true_distribution):
    k = max(len(distribution), len(true_distribution))
    distribution = np.pad(distribution,(0, k-len(distribution)),'constant')
    true_distribution = np.pad(true_distribution,(0, k-len(true_distribution)),'constant')
    return sum(abs(np.sort(distribution)-np.sort(true_distribution)))

# bound support size
# implementation based on the paper "Orlitsky, A., Suresh, A. T., & Wu, Y. (2016). Optimal prediction of the number of unseen species. Proceedings of the National Academy of Sciences, 113(47), 13283-13288."
def binomialcoefficient(n,p,i):
    output = (n-i)*log(1-p)
    for j in range(0,i):
        output += log((n-j)/(j+1))
    output += i*log(p)
    return exp(output)

def binomialtail(n,p,i):
    output = 0
    for j in range(i,n+1):
        output += binomialcoefficient(n,p,j)
    return output

def et_estimator(t,n,r,q,maxfreq):
    if (t > 0 and n > 0):
        logt = log(t)
        running = 0
        estimator = [0] * (maxfreq+1)
        for i in range(0,min(r+1, maxfreq+1)):
            running += logt
            if binomialtail(r,q,i+1) > 0:
                estimator[i] = ((-1)**(i))*exp(running + log(binomialtail(r,q,i+1)))
    else:
        estimator = [0] * (maxfreq+1)
    return estimator

def linear_estimator(profile,estimator,n,t):
    output = 0
    maxfreq = len(profile)
    for i in range(0,maxfreq):
        output += profile[i]*estimator[i]
    return float(min(max(output,0),n*t))

def bound_support(profile,t,n):
    tnew = max(t,1.001)
    binoptinot = int((log(n*tnew*tnew/(tnew-1))/(2*log(3))))
    maxfreq = len(profile)
    U = linear_estimator(profile, et_estimator(t,n,binoptinot, 1/float(t+1),maxfreq), n, t)
    S = sum(profile)
    return S+U

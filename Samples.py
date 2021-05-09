import numpy as np
from scipy.stats import nbinom
from scipy.stats import binom
from scipy.stats import logser
from scipy.stats import boltzmann

def get_samples(distribution, n):
    values = np.random.rand(n)
    values.sort()
    cumulative = np.cumsum(distribution)
    z = 0
    k = len(distribution)
    freq = [0] * k
    for x in values:
        while x > cumulative[z]:
            z+=1
        freq[z]+=1
    
    freq_positive = []
    D_freq_positive = {}
    for j in range(k):
        x = freq[j]
        if x>0:
            freq_positive.append(x)
            if x in D_freq_positive:
                D_freq_positive[x].append(j)
            else:
                D_freq_positive[x] = [j]

    return [freq_positive, D_freq_positive]

def get_distribution(dist_name, k):
    k = int(k)
    if dist_name == 'Uniform':
        raw_distribution = [1] * k
    elif dist_name == 'Two-step':
        raw_distribution = [1] * (k//2) + [4] * (k-k//2)
    elif dist_name == 'Three-step':
        raw_distribution = [1] * (k//3) + [3] * (k//3)+[9]*(k-2*k//3)
    elif dist_name == 'Subset-uniform':
        raw_distribution = [1] * (k//2) + [0] * (k-k//2)
    elif dist_name == 'Log-series':
        p = (k-2)/k
        raw_distribution = [logser.pmf(j, p) for j in range(1,k+1)]
    elif dist_name == 'Geometric':
        p = (k-1)*1.0/k
        raw_distribution = [(1-p)*p**i for i in range(k)]
    elif dist_name == 'Zipf-half':
        raw_distribution = [1/(float(i)**(0.5)) for i in range(1,k+1)]
    elif dist_name == 'Zipf-one':
        raw_distribution = [1/(float(i)**(1)) for i in range(1,k+1)]

    return raw_distribution

import numpy as np
import matplotlib.pyplot as plt
import ctypes
import platform
from ctypes import cdll
from collections import Counter
from Samples import *
from Functions import *

log = np.log
exp = np.exp

def generate_distribution(distribution_name,k):
    raw_distribution = get_distribution(distribution_name,k)
    sum_raw = sum(raw_distribution)
    distribution = sorted([float(y)/float(sum_raw) for y in raw_distribution])
    
    return np.array(distribution)

def generate_sample(distribution,n):
    [multiplicity, D_multiplicity] = get_samples(distribution,n)
    
    return [np.array(multiplicity), D_multiplicity]

def compute_empirical(multiplicity):
    n = np.sum(multiplicity)*1.0
    
    return np.array(multiplicity)/n

def D_compute_empirical(D_multiplicity):
    n = 0
    Key = D_multiplicity.keys()
    for j in Key:
        n += j*len(D_multiplicity[j])
    D_labeled_empirical = {}
    for j in Key:
        symbol_list = D_multiplicity[j]
        for s in symbol_list:
            D_labeled_empirical[s] = j/n

    return D_labeled_empirical

def compute_PML(multiplicity,max_EM_iter=30,support_bound=None):
    proFile_name = 'proFile'
    PMLFile_name = 'PMLFile'
    
    n = np.sum(multiplicity)
    multiplicity_large = []
    multiplicity_small = []
    for i in range(len(multiplicity)):
        v = multiplicity[i]
        if v> int(1.5*log(n)**2): #n/100: int(1.5*log(n)): #
            multiplicity_large.append(v)
        else:
            multiplicity_small.append(v)
    profile_counter = Counter(multiplicity_small)
    profile = [profile_counter[i] for i in range(1,max(multiplicity_small)+1)] #prevalences of positive integers
    
    with open(proFile_name,'w') as file:
        for prevalence in profile:
            file.write("%i " % prevalence)

    #approximately compute a PML distribution
    if support_bound == None:
        t = int(log(n))
        support_bound = int(bound_support(profile,t,n))
    if platform.system()== 'Windows':
        lib = cdll.LoadLibrary('./MCMCEM.dll')
    else:
        lib = cdll.LoadLibrary('./MCMCEM.so')
    lib.PMLplus(support_bound,max_EM_iter,int(n))

    PML_distribution = []
    with open(PMLFile_name,'r') as PML_file:
        for probability in PML_file:
            PML_distribution.append(float(probability.strip('\n')))

    probability_large = np.array(multiplicity_large)/n
    probability_small = np.array(PML_distribution)
    PML_distribution = np.concatenate([probability_small*(1-np.sum(probability_large)),probability_large])
    
    return PML_distribution

def Good_Turing(D_multiplicity, profile, k):
    n = 0
    Key = D_multiplicity.keys()
    D_labeled_distribution = {}
    total = 0
    Phi0 = 0
    for j in Key:
        Phi0 += len(D_multiplicity[j])
    Phi0 = k - Phi0
    for j in Key:
        symbol_list = D_multiplicity[j]
        if j<= len(profile)-1:
            if j>profile[j]:
                total += j*len(symbol_list)
                for s in symbol_list:
                    D_labeled_distribution[s] = j
            else:
                v = profile[j]/profile[j-1]*(j+1)
                total += v*len(symbol_list)
                for s in symbol_list:
                    D_labeled_distribution[s] = v
        else:
            total += j*len(symbol_list)
            for s in symbol_list:
                D_labeled_distribution[s] = j
    if Phi0>0:
        v = profile[0]/Phi0
        for s in range(k):
            if s not in D_labeled_distribution:
                D_labeled_distribution[s] = v
        total += profile[0]
    for s in D_labeled_distribution.keys():
        D_labeled_distribution[s] /= total

    return [v/total,D_labeled_distribution]

def weighted_median(real_list,n,j):
    #assume that real_list is sorted
    weight_list = [binomialcoefficient(n,x,j) for x in real_list] #iteration???
    total_weight = sum(weight_list)
    left_sum, right_sum = 0, total_weight
    for j in range(len(real_list)):
        right_sum = right_sum - weight_list[j]
        if right_sum>total_weight/2 or left_sum>total_weight/2:
            left_sum += weight_list[j]
        else:
            return real_list[j]

def to_labled(unlabled_distribution, D_multiplicity,k,v1):
    n = 0
    Key = D_multiplicity.keys()
    for j in Key:
        n += j*len(D_multiplicity[j])
    Additional = []

    for j in range(1,int(log(n)**2)):
        Additional += int(n/(j*log(n)**4))*[j/n]

    Additional = np.array(Additional)
    unlabled_distribution = np.sort(np.concatenate([unlabled_distribution*(1-np.sum(Additional)),Additional]))
    D_labeled_distribution = {}
    for j in Key:
        symbol_list = D_multiplicity[j]
        if j<log(n)**2:
            m = weighted_median(unlabled_distribution,n,j)
            for s in symbol_list:
                D_labeled_distribution[s] = m
        else:
            v = j/n
            for s in symbol_list:
                D_labeled_distribution[s] = v

    for s in range(k):
        if s not in D_labeled_distribution:
            D_labeled_distribution[s] = v1
    return D_labeled_distribution

def labeled_L1(D_dist_estimate,distribution):
    L1 = 0
    for j in range(distribution.size):
        if j in D_dist_estimate:
            L1 += abs(D_dist_estimate[j]-distribution[j])
        else:
            L1 += abs(distribution[j])

    return L1

def functions(distribution,function_name='entropy',para=None):
    if function_name == 'entropy':
        return entropy(distribution)
    elif function_name == 'support_size': #para = minimum probability
        return support_size(distribution,para)
    elif function_name == 'support_coverage': #para = m
        return support_coverage(distribution,para)
    elif function_name == 'distance_to_uniformity': #para = k
        return distance_to_uniformity(distribution,para)
    elif function_name == 'renyi_entropy': #para = alpha
        return renyi_entropy(distribution,para)
    elif function_name == 'sorted_L1_distance': #para = true_distribution
        return sorted_L1_distance(distribution,para)

def experiment(function_name,distribution_name):
    #r controls the number of MC simulations, and max_EM_iter reflects the number of EM iterations
    r, max_EM_iter = 1, 30
    k = 5000
    num = 20 #number of different sample sizes
    nlist = range(2*k,20*k+1,5000) #at least 4 values
    m, alpha = 2*k, 0.5 #m and alpha are parameters for support coverage and Renyi entropy, respectively
    true, empirical, empirical_nlogn, PML, GT = [], [], [], [], []
    distribution = generate_distribution(distribution_name,k)
    support_bound = None
    
    if function_name == 'support_size': para = min(distribution)
    elif function_name == 'distance_to_uniformity': para, support_bound = k, k
    elif function_name == 'renyi_entropy': para = alpha
    elif function_name == 'support_coverage': para = m
    elif function_name == 'sorted_L1_distance': para = distribution
    else: para = None

    true_value = functions(distribution,function_name,para)
    for n in nlist:
        emp_est, emp_nlogn_est, PML_est, GT_est = 0, 0, 0, 0
        for j in range(r):
            [multiplicity, D_multiplicity] = generate_sample(distribution,n)
            empirical_distribution = compute_empirical(multiplicity)
            m = int(n*log(n))
            [multiplicity_m, D_multiplicity_m] = generate_sample(distribution,m)
            empirical_dist_nlogn = compute_empirical(multiplicity_m)
            PML_distribution = compute_PML(multiplicity,max_EM_iter,support_bound)
            if function_name == 'L1_distance':
                profile_counter = Counter(multiplicity)
                profile = [profile_counter[i] for i in range(1,max(multiplicity)+1)]
                [v,labeled_GT] = Good_Turing(D_multiplicity,profile,k)
                labeled_empirical = D_compute_empirical(D_multiplicity)
                labeled_empirical_nlogn = D_compute_empirical(D_multiplicity_m)
                labeled_PML = to_labled(PML_distribution,D_multiplicity,k,v)
            if function_name != 'L1_distance':
                emp_est += abs(functions(empirical_distribution,function_name,para)-true_value)
                emp_nlogn_est += abs(functions(empirical_dist_nlogn,function_name,para)-true_value)
                PML_est += abs(functions(PML_distribution,function_name,para)-true_value)
            else:
                GT_est += labeled_L1(labeled_GT,distribution)
                emp_est += labeled_L1(labeled_empirical,distribution)
                emp_nlogn_est += labeled_L1(labeled_empirical_nlogn,distribution)
                PML_est += labeled_L1(labeled_PML,distribution)

        if function_name != 'L1_distance':
             true.append(functions(distribution,function_name,para))
        else:
             true.append(0.0)
             GT.append(GT_est/r)

        empirical.append(emp_est/r)
        empirical_nlogn.append(emp_nlogn_est/r)
        PML.append(PML_est/r)
       
    print("\n Done.\n")

    fig = plt.figure()
    xnew = np.linspace(min(nlist),max(nlist),200,endpoint=True)

    legend = []
    if function_name == 'L1_distance':
        plt.plot(nlist,GT,marker='<',color='darkviolet',linestyle = '-',markerfacecolor='None',markersize=6)
        legend += ['Good-Turing']
    plt.plot(nlist,PML, marker='D',color='c',linestyle = '-',markerfacecolor='None',markersize=6)
    legend += ['PML']
    plt.plot(nlist,empirical,marker='s',color='gold',linestyle = '-',markerfacecolor='None',markersize=6)
    legend += ['Empirical']
    plt.plot(nlist,empirical_nlogn,marker='>',color='orangered',linestyle = '-',markerfacecolor='None',markersize=6)
    legend += ['Empirical-nlogn']

    plt.legend(legend, fontsize=18, loc='upper right')

    if function_name == 'renyi_entropy':
        function_name  = str(alpha)+'-'+function_name

    plt.title(distribution_name.replace("_", " ")+' distribution', fontsize=25)
    plt.xlabel('Sample size (n)', fontsize=16)
    plt.ylabel('Average error', fontsize=16)
    save_name = function_name+'_'+distribution_name+'.jpg'
    fig.savefig(save_name)
    #plt.show()

def main():
    #Options for distribution_name: 'Uniform','Two-step','Three-step','Geometric','Zipf-one','Zipf-half','Log-series'
    #Options for function_name: 'L1_distance', 'support_coverage','support_size','entropy', 'renyi_entropy','distance_to_uniformity','sorted_L1_distance',
    D_list = ['Two-step']
    for function_name in ['L1_distance']:
        for distribution_name in D_list:
            print("Performing experiments for "+function_name+" and "+distribution_name+".")
            experiment(function_name,distribution_name)

if __name__=="__main__":
    main()

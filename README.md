# Code for the paper "The Broad Optimality of Profile Maximum Likelihood" by Hao & Orlitsky (NeurIPS 2019)
Paper is available at https://proceedings.neurips.cc/paper/2019/hash/f9fd5ec4c141a95257aa99ef1b590672-Abstract.html

# Abstract
We study three fundamental statistical-learning problems: distribution estimation, property estimation, and property testing. We establish the profile maximum likelihood (PML) estimator as the first unified sample-optimal approach to a wide range of learning tasks. In particular, for every alphabet size k and desired accuracy ε:

"Distribution estimation"  Under ℓ1 distance, PML yields optimal Θ(k/(ε2logk)) sample complexity for sorted-distribution estimation, and a PML-based estimator empirically outperforms the Good-Turing estimator on the actual distribution;

"Additive property estimation"  For a broad class of additive properties, the PML plug-in estimator uses just four times the sample size required by the best estimator to achieve roughly twice its error, with exponentially higher confidence;

"α-Rényi entropy estimation"  For integer α>1, the PML plug-in estimator has optimal k1−1/α sample complexity; for non-integer α>3/4, the PML plug-in estimator has sample complexity lower than the state of the art;

"Identity testing"  In testing whether an unknown distribution is equal to or at least ε far from a given distribution in ℓ1 distance, a PML-based tester achieves the optimal sample complexity up to logarithmic factors of k.

With minor modifications, most of these results also hold for a near-linear-time computable variant of PML.


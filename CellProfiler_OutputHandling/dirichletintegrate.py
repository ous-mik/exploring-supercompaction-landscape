from numpy import *
from scipy.integrate import quadrature, romberg, fixed_quad
from scipy.special import gammaln, betaln, digamma, polygamma, betainc, gamma
import pdb
from hypergeom import hyper3F2Z1, hyper3F2aZ1
# import polyafit

# This file is copied from the source code of CellProfiler Analyst
# Copyright (c) Broad Institute
# Licensed under the BSD 3-Clause License. 

def beta_enriched(prior, posterior):
    # def f(x):
    #     return beta.cdf(x, prior[0], prior[1]) * beta.pdf(x, posterior[0], posterior[1])
    # def g(x):
    #     return beta.pdf(x, posterior[0], posterior[1])
    # def h(x):
    #     return pdf_cdf_prod(x, prior, posterior)
    # # compute by integration
    # splits = integrate_splits(prior, posterior)
    # v = integrate(f, splits) / integrate(g, splits)
    
    # use closed form
    a = prior[0]
    b = prior[1]
    c = posterior[0]
    d = posterior[1]
    # See Integration.mathetmatica
    # This would be better if we computed the log of the
    # hypergeometric function, but I don't think that's generally
    # possible.
    hyper = hyper3F2aZ1(a, 1-b, a+c, a+c+d)
    scale = exp(gammaln(a) + gammaln(a+c) + gammaln(d) - gammaln(1+a) - gammaln(a+c+d) - betaln(a,b) - betaln(c,d))
    if isnan(hyper * scale):
        # This can happen if hyper and scale are 0 and inf (or vice versa).
        if prior[0] / sum(prior) > posterior[0] / sum(posterior):
            return 0.0
        return 1.0
    return clip(hyper * scale, 0, 1)

def score(prior, counts):
    ''' score a well based on the prior fit to the data and the observed counts '''
    assert prior.shape==counts.shape, "dirichletintegrate.score: array shapes do not match: "+str(prior.shape)+' and '+str(counts.shape)
    K = len(prior)
    posterior = prior + counts
    def score_idx(idx):
        prior_a = prior[idx]
        prior_b = sum(prior) - prior_a
        posterior_a = posterior[idx]
        posterior_b = sum(posterior) - posterior_a
        return beta_enriched((prior_a, prior_b), (posterior_a, posterior_b))
    return [score_idx(i) for i in range(K)]

#counts = [[178,82,102],[184,280,260]]
# counts = [[178,82,102],[204,81,118],[76,110,81],[175,117,87]]
#counts = [[178,82,102],[82,102,178],[102,178,82]]
# alpha, converged = polyafit.fit_betabinom_minka_alternating(counts)
#counts = transpose(counts)
#scores = array(score(alpha, array(counts)))
#print(scores)

# major_counts = [[57,1269,160], [16,6,9], [20,1481,34], [28,473,240], [620,609,1342], [394,129,189], [484,3154,1250], [474,849,942],
#                 [326,1336,561], [101,1338,534], [738,59,175], [597,44,101], [205,515,829], [130,358,148], [228,2476,400], 
#                 [143,3112,198], [14,284,219], [877,467,1732], [807,137,228], [470,101,173], [197,1050,893], [728,980,1054], 
#                 [230,1788,292], [602,2012,764], [541,1620,777], [447,149,265], [251,227,168], [235,2351,331], [170,235,43], 
#                 [419,1071,664], [258,1989,270], [820,362,936], [1134,130,442], [695,325,1428], [817,21,79], [432,1015,692], 
#                 [333,139,531], [89,1856,585], [483,372,1018], [1516,110,357], [91,2880,163], [949,316,472], [196,302,666], 
#                 [93,867,294], [192,473,501], [1401,82,280], [116,737,567], [135,356,180], [509,75,142], [324,926,402], 
#                 [239,969,140], [26,1222,575], [117,1084,189], [273,51,106], [201,155,736], [6,719,279], [64,528,338],
#                 [76,209,384], [20,779,148], [266,111,160], [20,898,184], [110,83,158], [1,0,0], [0,0,0], [934,373,402], 
#                 [591,85,136], [41,1235,85], [67,1370,288], [162,230,415], [3,66,5], [5,1,0], [49,1399,114], [65,752,236], 
#                 [985,127,292], [440,133,154], [304,2127,746], [444,184,710], [95,797,542], [665,27,166], [57,604,272], 
#                 [86,156,64], [643,927,970], [1469,521,777], [65,1904,223], [56,608,302], [8,582,22], [46,1956,238], [41,433,211]]

# alpha, converged = polyafit.fit_betabinom_minka_alternating(major_counts)
# print(alpha)
# counts = [57,1269,160]
# scores = array(score(alpha, array(counts)))
# print(scores)
# logScores = []
# for eachCategory in range(0,3):
#     logScores.append([log10(scores[eachCategory]) - (log10(1 - scores[eachCategory]))])
# print(logScores)

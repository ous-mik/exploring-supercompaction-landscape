from scipy.special import gamma, gammaln
from numpy import *
import pdb

def pochdivpochgen(a, b):
    'generates 100 terms at a time of pochhammer(a)/pochhammer(b)'
    count = 0
    #oldval = 0.0
    out = zeros(101, float64)
    out[0] = 1
    while True:
        out[1:] = (a + arange(count,count+100,dtype=float64)) / (b + arange(count,count+100,dtype=float64))
        out = cumprod(out)
        yield out[1:]
        out[0] = out[-1]
        count += 100

def hyper3F2Z1(a1, a2, a3, b1, b2, tol=1e-15):
    '''the hypergeometric function 3(F^~)2(a,b,z) with z=1.
    The terms are grouped under the assumption that the arguments will be:
    {a, 1-b, a+c} and {1+a, a+c+d} w/ a,b,c,d >= 0, in which case the 
    groupings result in terms all < 1 (after enough iterations)'''

    a1b1terms = pochdivpochgen(a1, b1)
    a3b2terms = pochdivpochgen(a3, b2)
    a2kterms = pochdivpochgen(a2, 1)
    
    terms = []
    while True:
        termprod = next(a1b1terms) * next(a3b2terms) * next(a2kterms)
        terms.append(termprod)
        if abs(termprod[-1]) < tol:
            break

    # sum in reverse, for better accuracy
    terms.reverse()
    #     print 'L', len(terms)
    return 1.0+sum([sum(t[::-1]) for t in terms]), terms[0][-1]

def hyp2f1mine(a, b, c):
    'hyp2f1 from scipy gives nan if any of the arguments are too large'
    return exp(gammaln(c) + gammaln(c-a-b) - gammaln(c - a) - gammaln(c - b))

def hyper3F2aZ1(a1, a2, a3, b2, tol=1e-10):
    '''same has hyper3F2Z1 but with b1 = a1+1'''
    
    if (a2 > -2):
        # force a2 no higher than -2
        # use identity for hypergeom F:
        #      (ai - 1) * F = (ai - aj - 1) * F(ai -1) + aj * F(ai-1, aj+1)
        # in this case, i = 2, j = 1
        toladjust = tol / max(1, abs((a2 - a1 - 1) / (a2 - 1)))
        temp = ((a2 - a1 -1) * hyper3F2aZ1(a1, a2 - 1, a3, b2, tol=toladjust) + a1 * hyp2f1mine(a2 - 1, a3, b2)) / (a2 - 1)
        return temp
    if (a2 < -10):
        # force a2 no lower than -10
        # use identity for hypergeom F:
        #      (ai - aj) * F = ai * F(ai + 1) - aj * F(aj + 1)
        # in this case, i = 1, j = 2
        # no adjustment for tolerance because a1 > 0, a2 < 0
        
        a2new = arange(a2, -9)
        A = [a1 * hyp2f1mine(a2n, a3, b2) for a2n in a2new]
        B = a2new
        C = (a1 - a2new)
        # compute last entry
        B[-1] = a2new[-1] * hyper3F2aZ1(a1, a2new[-1]+1, a3, b2)
        for idx in arange(len(B)-1, 0, -1):
            B[idx-1] *= (A[idx] - B[idx]) / C[idx]

        # temp = (a1 * hyp2f1mine(a2, a3, b2) - a2 * hyper3F2aZ1(a1, a2 + 1, a3, b2)) / (a1 - a2)
        return (A[0] - B[0]) / C[0]
    else:
        return hyper3F2Z1(a1, a2, a3, a1+1, b2, tol=tol)[0]

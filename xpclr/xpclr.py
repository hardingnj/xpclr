# functions to run an XPCLR scan

# TO DO:
#1 add code to read in from text file as GenotypedChunkedArray. Talk to AM.

#2 run through some kind of python analysis to id bottlenecks

#3 Think about cacheing. a Euclidean vector would be ace.

import numpy as np
from scipy.spatial.distance import squareform
from scipy.optimize import minimize_scalar
from scipy.stats import binom
from scipy.integrate import romberg
from bisect import bisect_right, bisect_left
import allel
from math import sqrt, pi
import warnings


## FUNCTIONS ###################################################################
def estimate_omega(q1, q2):
    """
    :param q1: array of alt allele freq in pop1
    :param q2: array of alt allele freq in pop2
    :return:
    """

    assert np.all(0 < q2) and np.all(q2 < 1), "No SNPs in p2 can be fixed."
    w = np.mean(((q1-q2)**2)/(q2*(1-q2)))

    return w


def determine_c(r, s, effective_pop_size=20000, min_rd=1e-7):
    """
    :param effective_pop_size: Ne of population
    :param r: recombination fraction
    :param s: selection coefficient
    :return: c which is a proxy for selection strength.
    c is [0,1] the closer to 1, the weaker the influence of selection.
    """
    if s <= 0:
        return 1.0
    else:
        return 1 - np.exp(-np.log(2*effective_pop_size)*max(r, min_rd)/s)


# This is the function that needs to be called alot and could be cythonized
# given a set of potential p1s, we wish to know the probability density at each
# value of p1.
def pdf(p1, data):

    """
    param p1: x values for which we want pdens
    :param data is an array of c, p2, var
    :return:
    """
    if isinstance(p1, float):
        p1 = np.array([p1])

    c, p2, var = data

    a_term = sqrt(2 * pi * var) ** -1

    # find all the values where p1 is greater than 1-c
    ixr = bisect_right(p1, 1 - c)

    # find all the values where p1 is less than c
    ixl = bisect_left(p1, c)

    # initialize
    r = np.zeros(p1.shape)

    # left hand side
    b_term_l = (c - p1[:ixl])/(c**2)
    c_term_l = ((p1[:ixl] - (c*p2))**2) / (2*(c**2)*var)
    r[:ixl] += a_term * b_term_l * np.exp(-c_term_l)

    # right hand side
    b_term_r = (p1[ixr:] + c - 1)/(c**2)
    c_term_r = ((p1[ixr:] + c - 1 - (c*p2))**2) / (2*(c**2)*var)
    r[ixr:] += a_term * b_term_r * np.exp(-c_term_r)

    return r


# This is a simple wrapper function to pdf that combines it with the binomial
# probability of n/x alleles.
def pdf_integral(p1, data):

    # calculate pdens for range of p1
    xj, nj, c, p2, var = data

    dens = pdf(p1, data=data[2:])
    return np.exp(np.log(dens) + binom.logpmf(xj, nj, p=p1))


# This is recoded from the implementation of xpclr. I do not think it represents
# what is discussed in the paper. Here we take the integral of eq 5 in the paper
# then take the ratio of it to eq 5 without the binomial component.
# additionally they neglect the probability of p1 being 0 or 1, I presume to
# allow the romberg integration to converge ok.
# This calculates the likelihood of a given SNP
def chen_likelihood(values):

    """
    :param values: is an array of nobserved alt alleles, total obs alleles, c,
    p2freq, and var.
    :return: The likelihood ratio of the two likelihoods.
    """

    with warnings.catch_warnings(record=True) as w:

        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")

        # These 2 function calls are major speed bottlenecks.
        i_likl = romberg(pdf_integral, a=0.001, b=0.999, args=(values,),
                         divmax=50, vec_func=True, tol=1.48e-6)

        i_base = romberg(pdf, a=0.001, b=0.999, args=(values[2:],),
                         divmax=50, vec_func=True, tol=1.48e-6)

        ratio = np.log(float(i_likl)) - np.log(float(i_base))

        if np.isnan(ratio) or ratio == -np.inf:
            ratio = -1800

        if w:
            print(w[-1].message)
            print(i_likl, i_base)

    return ratio


# similar to above but I think this is what it says in the paper
# instead of normalizing by the other distribution we take the integral over
# 0-1. Not 0.001 - 0.999 and normalize. Essentially they say af of pop1, cannot
# be 0 or 1. I think this is a convenience hack for romberg integration? NOTE
# does not work with Romberg integration!
def harding_likelihood(values):

    print("Not to be used as Romberg does not work on non cont diffs")
    with warnings.catch_warnings(record=True) as w:

        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")

        cl = np.log(romberg(pdf_integral, a=0., b=1., args=(values,),
                            divmax=50, vec_func=True, tol=1.48e-6))

        if np.isnan(cl) or cl == -np.inf:
            cl = -1800

        if w:
            print(w[-1].message)

    return cl

calculate_likelihood = chen_likelihood


def calculate_cl_romberg(sc, indata):

    """
    This function calculates the CLR for a given coeficient/position
    combination.
    :param sc: selection coefficient 0-1
    :param indata: matrix of x, n, rd, p2freq, weight, omega for N snps
    :return:
    """

    # selection coef must be bounded by 0 and 1
    # NB: XPCLR Chen et al, limit from 0 > 0.1
    if not (1 > sc >= 0):
        return np.inf

    # this will hold the marginal likelihoods for each SNP
    marginall = np.empty(indata.shape[0])

    # loop through SNPs
    for j in range(indata.shape[0]):

        xj, nj, rd, p2freq, omega, weight = indata[j]
        var = omega * (p2freq * (1 - p2freq))
        c = determine_c(rd, sc)

        # allow this function to change
        cl = calculate_likelihood(values=np.array([xj, nj, c, p2freq, var]))

        marginall[j] = weight * cl

    # take the product of the likelihoods
    ml = marginall.sum()

    # return as +ve = wrong but handy for minimize function.
    return -ml


def compute_xpclr(dat):

    # values of dist cannot be 0. Equiv to ~10 bp

    def fx(s_coef):
        return calculate_cl_romberg(s_coef, dat)

    res = minimize_scalar(fx, method="brent",
                          options={"xtol": 0.0005})

    # maximum Likelihood selec coef and recombination pos
    # removed the np.round on this line, as if rec dist up or down can go
    # beyond boundary
    maxli_sc = res.x

    # likelihood of above sc
    maximum_li = res.fun

    if maxli_sc > 0:
        null_model_li = calculate_cl_romberg(0, dat)
    else:
        null_model_li = res.fun

    if null_model_li < maximum_li:
        print("Convergence failed.", res, null_model_li)
        return np.repeat(np.nan, 3)

    # return maxL, L of null model, maxL sel coef and maxL rp
    # if null model is maxL model, sc=0, and rp is meaningless
    return -maximum_li, -null_model_li, maxli_sc


def determine_window(pos, start_, stop_, maximum_size):

    start_ix = int(np.searchsorted(a=pos, v=start_))
    stop_ix = int(np.searchsorted(a=pos, v=stop_))

    if (stop_ix - start_ix) > maximum_size:
        # window is too large. I require that the window is symmetrical
        ix = np.sort(np.random.choice(range(start_ix, stop_ix),
                                      size=maximum_size, replace=False))

    else:
        ix = np.arange(start_ix, stop_ix).astype("int")

    return ix, (stop_ix - start_ix)


def determine_weights(genotypes, ldcutoff, isphased=False):

    if isphased:
        d = genotypes.to_haplotypes()
    else:
        d = genotypes.to_n_alt()

    ld = allel.stats.ld.rogers_huff_r(np.array(d))

    above_cut = np.abs(squareform(ld)) > ldcutoff
    return 1/(1 + np.sum(above_cut, axis=1))


def xpclr_scan(gt1, gt2, bpositions, windows, geneticd=None, ldcutoff=0.95,
               phased=False, maxsnps=200, minsnps=10, rrate=1e-8,
               verbose=False):

    if geneticd is None:
        geneticd = bpositions * rrate
        if verbose:
            print("No genetic distance provided; using rrate of {0}"
                  .format(rrate))

    ac1 = gt1.count_alleles()
    ac2 = gt2.count_alleles()
    w = estimate_omega(q1=ac1.to_frequencies()[:, 1],
                       q2=ac2.to_frequencies()[:, 1])
    if verbose:
        print("omega: {0}".format(w))

    count_calls = ac1.sum(axis=1)[:]
    count_alt = ac1[:, 1]
    p2_freqs = ac2.to_frequencies()[:, 1]

    li_data = np.zeros((windows.shape[0], 3))
    nsnp = np.zeros(windows.shape[0], dtype="int")
    nsnp_avail = np.zeros(windows.shape[0], dtype="int")
    ixspan = np.zeros(windows.shape, dtype="int")

    for i, (start, end) in enumerate(windows):

        if verbose and not i % 10:
            print("Processing window {0}/{1}...".
                  format(i + 1, windows.shape[0]))

        ix, n_avail = determine_window(bpositions, start, end, maxsnps)
        nsnp[i] = ix.size
        nsnp_avail[i] = n_avail

        ixspan[i] = np.take(bpositions, (ix[0], ix[-1]))
        if nsnp[i] < minsnps:
            # if not enough data in window, skip
            li_data[i] = np.repeat(np.nan, 3)
            continue

        weights = determine_weights(gt1.take(ix, axis=0), ldcutoff=ldcutoff,
                                    isphased=phased)

        dq = np.take(geneticd, ix)
        distance = np.abs(dq - dq.mean())

        # combine_arrays into single array for easier passing
        window_data = np.vstack((count_alt.take(ix, axis=0),
                                 count_calls.take(ix, axis=0),
                                 distance, p2_freqs.take(ix),
                                 np.repeat(w, distance.size),
                                 weights)).T

        li_data[i] = compute_xpclr(window_data)

    if verbose:
        print("...done")

    # modelL, nullL, selcoef, n snps, actual window edges.
    return li_data.T[0], li_data.T[1], li_data.T[2], nsnp, nsnp_avail, ixspan
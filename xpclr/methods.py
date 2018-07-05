# functions to run an XPCLR scan

import numpy as np
from scipy.spatial.distance import squareform
from scipy.stats import binom
from scipy.integrate import romberg, quad
from bisect import bisect_right, bisect_left
import allel
from math import sqrt, pi
import warnings
from functools import lru_cache
import logging
logger = logging.getLogger(__name__)


# FUNCTIONS ###################################################################
def estimate_omega(q1, q2):
    """
    :param q1: array of alt allele freq in pop1
    :param q2: array of alt allele freq in pop2
    :return:
    """

    assert np.all(0 < q2) and np.all(q2 < 1), "No SNPs in p2 can be fixed."
    w = np.mean(((q1-q2)**2)/(q2*(1-q2)))

    return w


def determine_var(w, q2):

    """
    :param w: omega as estimated
    :param q2: allele freq of SNP in p2
    :return: sigma2, estimate of variation
    """

    return w * (q2 * (1 - q2))


def determine_c(r, s, effective_pop_size=20000, min_rd=1e-7, sf=5):
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
        c = 1 - np.exp(-np.log(2*effective_pop_size)*max(r, min_rd)/s)
        return np.round(c, sf)


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

    return dens * binom.pmf(xj, nj, p=p1)


# This is recoded from the implementation of xpclr. I do not think it represents
# what is discussed in the paper. Here we take the integral of eq 5 in the paper
# then take the ratio of it to eq 5 without the binomial component.
# additionally they neglect the probability of p1 being 0 or 1, I presume to
# allow the romberg integration to converge ok.
# This calculates the likelihood of a given SNP
@lru_cache(maxsize=2**16)
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
        # to do: http://docs.scipy.org/doc/scipy/reference/tutorial/integrate
        # .html#faster-integration-using-ctypes
        i_likl = quad(pdf_integral, a=0.001, b=0.999, args=(values,),
                      epsrel=0.001, epsabs=0, full_output=1)

        i_base = quad(pdf, a=0.001, b=0.999, args=(values[2:],),
                      epsrel=0.001, epsabs=0, full_output=1)

        if w:
            logger.warning(w[-1].message, i_likl, i_base)

    like_i, like_b = i_likl[0], i_base[0]

    if like_i == 0.0:
        ratio = -1800
    elif like_b == 0.0:
        ratio = -1800
    else:
        ratio = np.log(like_i) - np.log(like_b)

    return ratio


# similar to above but I think this is what it says in the paper
# instead of normalizing by the other distribution we take the integral over
# 0-1. Not 0.001 - 0.999 and normalize. Essentially they say af of pop1, cannot
# be 0 or 1. I think this is a convenience hack for romberg integration? NOTE
# does not work with Romberg integration!
def harding_likelihood(values):

    logger.warning("Not to be used as Romberg does not work on non cont diffs")
    with warnings.catch_warnings(record=True) as w:

        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")

        cl = np.log(romberg(pdf_integral, a=0., b=1., args=(values,),
                            divmax=50, vec_func=True, tol=1.48e-6))

        if np.isnan(cl) or cl == -np.inf:
            cl = -1800

        if w:
            logger.warning(w[-1].message)

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
        var = determine_var(w=omega, q2=p2freq)
        c = determine_c(rd, sc, sf=5)

        # allow this function to change
        cl = calculate_likelihood(values=(xj, nj, c, p2freq, var))

        marginall[j] = weight * cl

    # take the product of the likelihoods
    ml = marginall.sum()

    # return as +ve = wrong but handy for minimize function.
    return -ml


def compute_xpclr(dat, selection_cfs):

    # values of dist cannot be 0. Equiv to ~10 bp

    lliks = []
    maximum_li = np.inf
    maxli_sc = 0.0

    for s_coef in selection_cfs:

        ll = calculate_cl_romberg(s_coef, dat)
        lliks.append(ll)
        if ll < maximum_li:
            maximum_li = ll
            maxli_sc = s_coef
        else:
            break

    # maximum Likelihood selec coef and recombination pos
    # removed the np.round on this line, as if rec dist up or down can go
    # beyond boundary
    null_model_li = lliks[0]

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
        d = genotypes.to_n_alt(fill=0)

    # nans are possible, but rare, ie where only alts in A are at positions
    # missing in B. We consider these sites in LD and they are dropped.
    ld = allel.stats.ld.rogers_huff_r(d[:], fill=1.0)

    above_cut = squareform(ld**2) > ldcutoff

    # add one as self ld reported as 0
    return 1/(1 + np.sum(above_cut, axis=1))


def xpclr_scan(gt1, gt2, bpositions, windows, geneticd=None, ldcutoff=0.95,
               phased=False, maxsnps=200, minsnps=10, rrate=1e-8,
               sel_coefs=(0.0, 0.00001, 0.00005, 0.0001, 0.0002, 0.0004, 0.0006,
                          0.0008, 0.001, 0.003, 0.005, 0.01, 0.05, 0.08, 0.1,
                          0.15)):

    if geneticd is None:
        geneticd = bpositions * rrate
        logger.info("No genetic distance provided; using rrate of {0}/bp".format(rrate))

    assert minsnps >= 2, "Minimum SNPs cannot be set at any fewer than 2"

    ac1 = gt1.count_alleles()
    ac2 = gt2.count_alleles()
    w = estimate_omega(q1=ac1.to_frequencies()[:, 1],
                       q2=ac2.to_frequencies()[:, 1])
    logger.info("Omega estimated as : {0:3f}".format(w))

    count_calls = ac1.sum(axis=1)[:]
    count_alt = ac1[:, 1]
    p2_freqs = ac2.to_frequencies()[:, 1]

    li_data = np.zeros((windows.shape[0], 3))
    nsnp = np.zeros(windows.shape[0], dtype="int")
    nsnp_avail = np.zeros(windows.shape[0], dtype="int")
    ixspan = np.zeros(windows.shape, dtype="int")

    for i, (start, end) in enumerate(windows):

        if 0 == (i % 10):
            logger.debug("Processing window {0}/{1}...".format(i + 1, windows.shape[0]))

        ix, n_avail = determine_window(bpositions, start, end, maxsnps)

        nsnp[i] = ix.size
        nsnp_avail[i] = n_avail

        if nsnp[i] < minsnps:
            # if not enough data in window, skip
            li_data[i] = np.repeat(np.nan, 3)
            continue

        ixspan[i] = np.take(bpositions, (ix[0], ix[-1]))

        weights = determine_weights(gt2.take(ix, axis=0), ldcutoff=ldcutoff,
                                    isphased=phased)

        dq = np.take(geneticd, ix)
        distance = np.abs(dq - dq.mean())

        # combine_arrays into single array for easier passing
        window_data = np.vstack((count_alt.take(ix, axis=0),
                                 count_calls.take(ix, axis=0),
                                 distance, p2_freqs.take(ix),
                                 np.repeat(w, distance.size),
                                 weights)).T

        li_data[i] = compute_xpclr(window_data, sel_coefs)

    # modelL, nullL, selcoef, n snps, actual window edges.
    return li_data.T[0], li_data.T[1], li_data.T[2], nsnp, nsnp_avail, ixspan

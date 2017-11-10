# -*- coding: utf-8 -*-
"""
Complexity estimates for solving LWE.

.. moduleauthor:: Martin R. Albrecht <martinralbrecht@googlemail.com>

"""

from functools import partial
from collections import OrderedDict

from scipy.optimize import newton

from sage.modules.all import vector
from sage.calculus.var import var
from sage.functions.log import exp, log
from sage.functions.other import ceil, sqrt, floor, binomial, erf
from sage.interfaces.magma import magma
from sage.matrix.all import Matrix
from sage.misc.all import cached_function
from sage.misc.all import set_verbose, get_verbose, prod
from sage.arith.srange import srange
from sage.numerical.optimize import find_root
from sage.rings.all import QQ, RR, ZZ, RealField, PowerSeriesRing, RDF
from sage.symbolic.all import pi, e
from sage.rings.infinity import PlusInfinity

from sage.crypto.lwe import LWE, Regev, LindnerPeikert

# config

oo = PlusInfinity()
tau_default = 0.3  # τ used in uSVP
tau_prob_default = 0.1  # probability of success for given τ
cfft = 1  # convolutions mod q
enable_LP_estimates =  True  # enable LP estimates
enable_fplll_estimates = False  # enable fplll estimates


# Utility Classes #

class OutOfBoundsError(ValueError):
    """
    Used to indicate a wrong value, for example delta_0 < 1.
    """
    pass


class InsufficientSamplesError(ValueError):
    """
    Used to indicate the number of samples given is too small, especially
    useful for #samples <= 0.
    """
    pass


# Utility Functions #

def binary_search_minimum(f, start, stop, param, extract=lambda x: x, *arg, **kwds):
    """
    Return minimum of `f` if `f` is convex.

    :param start: start of range to search
    :param stop:  stop of range to search (exclusive)
    :param param: the parameter to modify when calling `f`
    :param extract: comparison is performed on `extract(f(param=?, *args, **kwds))`

    """
    return binary_search(f, start, stop, param, better=lambda x, best: extract(x)<=extract(best), *arg, **kwds)


def binary_search(f, start, stop, param, better=lambda x, best: x<=best, *arg, **kwds):
    """
    Searches for the `best` value in the interval [start,stop]
    depending on the given predicate `better`.

    :param start: start of range to search
    :param stop:  stop of range to search (exclusive)
    :param param: the parameter to modify when calling `f`
    :param better: comparison is performed by evaluating `better(current, best)`

    """
    kwds[param] = stop
    D = {}
    D[stop] = f(*arg, **kwds)
    best = D[stop]
    b = ceil((start+stop)/2)
    direction = 0
    while True:
        if b == start:
            best = D[start]
            break
        if b not in D:
            kwds[param] = b
            D[b] = f(*arg, **kwds)
        if not better(D[b], best):
            if direction == 0:
                start = b
                b = ceil((stop+b)/2)
            else:
                stop = b
                b = floor((start+b)/2)
        else:
            best = D[b]
            if b-1 not in D:
                kwds[param] = b-1
                D[b-1] = f(*arg, **kwds)
            if better(D[b-1], best):
                stop = b
                b = floor((b+start)/2)
                direction = 0
            else:
                if b+1 not in D:
                    kwds[param] = b+1
                    D[b+1] = f(*arg, **kwds)
                if not better(D[b+1], best):
                    break
                else:
                    start = b
                    b = ceil((stop+b)/2)
                    direction = 1
    return best


def cost_str(d, keyword_width=None, newline=None, round_bound=2048):
    """
    Return string of key,value pairs as a string "key0: value0, key1: value1"

    :param d:        report dictionary
    :keyword_width:  keys are printed with this width

    EXAMPLE:

    By default dicts are unordered, hence the order of the output of this function is undefined::

        sage: s = {"delta":5, "bar":2}
        sage: print cost_str(s)
        bar:         2,  delta:         5

    Use `OrderedDict` if you require ordered output::

        sage: s = OrderedDict([(u"delta", 5), ("bar",2)])
        sage: print cost_str(s)
        delta:         5,  bar:         2

    """
    if d is None:
        return
    s = []
    for k in d:
        v = d[k]
        if keyword_width:
            fmt = u"%%%ds" % keyword_width
            k = fmt % k
        if ZZ(1)/round_bound < v < round_bound or v == 0 or ZZ(-1)/round_bound > v > -round_bound:
            try:
                s.append(u"%s: %9d" % (k, ZZ(v)))
            except TypeError:
                if v < 2.0 and v >= 0.0:
                    s.append(u"%s: %9.7f" % (k, v))
                else:
                    s.append(u"%s: %9.4f" % (k, v))
        else:
            t = u"≈%s2^%.1f" % ("-" if v < 0 else "", log(abs(v), 2).n())
            s.append(u"%s: %9s" % (k, t))
    if not newline:
        return u",  ".join(s)
    else:
        return u"\n".join(s)


def cost_reorder(d, ordering):
    """
    Return a new ordered dict from the key:value pairs in `d` but reordered such that the keys in
    ordering come first.

    :param d:        input dictionary
    :param ordering: keys which should come first (in order)


    EXAMPLE::

        sage: d = OrderedDict([("a",1),("b",2),("c",3)]); d
        OrderedDict([('a', 1), ('b', 2), ('c', 3)])

        sage: cost_reorder(d, ["b","c","a"])
        OrderedDict([('b', 2), ('c', 3), ('a', 1)])

    """
    keys = list(d)
    for key in ordering:
        keys.pop(keys.index(key))
    keys = list(ordering) + keys
    r = OrderedDict()
    for key in keys:
        r[key] = d[key]
    return r


def cost_filter(d, keys):
    """
    Return new ordered dict from the key:value pairs in `d` restricted to the keys in `keys`.

    :param d:    input dictionary
    :param keys: keys which should be copied (ordered)
    """
    r = OrderedDict()
    for key in keys:
        r[key] = d[key]
    return r


def cost_repeat(d, times, repeat=None):
    """
    Return a report with all costs multiplied by `times`.

    :param d:     a cost estimate
    :param times: the number of times it should be run
    :param repeat: toggle which fields ought to be repeated and which shouldn't
    :returns:     a new cost estimate

    We maintain a local dictionary which decides if an entry is multiplied by `times` or not.
    For example, δ would not be multiplied but "\#bop" would be. This check is strict such that
    unknown entries raise an error. This is to enforce a decision on whether an entry should be
    multiplied by `times` if the function `report` reports on is called `times` often.

    EXAMPLE::

        sage: n, alpha, q = unpack_lwe(Regev(128))
        sage: print cost_str(cost_repeat(sis(n, alpha, q), 2^10))
        bkz2:   ≈2^85.4,  oracle:   ≈2^36.5,  δ_0: 1.0089681, ...
        sage: print cost_str(cost_repeat(sis(n, alpha, q), 1))
        bkz2:   ≈2^75.4,  oracle:   ≈2^26.5,  δ_0: 1.0089681, ...

    """

    do_repeat = {
        u"bop": True,
        u"rop": True,
        u"oracle": True,
        u"bkz2": True,
        u"lp": True,
        u"ds": True,
        u"fplll": True,
        u"sieve": True,
        u"enum": True,
        u"enumop": True,
        u"log(eps)": False,
        u"quantum_sieve": True,

        u"mem": False,
        u"delta_0": False,
        u"beta": False,
        u"k": False,
        u"ε": False,
        u"D_reg": False,
        u"t": False,
        u"Pr[⊥]": False,  # we are leaving probabilities alone
        u"m": False,
        u"dim": False,
        u"|v|": False,
        u"amplify": False,
        u"repeat": False,  # we deal with it below
        u"c": False,
    }

    if repeat is not None:
        for key in repeat:
            do_repeat[key] = repeat[key]

    ret = OrderedDict()
    for key in d:
        try:
            if do_repeat[key]:
                ret[key] = times * d[key]
            else:
                ret[key] = d[key]
        except KeyError:
            raise NotImplementedError(u"You found a bug, this function does not know about '%s' but should."%key)
    ret[u"repeat"] = times * ret.get("repeat", 1)
    return ret


def cost_combine(left, right, base=None):
    """Combine ``left`` and ``right``.

    :param left: cost dictionary
    :param right: cost dictionary
    :param base: add entries to ``base``

    """
    if base is None:
        cost = OrderedDict()
    else:
        cost = base
    for key in left:
        cost[key] = left[key]
    for key in right:
        cost[key] = right[key]
    return cost


def stddevf(sigma):
    """
    σ → std deviation

    :param sigma: Gaussian width parameter σ

    EXAMPLE::

        sage: n = 64.0
        sage: stddevf(n)
        25.532...
    """
    return sigma/sqrt(2*pi).n()


def sigmaf(stddev):
    """
    std deviation → σ

    :param sigma: standard deviation

    EXAMPLE::

        sage: n = 64.0
        sage: sigmaf(stddevf(n))
        64.000...
    """
    RR = stddev.parent()
    return RR(sqrt(2*pi))*stddev


def alphaf(sigma, q, sigma_is_stddev=False):
    """
    σ, q → α

    :param sigma: Gaussian width parameter (or standard deviation if `sigma_is_stddev` is `True`)
    :param q: modulus
    :param sigma_is_stddev: if `True` then `sigma` is interpreted as the standard deviation
    :returns: α = σ/q or σ·sqrt(2π)/q depending on `sigma_is_stddev`
    :rtype: real number

    """
    if sigma_is_stddev is False:
        return RR(sigma/q)
    else:
        return RR(sigmaf(sigma)/q)


def amplify(target_success_probability, success_probability, majority=False):
    """
    Return the number of trials needed to amplify current `success_probability` to
    `target_success_probability`

    :param target_success_probability: 0 < real value < 1
    :param success_probability:        0 < real value < 1
    :param majority: if `True` amplify a deicsional problem, not a computational one
       if `False` then we assume that we can check solutions, so one success suffices

    :returns: number of required trials to amplify
    """
    prec = max(53,
               2*ceil(abs(log(success_probability, 2))),
               2*ceil(abs(log(1-success_probability, 2))),
               2*ceil(abs(log(target_success_probability, 2))),
               2*ceil(abs(log(1-target_success_probability, 2))))
    RR = RealField(prec)

    if target_success_probability < success_probability:
        return RR(1)

    success_probability = RR(success_probability)
    target_success_probability = RR(target_success_probability)

    if majority:
        eps = success_probability/2
        repeat = ceil(2*log(2 - 2*target_success_probability)/log(1 - 4*eps**2))
    else:
        # target_success_probability = 1 - (1-success_probability)^trials
        repeat = ceil(log(1-target_success_probability)/log(1 -success_probability))

    return repeat


def rinse_and_repeat(f, n, alpha, q, success_probability=0.99,
                     optimisation_target=u"bkz2",
                     decision=True,
                     samples=None,
                     *args, **kwds):
    """Find best trade-off between success probability and running time.

    :param f: a function returning a cost estimate
    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  target success probability
    :param optimisation_target:  what value out to be minimized
    :param decision:             ``True`` if ``f`` solves Decision-LWE, ``False`` for Search-LWE.
    :param samples:              the number of available samples

    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)

    best = None
    step_size = 32
    i = floor(-log(success_probability, 2))
    has_solution = False
    while True:
        prob = min(2**-i, success_probability)
        try:
            current = f(n, alpha, q,
                        optimisation_target=optimisation_target,
                        success_probability=prob,
                        samples=samples,
                        *args, **kwds)
            repeat = amplify(success_probability, prob, majority=decision)
            do_repeat = None if samples is None else {"oracle": False}
            current = cost_repeat(current, repeat, do_repeat)
            has_solution = True
        except (OutOfBoundsError, InsufficientSamplesError) as err:
            key = list(best)[0] if best is not None else optimisation_target
            current = OrderedDict()
            current[key] = oo
            if get_verbose() >= 2:
                print err
        current["log(eps)"] = -i

        if get_verbose() >= 2:
            print cost_str(current)

        key = list(current)[0]
        if best is None:
            best = current
            i += step_size
            continue

        if key not in best or current[key] < best[key]:
            best = current
            i += step_size
        else:
            # we go back
            i = -best["log(eps)"] - step_size
            i += step_size/2
            if i <= 0:
                i = step_size/2
            # and half the step size
            step_size = step_size/2

        if step_size == 0:
            break

    if not has_solution:
        raise RuntimeError("No solution found for chosen parameters.")

    return best


@cached_function
def distinguish_required_m(sigma, q, success_probability, other_sigma=None):
    RR = sigma.parent()
    if other_sigma is not None:
        sigma = RR(sqrt(sigma**2 + other_sigma**2))
    adv = RR(exp(-RR(pi)*(RR(sigma/q)**2)))
    return RR(success_probability)/RR(adv**2)


def uniform_variance_from_bounds(a, b, h=None):
    """
    Variance for uniform distribution from bounds.

    :param a:
    :param b:
    :param h:        number of non-zero components.
    :returns:
    :rtype:

    """
    assert a < 0 and b > 0 and abs(a) == abs(b)
    if h is None:
        n = b - a + 1
        return (n**2 - 1)/ZZ(12)
    else:
        # sage: var("i,a,b")
        # sage: p = 1/(b-a)
        # sage: sum(i^2*p, i, a, b)
        return (2*a**3 - 2*b**3 - 3*a**2 - 3*b**2 + a - b)/(6*ZZ(a - b))


def unpack_lwe(lwe):
    """
    Return n, α, q given an LWE instance object.

    :param lwe: LWE object
    :returns: n, α, q
    :rtype: tuple

    """
    n = lwe.n
    q = lwe.K.order()
    try:
        alpha = alphaf(sigmaf(lwe.D.sigma), q)
    except AttributeError:
        # older versions of Sage use stddev, not sigma
        alpha = alphaf(sigmaf(lwe.D.stddev), q)
    return n, alpha, q


def unpack_lwe_dict(lwe):
    """
    Return dictionary consisting of n, α, q and samples given an LWE instance object.

    :param lwe: LWE object
    :returns: "n": n, "alpha": α, "q": q, "samples": samples
    :rtype: dictionary

    """
    n, alpha, q = unpack_lwe(lwe)
    samples = lwe.m
    return {"n": n, "alpha": alpha, "q": q, "samples": samples}


def preprocess_params(n, alpha, q, success_probability=None, prec=None, samples=None):
    """
    Check if parameters n, α, q are sound and return correct types.
    Also, if given, the soundness of the success probability and the
    number of samples is ensured.
    """
    if n < 1:
        raise ValueError("LWE dimension must be greater than 0.")
    if alpha <= 0:
        raise ValueError("Fraction of noise must be > 0.")
    if q < 1:
        raise ValueError("LWE modulus must be greater than 0.")
    if samples is not None and samples < 1:
        raise ValueError("Given number of samples must be greater than 0.")
    if prec is None:
        prec = 128
    RR = RealField(prec)
    n, alpha, q =  ZZ(n), RR(alpha), ZZ(q),

    samples = ZZ(samples)

    if success_probability is not None:
        if success_probability >= 1 or success_probability <= 0:
            raise ValueError("success_probability must be between 0 and 1.")
        return n, alpha, q, RR(success_probability)
    else:
        return n, alpha, q


################################
# Section 2                    #
################################


def switch_modulus(n, alpha, q, s_variance, h=None):
    """
    Return modulus switched parameters.

    :param n:        the number of variables in the LWE instance
    :param alpha:    noise size
    :param q:        modulus
    :param s_var:    the variance of the secret
    :param h:        number of non-zero components.

    If ``h`` is given, then ``s_var`` refers to the variance of non-zero components.

    EXAMPLE::

       sage: switch_modulus(128, 0.01, 65537, uniform_variance_from_bounds(0,1))
       (128, 0.0141421356237310, 410.000000000000)

       sage: switch_modulus(128, 0.001, 65537, uniform_variance_from_bounds(0,1))
       (128, 0.00141421356237310, 4094.00000000000)

       sage: switch_modulus(128, 0.001, 65537, uniform_variance_from_bounds(-5,5))
       (128, 0.00141421356237310, 25889.0000000000)

    """
    if h is not None:
        length = h
    else:
        length = n
    p = RR(ceil(sqrt(2*pi*s_variance*length/ZZ(12)) / alpha))

    if p < 32:  # some random point
        # we can't pretend everything is uniform any more, p is too small
        p = RR(ceil(sqrt(2*pi*s_variance*length*2/ZZ(12)) / alpha))
    beta = RR(sqrt(2)*alpha)
    return n, beta, p


# Lattice Reduction

def bkz_svp_repeat(n, k):
    """Return number of SVP calls in BKZ-k

    :param n: dimension
    :param k: block size

    .. note :: loosely based on experiments in [PhD:Chen13]

    """
    return 8*n


def _delta_0f(k):
    """
    Compute `δ_0` from block size `k` without enforcing `k` in ZZ.

    δ_0 for k<=40 were computed as follows:

    ```
    # -*- coding: utf-8 -*-
    from fpylll import BKZ, IntegerMatrix

    from multiprocessing import Pool
    from sage.all import mean, sqrt, exp, log, cputime

    d, trials = 320, 32

    def f((A, beta)):

        par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, flags=BKZ.AUTO_ABORT)
        q = A[-1, -1]
        d = A.nrows
        t = cputime()
        A = BKZ.reduction(A, par, float_type="dd")
        t = cputime(t)
        return t, exp(log(A[0].norm()/sqrt(q).n())/d)

    if __name__ == '__main__':
        for beta in (5, 10, 15, 20, 25, 28, 30, 35, 40):
            delta_0 = []
            t = []
            i = 0
            while i < trials:
                threads = int(open("delta_0.nthreads").read()) # make sure this file exists
                pool = Pool(threads)
                A = [(IntegerMatrix.random(d, "qary", k=d//2, bits=50), beta) for j in range(threads)]
                for (t_, delta_0_) in pool.imap_unordered(f, A):
                    t.append(t_)
                    delta_0.append(delta_0_)
                i += threads
                print u"β: %2d, δ_0: %.5f, time: %5.1fs, (%2d,%2d)"%(beta, mean(delta_0), mean(t), i, threads)
            print
    ```

    """
    small = (( 2, 1.02190),  # noqa
             ( 5, 1.01862),  # noqa
             (10, 1.01616),
             (15, 1.01485),
             (20, 1.01420),
             (25, 1.01342),
             (28, 1.01331),
             (40, 1.01295))

    if k <= 2:
        return RR(1.0219)
    elif k < 40:
        for i in range(1, len(small)):
            if small[i][0] > k:
                return RR(small[i-1][1])
    elif k == 40:
        return RR(small[-1][1])
    else:
        return RR(k/(2*pi*e) * (pi*k)**(1/k))**(1/(2*(k-1)))


def delta_0f(k):
    """
    Compute `δ_0` from block size `k`.
    """
    k = ZZ(round(k))
    return _delta_0f(k)


def k_chen_secant(delta):
    """
    Estimate required blocksize `k` for a given root-hermite factor δ based on [PhD:Chen13]_

    :param delta: root-hermite factor

    EXAMPLE::

        sage: 50 == k_chen(1.0121)
        True
        sage: 100 == k_chen(1.0093)
        True
        sage: k_chen(1.0024) # Chen reports 800
        808

    .. [PhD:Chen13] Yuanmi Chen. Réduction de réseau et sécurité concrète du chiffrement
                    complètement homomorphe. PhD thesis, Paris 7, 2013.
    """
    # newton() will produce a "warning", if two subsequent function values are
    # indistinguishable (i.e. equal in terms of machine precision). In this case
    # newton() will return the value k in the middle between the two values
    # k1,k2 for which the function values were indistinguishable.
    # Since f approaches zero for k->+Infinity, this may be the case for very
    # large inputs, like k=1e16.
    # For now, these warnings just get printed and the value k is used anyways.
    # This seems reasonable, since for such large inputs the exact value of k
    # doesn't make such a big difference.
    try:
        k = newton(lambda k: RR(_delta_0f(k) - delta), 100, fprime=None, args=(), tol=1.48e-08, maxiter=500)
        k = ceil(k)
        if k < 40:
            # newton may output k < 40. The old k_chen method wouldn't do this. For
            # consistency, call the old k_chen method, i.e. consider this try as "failed".
            raise RuntimeError("k < 40")
        return k
    except (RuntimeError, TypeError):
        # if something fails, use old k_chen method
        if get_verbose() >= 2:
            print "secant method failed, using k_chen_old(delta) instead!"
        k = k_chen_old(delta)
        return k


def k_chen_find_root(delta):
    # handle k < 40 separately
    k = ZZ(40)
    if delta_0f(k) < delta:
        return k

    try:
        k = find_root(lambda k: RR(_delta_0f(k) - delta), 40, 2**16, maxiter=500)
        k = ceil(k)
    except RuntimeError:
        # finding root failed; reasons:
        # 1. maxiter not sufficient
        # 2. no root in given interval
        k = k_chen_old(delta)
    return k


def k_chen(delta):
    # TODO: decide for one strategy (secant, find_root, old) and it's handling of errors.
    k = k_chen_find_root(delta)
    return k


def k_chen_old(delta):
    """
    Estimate required blocksize `k` for a given root-hermite factor δ based on [PhD:Chen13]_

    :param delta: root-hermite factor

    EXAMPLE::

        sage: 50 == k_chen(1.0121)
        True
        sage: 100 == k_chen(1.0093)
        True
        sage: k_chen(1.0024) # Chen reports 800
        808

    .. [PhD:Chen13] Yuanmi Chen. Réduction de réseau et sécurité concrète du chiffrement
                    complètement homomorphe. PhD thesis, Paris 7, 2013.
    """
    k = ZZ(40)

    while delta_0f(2*k) > delta:
        k *= 2
    while delta_0f(k+10) > delta:
        k += 10
    while True:
        if delta_0f(k) < delta:
            break
        k += 1

    return k


def bkz_runtime_delta_LP(delta, n):
    """
    Runtime estimation assuming the Lindner-Peikert model.
    """
    return RR(1.8/log(delta, 2) - 110 + log(2.3*10**9, 2))


def bkz_runtime_k_sieve_bdgl16_small(k, n):
    u"""
    Runtime estimation given `k` and assuming sieving is used to realise the SVP oracle.

    For small `k` we use estimates based on experiments in [BDGL16]

    :param k: block size
    :param n: lattice dimension

    ..  [BDGL16] Becker, A., Ducas, L., Gama, N., & Laarhoven, T.  (2016).  New directions in
        nearest neighbor searching with applications to lattice sieving.  In SODA 2016, (pp. 10–24).

    """
    return RR(0.387*k + 16.4 + log(bkz_svp_repeat(n, k), 2))


def bkz_runtime_k_sieve_bdgl16_asymptotic(k, n):
    u"""
    Runtime estimation given `k` and assuming sieving is used to realise the SVP oracle.

    :param k: block size
    :param n: lattice dimension

    ..  [BDGL16] Becker, A., Ducas, L., Gama, N., & Laarhoven, T.  (2016).  New directions in
        nearest neighbor searching with applications to lattice sieving.  In SODA 2016, (pp. 10–24).
    """
    # we simply pick the same additive constant 16.4 as for the experimental result in [BDGL16]
    return RR(0.292*k + 16.4 + log(bkz_svp_repeat(n, k), 2))


def bkz_runtime_k_quantum_sieve(k, n):
    """
    Runtime estimation for quantum sieving.

    ..  [LaaMosPol14] Thijs Laarhoven, Michele Mosca, & Joop van de Pol.  Finding shortest lattice
        vectors faster using quantum search.  Cryptology ePrint Archive, Report 2014/907, 2014.
        https://eprint.iacr.org/2014/907.
    """
    return RR((0.265*k + 16.4 + log(bkz_svp_repeat(n, k), 2)))


bkz_runtime_k_sieve_asymptotic  = bkz_runtime_k_sieve_bdgl16_asymptotic
bkz_runtime_k_sieve_small       = bkz_runtime_k_sieve_bdgl16_small


def bkz_runtime_k_sieve(k, n):
    u"""
     Runtime estimation given `k` and assuming sieving is used to realise the SVP oracle.

    :param k: block size
    :param n: lattice dimension

     """
    if k <= 90:
        return bkz_runtime_k_sieve_small(k, n)
    else:
        return bkz_runtime_k_sieve_asymptotic(k, n)


def bkz_runtime_k_bkz2(k, n):
    """
    Runtime estimation given `k` and assuming [CheNgu12]_ estimates are correct.

    The constants in this function were derived as follows based on Table 4 in [CheNgu12]_::

        sage: dim = [100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210, 220, 230, 240, 250]
        sage: nodes = [39.0, 44.0, 49.0, 54.0, 60.0, 66.0, 72.0, 78.0, 84.0, 96.0, 99.0, 105.0, 111.0, 120.0, 127.0, 134.0]  # noqa
        sage: times = [c + log(200,2).n() for c in nodes]
        sage: T = zip(dim, nodes)
        sage: var("a,b,c,k")
        (a, b, c, k)
        sage: f = a*k*log(k, 2.0) + b*k + c
        sage: f = f.function(k)
        sage: f.subs(find_fit(T, f, solution_dict=True))
        k |--> 0.270188776350190*k*log(k) - 1.0192050451318417*k + 16.10253135200765

    .. [CheNgu12] Yuanmi Chen and Phong Q. Nguyen. BKZ 2.0: Better lattice security estimates (Full Version).
                  2012. http://www.di.ens.fr/~ychen/research/Full_BKZ.pdf


    """
    repeat = log(bkz_svp_repeat(n, k), 2)
    return RR(0.270188776350190*k*log(k) - 1.0192050451318417*k + 16.10253135200765 + repeat)


def bkz_runtime_delta_bkz2(delta, n):
    """
    Runtime estimation extrapolated from BKZ 2.0 timings.
    """
    k = k_chen(delta)
    return bkz_runtime_k_bkz2(k, n)


def bkz_runtime_k_fplll(k, n):
    """
    Runtime estimation extrapolated from fpLLL 4.0.4 experiments
    """
    repeat = log(bkz_svp_repeat(n, k), 2)
    return RR(0.013487467331762426*k**2 - 0.28245244492771304*k + 21.017892848466957 + repeat)


def bkz_runtime_delta(delta, n, log_repeat=0):
    """
    Runtime estimates for BKZ (2.0) given δ and n
    """
    if enable_LP_estimates:
        t_lp = bkz_runtime_delta_LP(delta, n) + log_repeat

    RR = delta.parent()

    k = k_chen(delta)
    t_sieve = RR(bkz_runtime_k_sieve(k, n) + log_repeat)
    t_bkz2  = RR(bkz_runtime_k_bkz2(k, n)  + log_repeat)
    t_fplll = RR(bkz_runtime_k_fplll(k, n) + log_repeat)
    t_quantum_sieve = RR(bkz_runtime_k_quantum_sieve(k, n) + log_repeat)

    r = OrderedDict()
    r[u"delta_0"] = delta
    r[u"bkz2"] = RR(2)**t_bkz2
    r[u"beta"] = k
    if enable_LP_estimates:
        r[u"lp"] = RR(2)**t_lp
    if enable_fplll_estimates:
        r[u"fplll"] = RR(2)**t_fplll
    r[u"sieve"] = RR(2)**t_sieve
    r[u"quantum_sieve"] = RR(2)**t_quantum_sieve
    return r


def lattice_reduction_opt_m(n, q, delta):
    """
    Return the (heuristically) optimal lattice dimension `m`

    :param n:     dimension
    :param q:     modulus
    :param delta: root Hermite factor `δ_0`

    """
    return ZZ(round(sqrt(n*log(q, 2)/log(delta, 2))))


def sieve_or_enum(func):
    """
    Take minimum of sieving or enumeration for lattice-based attacks.

    :param func: a lattice-reduction based estimator
    """
    def wrapper(*args, **kwds):
        from copy import copy
        kwds = copy(kwds)
        if "optimisation_target" in kwds:
            del kwds["optimisation_target"]

        a = func(*args, optimisation_target="bkz2", **kwds)
        b = func(*args, optimisation_target="sieve", **kwds)
        if a["bkz2"] <= b["sieve"]:
            return a
        else:
            return b
    return wrapper


def mitm(n, alpha, q, success_probability=0.99, secret_bounds=None, h=None, samples=None):
    """
    Return meet-in-the-middle estimates.

    :param n: dimension
    :param alpha: noise parameter
    :param q: modulus
    :param success_probability: desired success probability
    :param secret_bounds: tuple with lower and upper bound on the secret
    :param samples: the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    ret = OrderedDict()
    RR = alpha.parent()

    m = None
    if samples is not None:
        if not samples > n:
            raise InsufficientSamplesError("Number of samples: %d" % samples)
        m = samples - n

    t = ceil(2*sqrt(log(n)))
    if secret_bounds is None:
        # assert((2*t*alpha)**m * (alpha*q)**(n/2) <= 2*n)
        m_required = ceil((log(2*n) - log(alpha*q)*(n/2))/log(2*t*alpha))
        if m is not None and m < m_required:
            raise InsufficientSamplesError("Requirement not fulfilled. Number of samples: %d - %d < %d" % (
                samples, n, m_required))
        m = m_required
        if m*(2*alpha) > 1- 1/(2*n):
            raise ValueError("Cannot find m to satisfy constraints (noise too big).")
        ret["rop"] = RR((2*alpha*q+1)**(n/2) * 2*n * m)
        ret["mem"] = RR((2*alpha*q+1)**(n/2) * m)
    else:
        a, b = secret_bounds
        # assert((2*t*alpha)**m * (b-a+1)**(n/2) <= 2*n)
        m_required = ceil(log(2*n/((b-a+1)**(n/2)))/log(2*t*alpha))
        if m is not None and m < m_required:
            raise InsufficientSamplesError("Requirement not fulfilled. Number of samples: %d - %d < %d" % (
                samples, n, m_required))
        m = m_required
        if (m*(2*alpha) > 1- 1/(2*n)):
            raise ValueError("Cannot find m to satisfy constraints (noise too big).")
        ret["rop"] = RR((b-a+1)**(n/2) * 2*n * m)
        ret["mem"] = RR((b-a+1)**(n/2) * m)

    ret["oracle"] = n + m
    ret["bop"] = RR(log(q, 2) * ret["rop"])
    return cost_reorder(ret, ["bop", "oracle", "mem"])


# BKW


def bkw(n, alpha, q, success_probability=0.99, optimisation_target="bop", prec=None, search=False, samples=None):
    """
    Estimate the cost of running BKW to solve LWE

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param optimisation_target:  field to use to decide if parameters are better
    :param prec:                 precision used for floating point computations
    :param search:               if `True` solve Search-LWE, otherwise solve Decision-LWE
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    """
    if search:
        return bkw_search(n, alpha, q, success_probability, optimisation_target, prec, samples)
    else:
        return bkw_decision(n, alpha, q, success_probability, optimisation_target, prec, samples)


def bkw_decision(n, alpha, q, success_probability=0.99, optimisation_target="bop", prec=None, samples=None):
    """
    Estimate the cost of running BKW to solve Decision-LWE following [DCC:ACFFP15]_.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param optimisation_target:  field to use to decide if parameters are better
    :param prec:                 precision used for floating point computations
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    .. [DCC:ACFFP15] Albrecht, M. R., Cid, C., Jean-Charles Faugère, Fitzpatrick, R., &
                     Perret, L. (2015). On the complexity of the BKW algorithm on LWE.
                     Designs, Codes & Cryptography, Volume 74, Issue 2, pp 325-354
    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    sigma = alpha*q

    has_samples = samples is not None
    has_enough_samples = True

    RR = alpha.parent()

    def _run(t):
        a = RR(t*log(n, 2))  # target number of adds: a = t*log_2(n)
        b = RR(n/a)  # window width
        sigma_final = RR(n**t).sqrt() * sigma  # after n^t adds we get this σ

        m = distinguish_required_m(sigma_final, q, success_probability)

        tmp = a*(a-1)/2 * (n+1) - b*a*(a-1)/4 - b/6 * RR((a-1)**3 + 3/2*(a-1)**2 + (a-1)/2)
        stage1a = RR(q**b-1)/2 * tmp
        stage1b = m * (a/2 * (n + 2))
        stage1  = stage1a + stage1b

        nrops = RR(stage1)
        nbops = RR(log(q, 2) * nrops)
        ncalls = RR(a * ceil(RR(q**b)/RR(2)) + m)
        nmem = ceil(RR(q**b)/2) * a * (n + 1 - b * (a-1)/2)

        current = OrderedDict([(u"t", t),
                               (u"bop", nbops),
                               (u"oracle", ncalls),
                               (u"m", m),
                               (u"mem", nmem),
                               (u"rop", nrops),
                               (u"a", a),
                               (u"b", b),
                               ])

        if optimisation_target != u"oracle":
            current = cost_reorder(current, (optimisation_target, u"oracle", u"t"))
        else:
            current = cost_reorder(current, (optimisation_target, u"t"))

        return current

    best_runtime = None
    best_samples = None
    best = None
    may_terminate = False
    t = RR(2*(log(q, 2) - log(sigma, 2))/log(n, 2))
    while True:
        current = _run(t)

        if get_verbose() >= 2:
            print cost_str(current)

        # Usually, both the fewest samples required and the best runtime are
        # provided when choosing 't' such that the two steps are balanced. But,
        # better safe than sorry, both cases are searched for independently.
        # So, 'best_samples' and 'best_runtime' are only used to ensure termination,
        # such that 't' is increased until both the best runtime and the fewest
        # samples were seen. The result to be returned is hold in 'best'.

        if has_samples:
            has_enough_samples = current["oracle"] < samples
            if not best_samples:
                best_samples = current
            else:
                if best_samples["oracle"] > current["oracle"]:
                    best_samples = current
                    if has_enough_samples and (
                            best is None or best[optimisation_target] > current[optimisation_target]):
                        best = current
                else:
                    if may_terminate:
                        break
                    may_terminate = True

        if not best_runtime:
            best_runtime = current
        else:
            if best_runtime[optimisation_target] > current[optimisation_target]:
                best_runtime = current
                if has_enough_samples and (best is None or best[optimisation_target] > current[optimisation_target]):
                    best = current
            else:
                if not has_samples or may_terminate:
                    break
                may_terminate = True
        t += 0.05

    if best is None:
        raise InsufficientSamplesError("Too few samples (%d given) to achieve a success probability of %f." % (
            samples, success_probability))
    return best


def bkw_search(n, alpha, q, success_probability=0.99, optimisation_target="bop", prec=None, samples=None):
    """
    Estimate the cost of running BKW to solve Search-LWE following [C:DucTraVau15]_.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param optimisation_target:  field to use to decide if parameters are better
    :param prec:                 precision used for floating point computations
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    .. [EC:DucTraVau15] Duc, A., Florian Tramèr, & Vaudenay, S. (2015). Better algorithms for
                        LWE and LWR.
    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    sigma = stddevf(alpha*q)
    eps = success_probability

    has_samples = samples is not None
    has_enough_samples = True

    RR = alpha.parent()

    # "To simplify our result, we considered operations over C to have the same
    # complexity as operations over Z_q . We also took C_FFT = 1 which is the
    # best one can hope to obtain for a FFT."
    c_cost = 1
    c_mem = 1

    def _run(t):
        a = RR(t*log(n, 2))  # target number of adds: a = t*log_2(n)
        b = RR(n/a)  # window width
        epp = (1- eps)/a

        m = lambda j, eps: 8 * b * log(q/eps) * (1 -  (2 * pi**2 * sigma**2)/(q**2))**(-2**(a-j))  # noqa

        c1 = (q**b-1)/2 * ((a-1)*(a-2)/2 * (n+1) - b/6 * (a*(a-1) * (a-2)))
        c2 = sum([m(j, epp) * (a-1-j)/2 * (n+2) for j in range(a)])
        c3 = (2*sum([m(j, epp) for j in range(a)]) + cfft * n * q**b * log(q, 2)) * c_cost
        c4 = (a-1)*(a-2) * b * (q**b - 1)/2

        nrops = RR(c1 + c2 + c3 + c4)
        nbops = RR(log(q, 2) * nrops)
        ncalls = (a-1) * (q**b - 1)/2 + m(0, eps)
        nmem = ((q**b - 1)/2 * (a-1) * (n + 1 - b*(a-2)/2)) + m(0, eps) + c_mem * q**b

        current = OrderedDict([(u"t", t),
                               (u"bop", nbops),
                               (u"oracle", ncalls),
                               (u"m", m(0, eps)),
                               (u"mem", nmem),
                               (u"rop", nrops),
                               (u"a", a),
                               (u"b", b),
                               ])

        if optimisation_target != u"oracle":
            current = cost_reorder(current, (optimisation_target, u"oracle", u"t"))
        else:
            current = cost_reorder(current, (optimisation_target, u"t"))

        return current

    best_runtime = None
    best_samples = None
    best = None
    may_terminate = False
    t = RR(2*(log(q, 2) - log(sigma, 2))/log(n, 2))
    while True:
        current = _run(t)

        if get_verbose() >= 2:
            print cost_str(current)

        # Similar to BKW-Decision, both the fewest samples required and the
        # best runtime are provided when choosing 't' such that the two steps
        # are balanced. To be safe, every 't' is tested until the fewest number
        # of samples required is found.

        if has_samples:
            has_enough_samples = current["oracle"] < samples
            if not best_samples:
                best_samples = current
            else:
                if best_samples["oracle"] > current["oracle"]:
                    best_samples = current
                    if has_enough_samples and (
                            best is None or best[optimisation_target] > current[optimisation_target]):
                        best = current
                else:
                    if may_terminate:
                        break
                    may_terminate = True

        if not best_runtime:
            best_runtime = current
        else:
            if best_runtime[optimisation_target] > current[optimisation_target]:
                best_runtime = current
                if has_enough_samples and (best is None or best[optimisation_target] > current[optimisation_target]):
                    best = current
            else:
                if not has_samples or may_terminate:
                    break
                may_terminate = True
        t += 0.05

    if best is None:
        raise InsufficientSamplesError("Too few samples (%d given) to achieve a success probability of %f." % (
            samples, success_probability))
    return best


def _bkw_coded(n, alpha, q, t2, b, success_probability=0.99, ntest=None, secret_bounds=None, h=None):
    """
    Estimate complexity of Coded-BKW as described in [C:GuoJohSta15]

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param t2:                   number of coded BKW steps (≥ 0)
    :param b:                    table size (≥ 1)
    :param success_probability:  probability of success < 1.0, IGNORED
    :param ntest:                optional parameter ntest
    :returns: a cost estimate
    :rtype: OrderedDict

    .. note::

        You probably want to call bkw_coded instead.

    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    sigma = stddevf(alpha*q)  # [C:GuoJohSta15] use σ = standard deviation
    RR = alpha.parent()

    cost = OrderedDict()

    # Our cost is mainly determined by q**b, on the other hand there are
    # expressions in q**(l+1) below, hence, we set l = b - 1. This allows to
    # achieve the performance reported in [C:GuoJohSta15].

    b = ZZ(b)
    cost["b"] = b
    l = b - 1
    cost["l"] = l

    gamma = RR(1.2)  # TODO make this dependent on success_probability
    if secret_bounds:
        d = secret_bounds[1] - secret_bounds[0] + 1
    else:
        d = 3*sigma      # TODO make this dependent on success_probability

    cost["d"] = d
    cost[u"γ"] = gamma

    def N(i, sigma_set):
        """
        Return $N_i$ for the $i$-th $[N_i, b]$ linear code.

        :param i: index
        :param sigma_set: target noise level
        """
        return floor(b/(1-log(12*sigma_set**2/ZZ(2)**i, q)/2))

    def find_ntest(n, l, t1, t2, b):
        """
        If the parameter `ntest` is not provided, we use this function to estimate it.

        :param n:  dimension > 0
        :param l:  table size for hypothesis testing
        :param t1: number of normal BKW steps
        :param t2: number of coded BKW steps
        :param b:  table size for BKW steps

        """

        # there is no hypothesis testing because we have enough normal BKW
        # tables to cover all of of n
        if t1*b >= n:
            return 0

        # solve for nest by aiming for ntop == 0
        ntest = var("nest")
        sigma_set = sqrt(q**(2*(1-l/ntest))/12)
        ncod = sum([N(i, sigma_set) for i in range(1, t2+1)])
        ntop = n - ncod - ntest - t1*b
        try:
            ntest = round(find_root(0 == ntop, 0, n))
        except RuntimeError:
            # annoyingly we get a RuntimeError when find_root can't find a
            # solution, we translate to something more meaningful
            raise ValueError("Cannot find parameters for n=%d, l=%d, t1=%d, t2=%d, b=%d"%(n, l, t1, t2, b))
        return ZZ(ntest)

    # we compute t1 from N_i by observing that any N_i ≤ b gives no advantage
    # over vanilla BKW, but the estimates for coded BKW always assume
    # quantisation noise, which is too pessimistic for N_i ≤ b.
    t1 = 0
    if ntest is None:
        ntest_ = find_ntest(n, l, t1, t2, b)
    else:
        ntest_ = ntest
    sigma_set = sqrt(q**(2*(1-l/ntest_))/12)
    Ni = [N(i, sigma_set) for i in range(1, t2+1)]
    t1 = len([e for e in Ni if e <= b])

    # there is no point in having more tables than needed to cover n
    if b*t1 > n:
        t1 = n//b
    t2 -= t1

    cost["t1"] = t1
    cost["t2"] = t2

    # compute ntest with the t1 just computed
    if ntest is None:
        ntest = find_ntest(n, l, t1, t2, b)

    # if there's no ntest then there's no `σ_{set}` and hence no ncod
    if ntest:
        sigma_set = sqrt(q**(2*(1-l/ntest))/12)
        cost[u"σ_set"] = RR(sigma_set)
        ncod = sum([N(i, sigma_set) for i in range(1, t2+1)])
    else:
        ncod = 0

    ntot = ncod + ntest
    ntop = max(n - ncod - ntest - t1*b, 0)
    cost["ncod"] = ncod    # coding step
    cost["ntop"] = ntop    # guessing step, typically zero
    cost["ntest"] = ntest  # hypothesis testing

    # Theorem 1: quantization noise + addition noise
    if secret_bounds:
        s_var = uniform_variance_from_bounds(*secret_bounds, h=h)
        coding_variance = s_var * sigma_set**2 * ntot
    else:
        coding_variance = gamma**2 * sigma**2 * sigma_set**2 * ntot
    sigma_final = RR(sqrt(2**(t1+t2) * sigma**2 + coding_variance))
    cost[u"σ_final"] = RR(sigma_final)

    # we re-use our own estimator
    M = distinguish_required_m(sigmaf(sigma_final), q, success_probability)
    cost["m"] = M
    m = (t1+t2)*(q**b-1)/2 + M
    cost["oracle"] = RR(m)

    # Equation (7)
    n_ = n - t1*b
    C0 = (m-n_) * (n+1) * ceil(n_/(b-1))
    assert(C0 >= 0)
    cost["C0(gauss)"] = RR(C0)

    # Equation (8)
    C1 = sum([(n+1-i*b)*(m - i*(q**b - 1)/2) for i in range(1, t1+1)])
    assert(C1 >= 0)
    cost["C1(bkw)"] = RR(C1)

    # Equation (9)
    C2_ = sum([4*(M + i*(q**b - 1)/2)*N(i, sigma_set) for i in range(i, t2+1)])
    C2 = RR(C2_)
    for i in range(i, t2+1):
        C2 += RR(ntop + ntest + sum([N(j, sigma_set) for j in range(1, i+1)]))*(M + (i-1)*(q**b - 1)/2)
    assert(C2 >= 0)
    cost["C2(coded)"] = RR(C2)

    # Equation (10)
    C3 = M*ntop*(2*d + 1)**ntop
    assert(C3 >= 0)
    cost["C3(guess)"] = RR(C3)

    # Equation (11)
    C4_ = 4*M*ntest
    C4 = C4_ + (2*d+1)**ntop * (cfft * q**(l+1) * (l+1) * log(q, 2) + q**(l+1))
    assert(C4 >= 0)
    cost["C4(test)"] = RR(C4)

    C = (C0 + C1 + C2 + C3+ C4)/(erf(d/sqrt(2*sigma))**ntop)  # TODO don't ignore success probability
    cost["rop"] = RR(C)
    cost["bop"] = RR(C)*log(RR(q), RR(2))
    cost["mem"] = (t1+t2)*q**b

    cost = cost_reorder(cost, ["bop", "oracle", "m", "mem", "rop", "b", "t1", "t2"])
    return cost


def bkw_coded(n, alpha, q, success_probability=0.99, secret_bounds=None, h=None,
              cost_include=("bop", "oracle", "m", "mem", "rop", "b", "t1", "t2"), samples=None):
    """
    Estimate complexity of Coded-BKW as described in [C:GuoJohSta15]
    by optimising parameters.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0, IGNORED
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    EXAMPLE::


        sage: n, alpha, q = unpack_lwe(Regev(64))
        sage: print cost_str(bkw_coded(n, alpha, q))
        bop:   ≈2^53.1,  oracle:   ≈2^39.2,  m:   ≈2^30.2,  mem:   ≈2^40.2,  rop:   ≈2^49.5,  ...

    """
    bstart = ceil(log(q, 2))

    def _run(b=2):
        # the noise is 2**(t1+t2) * something so there is no need to go beyond, say, q**2
        return binary_search(_bkw_coded, 2, min(n//b, ceil(2*log(q, 2))), "t2",
                             lambda x, best: x["rop"]<=best["rop"] and (
                                 samples is None or best["oracle"]>samples or x["oracle"]<=samples),
                             n, alpha, q, b=b, t2=0,
                             secret_bounds=secret_bounds, h=h,
                             success_probability=success_probability)

    best = binary_search(_run, 2, 3*bstart, "b", lambda x, best: x["rop"]<=best["rop"] and (
        samples is None or best["oracle"]>samples or x["oracle"]<=samples), b=2)
    # binary search cannot "fail". It just outputs some X with X["oracle"]>samples.
    if samples is not None and best["oracle"] > samples:
        raise InsufficientSamplesError("Too few samples given (%d)." % samples)
    return best


# Dual Strategy

def _sis(n, alpha, q, success_probability=0.99, optimisation_target=u"bkz2", secret_bounds=None, h=None, samples=None):
    """Estimate cost of solving LWE by solving LWE.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param secret_bounds:        ignored
    :param h:                    ignored
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict

    """

    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    f = lambda eps: RR(sqrt(log(1/eps)/pi))  # noqa
    RR = alpha.parent()

    # we are solving Decision-LWE
    log_delta_0 = log(f(success_probability)/alpha, 2)**2 / (4*n*log(q, 2))
    delta_0 = RR(2**log_delta_0)
    m_optimal = lattice_reduction_opt_m(n, q, delta_0)
    if samples is None or samples > m_optimal:
        m = m_optimal
    else:
        if not samples > 0:
            raise InsufficientSamplesError("Number of samples: %d" % samples)
        m = samples
        log_delta_0 = log(f(success_probability)/alpha, 2)/m - RR(log(q, 2)*n)/(m**2)
        delta_0 = RR(2**log_delta_0)

    # check for valid delta
    if delta_0 < 1:
        raise OutOfBoundsError(u"Detected delta_0 = %f < 1. Too few samples?!" % delta_0)

    ret = bkz_runtime_delta(delta_0, m)
    ret[u"oracle"] = m
    ret[u"|v|"] = RR(delta_0**m * q**(n/m))
    ret[u"dim"] = m
    if optimisation_target != u"oracle":
        ret = cost_reorder(ret, [optimisation_target, u"oracle"])
    else:
        ret = cost_reorder(ret, [optimisation_target])
    return ret


sis = partial(rinse_and_repeat, _sis)


# Decoding

@cached_function
def gsa_basis(n, q, delta, m):
    """
    Create the basis lengths.

    :param n: determinant is q^n
    :param q:  determinant is q^n
    :param delta: root-Hermite factor
    :param m: lattice dimension

    .. note:: based on the GSA in [RSA:LinPei11]_

    .. [RSA:LinPei11] Richard Lindner and Chris Peikert. Better key sizes (and attacks) for LWE-based encryption.
                      In Aggelos Kiayias, editor, CT-RSA 2011, volume 6558 of LNCS, pages 319–339. Springer,
                      February 2011.
    """
    log_delta = RDF(log(delta))
    log_q = RDF(log(q))
    qnm = log_q*(n/m)
    qnm_p_log_delta_m = qnm + log_delta*m
    tmm1 = RDF(2*m/(m-1))
    b = [(qnm_p_log_delta_m - log_delta*(tmm1 * i)) for i in xrange(m)]
    b = [log_q - b[-1-i] for i in xrange(m)]
    b = map(lambda x: x.exp(), b)
    return b


def enum_cost(n, alpha, q, eps, delta_0, m=None, B=None, step=1, enums_per_clock=-15.1):
    """
    Estimates the runtime for performing enumeration.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param eps:
    :param delta_0:
    :param m:
    :param B:
    :param step:                 changes the increments for the values of d[i]
    :param enums_per_clock:      the log of the number of enumerations computed per clock cycle
    :returns: a cost estimate
    :rtype: OrderedDict
    """

    RR = alpha.parent()
    step = RDF(step)

    if B is None:
        if m is None:
            m = lattice_reduction_opt_m(n, q, delta_0)
        B = gsa_basis(n, q, delta_0, m)

    d = [RDF(1)]*m
    bd = [d[i] * B[i] for i in xrange(m)]
    scaling_factor = RDF(sqrt(pi) / (2*alpha*q))
    probs_bd = [RDF((bd[i]  * scaling_factor)).erf() for i in xrange(m)]
    success_probability = prod(probs_bd)

    if RR(success_probability).is_NaN():
        # try in higher precision
        step = RR(step)
        d = [RR(1)]*m
        bd = [d[i] * B[i] for i in xrange(m)]
        scaling_factor = RR(sqrt(pi) / (2*alpha*q))
        probs_bd = [RR((bd[i]  * scaling_factor)).erf() for i in xrange(m)]
        success_probability = prod(probs_bd)

        if success_probability.is_NaN():
            return OrderedDict([(u"delta_0", delta_0),
                                ("enum", oo),
                                ("enumop", oo)])

    # if m too small, probs_bd entries are of magnitude 1e-10 or
    # something like that. Therefore, success_probability=prod(probs_bd)
    # results in success_probability==0 and so, the loop never terminates.
    # To prevent this, success_probability should be calculated when needed,
    # i.e. at the end of each iteration and not be used to calculate things.
    # Then, a new problem arises: for step=1 this can take very long time,
    # starting at, for example, success_probability==1e-300.
    # Since this is a (rare) special case, an error is thrown.
    if not success_probability > 0:
        raise InsufficientSamplesError("success_probability == 0! Too few samples?!")

    bd = map(list, zip(bd, range(len(bd))))
    bd = sorted(bd)

    import bisect

    last_success_probability = success_probability

    while success_probability < RDF(eps):
        v, i = bd.pop(0)
        d[i] += step
        v += B[i]*step
        last_success_probability = success_probability
        success_probability /= probs_bd[i]
        probs_bd[i] = (v * scaling_factor).erf()
        success_probability *= probs_bd[i]
        bisect.insort_left(bd, [v, i])

        if success_probability == 0 or last_success_probability >= success_probability:
            return OrderedDict([(u"delta_0", delta_0),
                                ("enum", oo),
                                ("enumop", oo)])

    r = OrderedDict([(u"delta_0", delta_0),
                     ("enum", RR(prod(d))),
                     ("enumop", RR(prod(d)) / RR(2)**enums_per_clock)])

    return r


def _decode(n, alpha, q, success_probability=0.99,
            enums_per_clock=-15.1, optimisation_target="bkz2",
            secret_bounds=None, h=None, samples=None):
    """
    Estimates the optimal parameters for decoding attack

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param enums_per_clock:      the log of the number of enumerations computed per clock cycle
    :param optimisation_target:  lattice reduction estimate to use
    :param secret_bounds:        ignored
    :param h:                    ignored
    :param samples:              the number of available samples
    :returns: a cost estimate
    :rtype: OrderedDict
    """

    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    if samples is not None and not samples > 1:
        raise InsufficientSamplesError("Number of samples: %d" % samples)

    RR = alpha.parent()

    # Number of samples are'nt considered here, because only a "good" starting delta_0
    # is needed and there is no "good" value for number of samples known,
    # i.e. the given number of samples may not produce "good" results
    delta_0m1 = _sis(n, alpha, q, success_probability=success_probability)[u"delta_0"] - 1

    step = RR(1.05)
    direction = -1

    def combine(enum, bkz):
        current = OrderedDict()
        current["rop"]  = enum["enumop"] + bkz[optimisation_target]

        for key in bkz:
            current[key] = bkz[key]
        for key in enum:
            current[key] = enum[key]
        current[u"oracle"]  = m
        current = cost_reorder(current, ["rop", "oracle", optimisation_target])
        return current

    depth = 6
    current = None
    while True:
        delta_0 = 1 + delta_0m1

        if delta_0 >= 1.0219 and current is not None:  # LLL is enough
            break

        m_optimal = lattice_reduction_opt_m(n, q, delta_0)
        if samples is None or samples > m_optimal:
            m = m_optimal
        else:
            m = samples
        bkz = bkz_runtime_delta(delta_0, m)
        bkz["dim"] = m

        enum = enum_cost(n, alpha, q, success_probability, delta_0, m, enums_per_clock=enums_per_clock)
        current = combine(enum, bkz)

        # if lattice reduction is cheaper than enumration, make it more expensive
        if current[optimisation_target] < current["enumop"]:
            prev_direction = direction
            direction = -1
            if direction != prev_direction:
                step = 1 + RR(step-1)/2
            delta_0m1 /= step
        elif current[optimisation_target] > current["enumop"]:
            prev_direction = direction
            direction = 1
            delta_0m1 *= step
        else:
            break
        if direction != prev_direction:
            step = 1 + RR(step-1)/2
            depth -= 1
        if depth == 0:
            break

    return current


decode = partial(rinse_and_repeat, _decode, decision=False)


# uSVP

def kannan(n, alpha, q, tau=tau_default, tau_prob=tau_prob_default, success_probability=0.99,
           optimisation_target="bkz2", samples=None):
    """
    Estimate optimal parameters for using Kannan-embedding to solve CVP.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param success_probability:  probability of success < 1.0
    :param samples:              the number of available samples

    :returns: a cost estimate
    :rtype: OrderedDict
    """
    # TODO: Include the attack described in
    # [RSA:BaiGal14] Shi Bai and Steven D. Galbraith. An improved compression technique
    #                for signatures based on learning with errors. In Josh Benaloh,
    #                editor, Topics in Cryptology – CT-RSA 2014, volume 8366 of Lecture
    #                Notes in Computer Science, pages 28–47, San Francisco, CA, USA,
    #                February 25–28, 2014. Springer, Heidelberg, Germany.
    # The estimation of computational cost is the same as kannan with dimension samples=n+m.
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    RR = alpha.parent()

    log_delta_0 = log(tau*alpha*sqrt(e), 2)**2/(4*n*log(q, 2))
    delta_0 = RR(2**log_delta_0)

    m_optimal = lattice_reduction_opt_m(n, q, delta_0)
    if samples is None or samples > m_optimal:
        m = m_optimal
    else:
        if not samples > 0:
            raise InsufficientSamplesError("Number of samples: %d" % samples)
        m = samples
        delta_0 = RR((q**(1-n/m)*sqrt(1/(e)) / (tau*alpha*q))**(1.0/m))

    # check for valid delta
    if delta_0 < 1:
        raise OutOfBoundsError(u"Detected delta_0 = %f < 1. Too few samples?!" % delta_0)

    l2 = q**(1-n/m) * sqrt(m/(2*pi*e))
    if l2 > q:
        raise NotImplementedError(u"Case where λ_2 = q not implemented.")

    repeat = amplify(success_probability, tau_prob)

    r = bkz_runtime_delta(delta_0, m, log(repeat, 2.0))
    r[u"oracle"] = repeat*m if samples is None else m
    r[u"m"] = m
    r = cost_reorder(r, [optimisation_target, "oracle"])
    if get_verbose() >= 2:
        print cost_str(r)
    return r


# Gröbner bases

def gb_complexity(m, n, d, omega=2, call_magma=True, d2=None):
    """Estimate the complexity of computing a Gröbner basis.

    Estimation is done for `m` polynomials of degree `d` in `n`
    variables under the assumption that the system is semi-regular.

    If `d2` is not ``None`` then `n` polynomials of degree are added to
    the system. This is to encode restrictions on the solution. For
    example, if the solution is either `0` or `1`, then `(x_i)⋅(x_i+1)`
    would evaluate to zero on it for any `x_i`.

    :param m: number of polynomials (integer > 0)
    :param n: number of variables (integer > 0)
    :param d: degree of all input polynomials
    :param omega: linear algebra exponent, i.e. matrix-multiplication costs `O(n^ω)` operations.
    :param call_magma: use Magma to perform computation (highly recommended)
    :param d2: secondary degree (integer > 0 or ``None``)

    :returns: A dictionary containing the expected degree of regularity "Dreg"
    (assuming the system is semi-regular), the number of base ring operations "rop",
    and the memory requirements in base ring elements "mem".

    :rtype: OrderedDict
    """
    if m > n**d:
        m = n**d

    if call_magma:
        R = magma.PowerSeriesRing(QQ, 2*n)
        z = R.gen(1)
        coeff = lambda f, d: f.Coefficient(d)  # noqa
    else:
        R = PowerSeriesRing(QQ, "z", 2*n)
        z = R.gen()
        coeff = lambda f, d: f[d]  # noqa

    if d2 is None:
        s = (1-z**d)**m / (1-z)**n
    else:
        s = (1-z**d)**m * (1-z**d2)**n / (1-z)**n

    retval = OrderedDict([("Dreg", None), ("rop", None)])

    for dreg in xrange(2*n):
        if coeff(s, dreg) < 0:
            break
    else:
        return retval
    retval["Dreg"] = dreg
    retval["rop"] = RR(binomial(n + dreg, dreg)**omega)
    retval["mem"] = RR(binomial(n + dreg, dreg)**2)
    return retval


def arora_gb(n, alpha, q, success_probability=0.99, omega=2, call_magma=True, guess=0, d2=None, samples=None):

    if samples is not None:
        from warnings import warn
        warn("Given number of samples is ignored for arora_gb()!")

    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability,
                                                         prec=2*log(n, 2)*n)

    RR = alpha.parent()
    stddev = RR(stddevf(RR(alpha)*q))

    if stddev >= 1.1*sqrt(n):
        return None

    if d2 is True:
        d2 = 2*ceil(3*stddev)+1

    ps_single = lambda C: RR(1 - (2/(C*RR(sqrt(2*pi))) * exp(-C**2/2)))  # noqa

    m = floor(exp(RR(0.75)*n))
    d = n
    t = ZZ(floor((d-1)/2))
    C = t/stddev
    pred = gb_complexity(m, n-guess, d, omega, call_magma, d2=d2)
    pred["t"] = t
    pred["oracle"] = m
    pred[u"Pr[⊥]"] = RR(m*(1-ps_single(C)))
    pred["bop"] = log(q, 2) * pred["rop"]
    pred = cost_reorder(pred, ["t", "bop", "oracle", "Dreg"])

    if get_verbose() >= 2:
        print "PREDICTION:"
        print cost_str(pred)
        print
        print "ESTIMATION:"

    t = ceil(t/3)
    best = None
    stuck = 0
    for t in srange(t, n):
        d = 2*t + 1
        C = RR(t/stddev)
        if C < 1:  # if C is too small, we ignore it
            continue
        # Pr[success]^m = Pr[overall success]
        m = log(success_probability, 2) / log(ps_single(C), 2)
        if m < n:
            continue
        m = floor(m)

        current = gb_complexity(m, n-guess, d, omega, call_magma, d2=d2)

        if current["Dreg"] is None:
            continue

        current["t"] = t
        current[u"Pr[⊥]"] = RR(1-success_probability)
        current["rop"] *= RR((3*stddev)**guess)
        current["bop"] = log(q, 2) * current["rop"]
        current["oracle"] = m

        current = cost_reorder(current, ["bop", "oracle", "t", "Dreg"])

        if get_verbose() >= 2:
            print cost_str(current)

        if best is None:
            best = current
        else:
            if best["rop"] > current["rop"]:
                best = current
                stuck = 0
            else:
                stuck += 1
                if stuck >= 5:
                    break
    return best


# Exhaustive Search for Small Secrets

def small_secret_guess(f, n, alpha, q, secret_bounds, h=None, samples=None, **kwds):
    size = secret_bounds[1]-secret_bounds[0] + 1
    best = None
    step_size = 16
    fail_attempts, max_fail_attempts = 0, 5
    while step_size >= n:
        step_size /= 2
    i = 0
    while True:
        if i<0:
            break

        try:
            try:
                # some implementations make use of the secret_bounds parameter
                current = f(n-i, alpha, q, secret_bounds=secret_bounds, samples=samples, **kwds)
            except TypeError:
                current = f(n-i, alpha, q, samples=samples, **kwds)
        except (OutOfBoundsError, RuntimeError, InsufficientSamplesError) as err:
            if get_verbose() >= 2:
                print type(err).__name__, ":", err
            i += abs(step_size)
            fail_attempts += 1
            if fail_attempts > max_fail_attempts:
                break
            continue
        if h is None or i<h:
            repeat = size**i
        else:
            # TODO: this is too pessimistic
            repeat = (size)**h * binomial(i, h)
        do_repeat = None if samples is None else {"oracle": False}
        current = cost_repeat(current, repeat, do_repeat)

        key = list(current)[0]
        if best is None:
            best = current
            i += step_size
        else:
            if best[key] > current[key]:
                best = current
                i += step_size
            else:
                step_size = -1*step_size//2
                i += step_size

        if step_size == 0:
            break
    if best is None:
        raise RuntimeError("No solution could be found.")
    return best


# Modulus Switching

def decode_small_secret_mod_switch_and_guess(n, alpha, q, secret_bounds, h=None, samples=None, **kwds):
    """Solve LWE by solving BDD for small secret instances.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param secret_bounds:
    :param h:                    number of non-zero components in the secret
    :param samples:              the number of available samples

    """
    s_var = uniform_variance_from_bounds(*secret_bounds, h=h)
    n, alpha, q = switch_modulus(n, alpha, q, s_var, h=h)
    return small_secret_guess(decode, n, alpha, q, secret_bounds, h=h, samples=samples, **kwds)


def kannan_small_secret_mod_switch_and_guess(n, alpha, q, secret_bounds, h=None, samples=None, **kwds):
    """Solve LWE by Kannan embedding for small secret instances.

    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param secret_bounds:
    :param h:                    number of non-zero components in the secret.
    :param samples:              the number of available samples

    """
    s_var = uniform_variance_from_bounds(*secret_bounds, h=h)
    n, alpha, q = switch_modulus(n, alpha, q, s_var, h=h)
    return small_secret_guess(kannan, n, alpha, q, secret_bounds, h=h, samples=samples, **kwds)


# Bai's and Galbraith's uSVP Attack

def _bai_gal_small_secret(n, alpha, q, secret_bounds, tau=tau_default, tau_prob=tau_prob_default,
                          success_probability=0.99,
                          optimisation_target="bkz2",
                          h=None, samples=None):
    """
    :param n:                    dimension > 0
    :param alpha:                fraction of the noise α < 1.0
    :param q:                    modulus > 0
    :param tau:                  0 < τ ≤ 1.0
    :param success_probability:  probability of success < 1.0
    :param optimisation_target:  field to use to decide if parameters are better
    :param h:                    number of non-zero components in the secret

    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    RR = alpha.parent()

    stddev = stddevf(alpha*q)
    a, b = secret_bounds

    if h is None:
        xi = RR(2)/(b-a)
    else:
        assert -a == b  # TODO: don't be so lazy
        # we compute the stddev of |s| and xi to scale each component to σ on average
        variance = sum([ZZ(h)/n * i**2 for i in range(a, b+1) if i])
        s_stddev = variance.sqrt()
        xi = ZZ(1)/s_stddev

    num = (log(q/stddev) - log(tau*sqrt(4*pi*e)))**2 * log(q/stddev)
    den = n*(2*log(q/stddev)-log(xi))**2

    log_delta_0 = RR(num/den)

    delta_0 = RR(e**log_delta_0)
    m_prime_optimal = ceil(sqrt(n*(log(q)-log(stddev))/log_delta_0))
    if samples is None or samples > m_prime_optimal - n:
        m_prime = m_prime_optimal
    else:
        m = samples
        m_prime = m+n
        num = m_prime*(log(q/stddev) - log(2*tau*sqrt(pi*e))) +n*log(xi)-n*log(q/stddev)
        den = m_prime**2
        log_delta_0 = RR(num/den)
        delta_0 = RR(e**log_delta_0)

    m = m_prime - n

    l2 = RR((q**m * (xi*stddev)**n)**(1/m_prime) * sqrt(m_prime/(2*pi*e)))
    if l2 > q:
        raise NotImplementedError("Case λ_2 = q not implemented.")

    repeat = amplify(success_probability, tau_prob)

    r = bkz_runtime_delta(delta_0, m_prime, log(repeat, 2))
    r[u"oracle"] = repeat*m if samples is None else m
    r[u"m"] = m

    if optimisation_target != u"oracle":
        r = cost_reorder(r, [optimisation_target, u"oracle"])
    else:
        r = cost_reorder(r, [optimisation_target])

    if get_verbose() >= 2:
        print cost_str(r)
    return r


def bai_gal_small_secret(n, alpha, q, secret_bounds, tau=tau_default, tau_prob=tau_prob_default,
                         success_probability=0.99,
                         optimisation_target="bkz2",
                         h=None, samples=None):
    """
    Bai's and Galbraith's uSVP attack + small secret guessing [ACISP:BaiGal14]_

    :param n: dimension > 0
    :param alpha: fraction of the noise α < 1.0
    :param q: modulus > 0
    :param tau: 0 < τ ≤ 1.0
    :param success_probability: probability of success < 1.0
    :param optimisation_target: field to use to decide if parameters are better
    :param h: number of non-zero components in the secret
    :param samples: the number of available samples

    .. [ACISP:BaiGal14] Bai, S., & Galbraith, S. D. (2014). Lattice decoding attacks on binary
       LWE. In W. Susilo, & Y. Mu, ACISP 14 (pp.  322–337).

    """
    return small_secret_guess(_bai_gal_small_secret, n, alpha, q, secret_bounds,
                              tau=tau, tau_prob=tau_prob,
                              success_probability=0.99,
                              optimisation_target=optimisation_target,
                              h=h, samples=samples)


# Small, Sparse Secret SIS

def success_probability_drop(n, h, k, fail=0):
    """
    Probability ``k`` randomly sampled components have at most
    ``fail`` non-zero components amongst them.

    :param n: dimension of LWE samples
    :param h: number of non-zero components
    :param k: number of components to ignore
    :param fail: we tolerate ``fail`` number of non-zero components
        amongst the ``k`` ignored components
    """

    N = n         # population size
    K = n-h       # number of success states in the population
    n = k         # number of draws
    k = n - fail  # number of observed successes
    return (binomial(K, k)*binomial(N-K, n-k)) / binomial(N, n)


def drop_and_solve(f, n, alpha, q, secret_bounds=None, h=None,
                   success_probability=0.99,
                   optimisation_target=u"bkz2", postprocess=False, **kwds):
    """
    Solve instances of dimension ``n-k`` with increasing ``k`` using
    ``f`` and pick parameters such that cost is reduced.

    :param n: dimension
    :param alpha: noise parameter
    :param q: modulus q
    :param secret_bounds: lower and upper bound on the secret
    :param h: number of non-zero components of the secret
    :param success_probability: target success probability
    :param optimisation_target: what value out to be minimized
    """
    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)

    RR = alpha.parent()

    best = None

    # too small a step size leads to an early abort, too large a step
    # size means stepping over target
    step_size = int(n/32)

    k = 0
    while True:
        current = f(n-k, alpha, q,
                    success_probability=max(1-1/RR(2)**80, success_probability),
                    optimisation_target=optimisation_target,
                    h=h, secret_bounds=secret_bounds, **kwds)

        cost_lat  = current[optimisation_target]
        cost_post = 0
        probability = success_probability_drop(n, h, k)
        if postprocess:
            repeat = current["repeat"]
            dim    = current["dim"]
            for i in range(1, k):
                # compute inner products with rest of A
                cost_post_i = 2 * repeat * dim * k
                # there are (k)(i) positions and max(s_i)-min(s_i) options per position
                # for each position we need to add/subtract the right elements
                cost_post_i += repeat * binomial(k, i) * (secret_bounds[1]-secret_bounds[0])**i  * i
                if cost_post + cost_post_i >= cost_lat:
                    postprocess = i
                    break
                cost_post += cost_post_i
                probability += success_probability_drop(n, h, k, i)

        current["rop"] = cost_lat + cost_post
        current = cost_repeat(current, 1/probability)
        current["k"] = k
        current["postprocess"] = postprocess
        current = cost_reorder(current, ["rop"])

        key = list(current)[0]
        if best is None:
            best = current
            k += step_size
            continue

        if current[key] < best[key]:
            best = current
            k += step_size
        else:
            # we go back
            k = best["k"] - step_size
            k += step_size/2
            if k <= 0:
                k = step_size/2
            # and half the step size
            step_size = step_size/2

        if step_size == 0:
            break

    return best


sis_drop_and_solve = partial(drop_and_solve, sis)
decode_drop_and_solve = partial(drop_and_solve, decode)
bai_gal_drop_and_solve = partial(drop_and_solve, bai_gal_small_secret)


def sis_small_secret_mod_switch(n, alpha, q, secret_bounds, h=None,
                                success_probability=0.99,
                                optimisation_target=u"bkz2",
                                c=None,
                                use_lll=False,
                                samples=None):
    """
    Estimate cost of solveing LWE by finding small `(y,x/c)` such that
    `y A = c x`.

    :param n:                   dimension
    :param alpha:               noise parameter
    :param q:                   modulus q
    :param secret_bounds:       lower and upper bound on the secret
    :param h:                   number of non-zero components of the secret
    :param success_probability: target success probability
    :param optimisation_target: what value out to be minimized
    :param c:                   explicit constant `c`
    :param use_lll:             use LLL calls to produce more small vectors

    """

    if samples is not None:
        from warnings import warn
        warn("Given number of samples is ignored for sis_small_secret_mod_switch()!")

    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability)
    RR = alpha.parent()

    assert(secret_bounds[0] == -1 and secret_bounds[1] == 1)

    if h is None:
        B = ZZ(secret_bounds[1] - secret_bounds[0] + 1)
        h = ceil((B-1)/B * n)

    # stddev of the error
    e = stddevf(alpha*q).n()

    delta_0 = sis(n, alpha, q, optimisation_target=optimisation_target)["delta_0"]

    best = None

    if c is None:
        c = RR(e*sqrt(2*n - n)/sqrt(h))

    if use_lll:
        scale = 2
    else:
        scale = 1

    while True:

        m = lattice_reduction_opt_m(n, c*q, delta_0)

        # the vector found will have norm
        v = scale * delta_0**m * (q/c)**(n/m)

        # each component has stddev v_
        v_ = v/RR(sqrt(m))

        # we split our vector in two parts.
        # 1. v_r is multiplied with the error e (dimension m-n)
        # 2. v_l is the rounding noise (dimension n)

        # scale last q components down again.
        v_r = e*sqrt(m-n)*v_
        v_l = c*sqrt(h)*v_

        repeat = max(distinguish_required_m(v_r, q, success_probability, other_sigma=v_l), RR(1))

        ret = bkz_runtime_delta(delta_0, m, log(repeat, 2))

        if use_lll:
            if "lp" in ret:
                keys = ("bkz2", "sieve", "lp")
            else:
                keys = ("bkz2", "sieve")
            ret_lll = bkz_runtime_delta(delta_0, m)
            for key in keys:
                # CN11: LLL: n^3 log^2 B
                ret_lll[key] += (n**3 * log(q, 2)**2) * repeat
            ret = ret_lll

        ret[u"oracle"] = m * repeat
        ret[u"repeat"] = repeat
        ret[u"dim"] = m
        ret[u"c"] = c
        ret = cost_reorder(ret, [optimisation_target, u"oracle"])

        if get_verbose() >= 2:
            print cost_str(ret)

        if best is None:
            best = ret

        if ret[optimisation_target] > best[optimisation_target]:
            break

        best = ret
        delta_0 = delta_0 + RR(0.00005)

    return best


def applebaum_transform(n, alpha, q, m, secret_bounds, h=None):
    """
    Swap part of the error vector with the secret as suggested in [EPRINT:GenHalSma12]_

    ..  [EPRINT:GenHalSma12] Gentry, C., Halevi, S., & Smart, N.  P.  (2012).  Homomorphic
        evaluation of the AES circuit.

    :param n:                    dimension
    :param alpha:                noise parameter
    :param q:                    modulus q
    :param m:                    number of samples in lattice.
    :param secret_bounds:        lower and upper bound on the secret
    :param h:                    number of non-zero components of the secret
    :param success_probability:  target success probability
    """
    # we update alpha to reflect that we're replacing part of the error by the secret
    s_var = uniform_variance_from_bounds(*secret_bounds, h=h)
    e_var = stddevf(alpha*q)**2

    stddev_ = ((n*s_var + (m-n)*e_var)/m).sqrt()
    alpha_ = alphaf(stddev_, q, sigma_is_stddev=True)
    alpha = alpha_
    return n, alpha, q


def applebaum(f, n, alpha, q, m, secret_bounds, h=None,
              success_probability=0.99,
              optimisation_target=u"bkz2", **kwds):
    """Run ``f`` after transforming LWE instance using ``applebaum_transform``.

    :param f:                    LWE solving cost function
    :param n:                    dimension
    :param alpha:                noise parameter
    :param q:                    modulus q
    :param m:                    number of samples in lattice.
    :param secret_bounds:        lower and upper bound on the secret
    :param h:                    number of non-zero components of the secret
    :param success_probability:  target success probability
    :param optimisation_target:  use this field to optimise

    """
    n, alpha, q = applebaum_transform(n, alpha, q, m, secret_bounds, h)
    return f(n, alpha, q, success_probability=success_probability, optimisation_target=optimisation_target)


sis_applebaum = partial(applebaum, sis)
decode_applebaum =  partial(applebaum, decode)


# BKW for Small Secrets


def bkw_small_secret_variances(q, a, b, kappa, o, RR=None):
    """
    Helper function for small secret BKW variant.

    :param q:
    :param a:
    :param b:
    :param kappa:
    :param o:
    :param RR:
    :returns:
    :rtype:

    """
    if RR is None:
        RR = RealField()
    q = RR(q)
    a = RR(a).round()
    b = RR(b)
    n = a*b
    kappa = RR(kappa)
    T = RR(2)**(b*kappa)
    n = RR(o)/RR(T*(a+1)) + RR(1)

    U_Var = lambda x: (x**2 - 1)/12  # noqa
    red_var   = 2*U_Var(q/(2**kappa))

    if o:
        c_ = map(RR, [0.0000000000000000,
                      0.4057993538687922, 0.6924478992819291, 0.7898852691349439,
                      0.8441959360364506, 0.8549679124679972, 0.8954469872316165,
                      0.9157093365103325, 0.9567635780119543, 0.9434245442818547,
                      0.9987153221343770])

        M = Matrix(RR, a, a)  # rows are tables, columns are entries those tables
        for l in range(M.ncols()):
            for c in range(l, M.ncols()):
                M[l, c] = U_Var(q)

        for l in range(1, a):
            for i in range(l):
                M[l, i] = red_var + sum(M[i+1:l].column(i))

                bl = b*l
                if round(bl) < len(c_):
                    c_tau = c_[round(bl)]
                else:
                    c_tau = RR(1)/RR(5)*RR(sqrt(bl)) + RR(1)/RR(3)

                f = (c_tau*n**(~bl) + 1 - c_tau)**2
                for i in range(l):
                    M[l, i] = M[l, i]/f

        v = vector(RR, a)
        for i in range(a):
            v[i] = red_var + sum(M[i+1:].column(i))
    else:
        v = vector(RR, a)
        for i in range(a)[::-1]:
            v[i] = 2**(a-i-1) * red_var

    if get_verbose() >= 3:
        print log(red_var, 2).str(), [RealField(14)(log(x, 2)).str() for x in v]

    return v


def bkw_small_secret(n, alpha, q, success_probability=0.99, secret_bounds=(0, 1), t=None, o=0, samples=None):  # noqa
    """
    :param n:               number of variables in the LWE instance
    :param alpha:           standard deviation of the LWE instance
    :param q:               size of the finite field (default: n^2)
    :param secret_bounds:   minimum and maximum value of secret
    :param samples:         the number of available samples
    """

    def sigma2f(kappa):
        v = bkw_small_secret_variances(q, a, b, kappa, o, RR=RR)
        return sigmaf(sum([b * e * secret_variance for e in v], RR(0)).sqrt())

    def Tf(kappa):
        return min(q**b, ZZ(2)**(b*kappa))/2

    def ops_tf(kappa):
        T = Tf(kappa)
        return T * (a*(a-1)/2 * (n+1) - b*a*(a-1)/4 - b/6 * ((a-1)**3 + 3/2*(a-1)**2 + 1/RR(2)*(a-1)))

    def bkwssf(kappa):
        ret = OrderedDict()
        ret[u"κ"] = kappa
        m = distinguish_required_m(sigma_final, q, success_probability, sigma2f(kappa))
        ret["m"] = m
        ropsm = (m + o)  * (a/2 * (n + 2))
        ropst = ops_tf(kappa)
        ret["rop"] = ropst + ropsm
        ret["bop"] = log(q, 2) * ret["rop"]
        T = Tf(kappa)
        ret["mem"] = T * a * (n + 1 - b * (a-1)/2)
        ret["oracle"] = T * a + ret["m"] + o
        return ret

    n, alpha, q, success_probability = preprocess_params(n, alpha, q, success_probability, prec=4*n)
    RR = alpha.parent()
    sigma = alpha*q

    has_samples = samples is not None

    if o is None:
        best = bkw_small_secret(n, alpha, q, success_probability, secret_bounds, t=t, o=0)
        o = best["oracle"]/2
        while True:
            try:
                current = bkw_small_secret(n, alpha, q, success_probability, secret_bounds, t=t, o=o)
            except InsufficientSamplesError:
                break
            if best is None or (current["bop"] < best["bop"] and (not has_samples or current["oracle"] <= samples)):
                best = current
            if current["bop"] > best["bop"]:
                break
            if get_verbose() >= 2:
                print cost_str(current)

            o = o/2
        if has_samples and best["oracle"] > samples:
            raise InsufficientSamplesError("No solution could be found with given samples (%d)" % samples)
        return best

    if t is None:
        t = RR(2*(log(q, 2) - log(sigma, 2))/log(n, 2))
        best = None
        while True:
            try:
                current = bkw_small_secret(n, alpha, q, success_probability, secret_bounds, t=t, o=o)
            except InsufficientSamplesError:
                break
            if best is None or (current["bop"] < best["bop"] and (not has_samples or current["oracle"] <= samples)):
                best = current
            if current["bop"] > best["bop"]:
                break
            if get_verbose() >= 2:
                print cost_str(current)
            t += 0.01
        if has_samples and best["oracle"] > samples:
            raise InsufficientSamplesError("No solution could be found with given samples (%d)" % samples)
        return best

    secret_variance = uniform_variance_from_bounds(*secret_bounds)
    secret_variance = RR(secret_variance)

    a = RR(t*log(n, 2))  # the target number of additions: a = t*log_2(n)
    b = n/a  # window width b = n/a
    sigma_final = RR(n**t).sqrt() * sigma  # after n^t additions we get this stddev
    transformation_noise = sqrt(n * 1/RR(12) * secret_variance)
    kappa = ceil(log(round(q*transformation_noise/stddevf(sigma)), 2.0)) + 1

    if kappa > ceil(log(q, 2)):
        kappa = ceil(log(q, 2))

    best = None
    while kappa > 0:
        current = bkwssf(kappa)
        if best is None or (current["bop"] < best["bop"] and (not has_samples or current["oracle"] <= samples)):
            best = current
        if current["bop"] > best["bop"]:
            break
        kappa -= 1
    if has_samples and best["oracle"] > samples:
        raise InsufficientSamplesError("No solution could be found with given samples (%d)" % samples)

    best["o"] = o
    best["t"] = t
    best["a"] = a
    best["b"] = b
    best = cost_reorder(best, ["bop", "oracle", "t", "m", "mem"])
    return best


# Arora-GB for Small Secrets

def arora_gb_small_secret(n, alpha, q, secret_bounds, h=None, samples=None, **kwds):
    """FIXME! briefly describe function

    :param n:
    :param alpha:
    :param q:
    :param secret_bounds:
    :param h:
    :returns:
    :rtype:

    """
    a, b = secret_bounds
    s_var = uniform_variance_from_bounds(*secret_bounds, h=h)
    n, alpha, q = switch_modulus(n, alpha, q, s_var, h=h)
    return arora_gb(n, alpha, q, d2=b-a+1, **kwds)


# Toplevel function

def estimate_lwe(n, alpha, q, samples=None, skip=None, small=False, secret_bounds=None, h=None):
    """
    Estimate the complexity of solving LWE with the given parameters.

    :param n:
    :param alpha:
    :param q:
    :param skip:
    :param small:
    :param secret_bounds:
    :returns:
    :rtype:

    EXAMPLE::

        sage: n, alpha, q = unpack_lwe(Regev(64))
        sage: set_verbose(1)
        sage: d = estimate_lwe(n,alpha,q, skip="arora-gb")
          mitm   bop:  ≈2^133.2,  oracle:        52,    mem:  ≈2^129.6,    rop:  ≈2^129.6
           bkw   bop:   ≈2^53.1,  oracle:   ≈2^39.2,      m:   ≈2^30.2,    mem:   ≈2^40.2, ...
           sis  bkz2:   ≈2^34.5,  oracle:   ≈2^11.7,    δ_0: 1.0130466,      k:        40, ...
           dec   bop:   ≈2^33.9,  oracle:      1288,    δ_0: 1.0157998,   bkz2:   ≈2^32.8, ...
        kannan  bkz2:   ≈2^35.3,  oracle:   ≈2^12.9,    δ_0: 1.0171085,      k:        40, ...

    """
    if not small:
        algorithms = OrderedDict([("mitm", mitm),
                                  ("bkw", bkw_coded),
                                  ("sis", sis),
                                  ("dec", decode),
                                  ("kannan", kannan),
                                  ("arora-gb", arora_gb)])
    else:
        algorithms = OrderedDict([("mitm", mitm),
                                  ("bkw", bkw_coded),
                                  ("sis", partial(drop_and_solve, sis_small_secret_mod_switch)),
                                  ("dec", decode_small_secret_mod_switch_and_guess),
                                  ("kannan", kannan_small_secret_mod_switch_and_guess),
                                  ("baigal", bai_gal_small_secret),
                                  ("arora-gb", arora_gb_small_secret)])

    if skip is None:
        skip = []
    try:
        skip = [s.strip().lower() for s in skip.split(",")]
    except AttributeError:
        pass
    skip = [s.strip().lower() for s in skip]

    alg_width = max(len(key) for key in set(algorithms).difference(skip))
    cost_kwds = {"keyword_width": 5}

    results = OrderedDict()
    for alg in algorithms:
        if alg not in skip:
            algf = algorithms[alg]
            if alg in ("dec", "sis", "kannan", "baigal"):
                algf = sieve_or_enum(algf)
            try:
                if small:
                    tmp = algf(n, alpha, q, secret_bounds=secret_bounds, h=h, samples=samples)
                else:
                    tmp = algf(n, alpha, q, samples=samples)
                if tmp:
                    results[alg] = tmp
                    if get_verbose() >= 1:
                        print ("%%%ds" % alg_width) % alg,
                        print cost_str(results[alg], **cost_kwds)
            except Exception as e:
                print u"Algorithm '%s' failed with message '%s'"%(alg, e)

    return results


# Plots

def plot_costs(LWE, N, skip=None, filename=None, small=False, secret_bounds=None):
    plots = {}
    for n in N:
        lwe = LWE(n)
        r = estimate_lwe(skip=skip, small=small, secret_bounds=secret_bounds, **unpack_lwe_dict(lwe))
        if get_verbose() >= 1:
            print

        for key in r:
            value = r[key].values()[0]
            plots[key] = plots.get(key, tuple()) + ((n, log(value, 2)),)

    colors = ("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD")

    import matplotlib.pyplot as plt
    plt.clf()
    plt.figure(1)

    for i, plot in enumerate(plots):
        x, y = [x_ for x_, y_ in plots[plot]], [y_ for x_, y_ in plots[plot]]
        plt.plot(x, y, label=plot, color=colors[i], linewidth=1.5)

    plt.legend(loc=2)
    plt.xlabel("n")
    plt.ylabel("$\log_2$(bop)")
    if small:
        plt.title(u"%s (%d-%d), $s ← %s^n$"%(LWE.__name__, N[0], N[-1], secret_bounds))
    else:
        plt.title(u"%s (%d-%d)"%(LWE.__name__, N[0], N[-1]))
    if filename is None:
        if small:
            small_str = "-(%d,%d)"%(secret_bounds[0], secret_bounds[1])
        else:
            small_str = ""
        filename="%s%s-%d-%d.pdf"%(LWE.__name__, small_str, N[0], N[-1])
    plt.savefig(filename, dpi=128)


def plot_fhe_costs(L, N, skip=None, filename=None, small=False, secret_bounds=None):
    plots = {}
    for n in N:
        params = fhe_params(L, n)
        r = estimate_lwe(*params, skip=skip, small=small, secret_bounds=secret_bounds)
        if get_verbose() >= 1:
            print

        for key in r:
            value = r[key].values()[0]
            plots[key] = plots.get(key, tuple()) + ((n, log(value, 2)),)

    colors = ("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD")

    import matplotlib.pyplot as plt
    plt.clf()
    plt.figure(1)

    for i, plot in enumerate(plots):
        x, y = [x_ for x_, y_ in plots[plot]], [y_ for x_, y_ in plots[plot]]
        plt.plot(x, y, label=plot, color=colors[i], linewidth=1.5)

    plt.legend(loc=2)
    plt.xlabel("n")
    plt.ylabel("$\log_2$(bop)")
    if small:
        plt.title(u"FHE (%d-%d), $L=%d$, $s ← %s^n$"%(N[0], N[-1], L, secret_bounds))
    else:
        plt.title(u"FHE (%d-%d), $L=%d$"%(N[0], N[-1], L))
    if filename is None:
        if small:
            small_str = "-(%d,%d)"%(secret_bounds[0], secret_bounds[1])
        else:
            small_str = ""
        filename="FHE-%d%s-%d-%d.pdf"%(L, small_str, N[0], N[-1])
    plt.savefig(filename, dpi=128)


# LaTeX tables

dfs = "%4.0f"

latex_config = {
    "mitm":     OrderedDict([("bop", dfs), ("mem", dfs), ("oracle", dfs)]),
    "bkw":      OrderedDict([("bop", dfs), ("mem", dfs), ("oracle", dfs)]),
    "arora-gb": OrderedDict([("bop", dfs), ("mem", dfs), ("oracle", dfs)]),
    "sis":      OrderedDict([("bkz2", dfs), ("sieve", dfs), ("oracle", dfs), ("repeat", dfs)]),
    "kannan":   OrderedDict([("bkz2", dfs), ("sieve", dfs), ("oracle", dfs), ("repeat", dfs)]),
    "baigal":   OrderedDict([("bkz2", dfs), ("sieve", dfs), ("oracle", dfs), ("repeat", dfs)]),
    "dec":      OrderedDict([("rop", dfs), ("enum", dfs), ("oracle", dfs), ("repeat", dfs)]),
}


def latex_cost_header(cur):
    header = []
    header.append(r"\begin{scriptsize}")

    pretty_algorithm_names = {
        "mitm": "MitM",
        "bkw":  "Coded-BKW",
        "arora-gb": "Arora-GB",
        "sis":  "SIS",
        "kannan": "Kannan",
        "baigal": "Bai-Gal",
        "dec": "Dec"
    }

    pretty_column_names = {
        "oracle": "$\\Ldis$",
        "repeat": "g",
    }

    line = [r"\begin{tabular}{r"]
    for alg in cur:
        line.append("@{\hskip 8pt}")
        line.append("r" * len([key for key in latex_config[alg].keys() if key in cur[alg]]))
    line.append("}")

    header.append("".join(line))
    header.append(r"\toprule")

    line = ["    "]
    for alg in cur:
        count = len([key for key in latex_config[alg].keys() if key in cur[alg]])
        line.append(r"\multicolumn{%d}{c}{%s}"%(count, pretty_algorithm_names[alg]))
    line = " & ".join(line) + "\\\\"
    header.append(line)
    header.append(r"\midrule")

    line = [" $n$"]

    for alg in cur:
        line.extend([pretty_column_names.get(key, key) for key in latex_config[alg].keys() if key in cur[alg]])

    line = " & ".join(line) + "\\\\"
    header.append(line)
    header.append(r"\midrule")
    return header


def latex_cost_row(cur):
    line = []
    for alg in cur:
        cost = cur[alg]
        for col, format in latex_config[alg].iteritems():
            if (col == "repeat" and col in cost) or col != "repeat":
                line.append(format % log(cost[col], 2))
    return line


def latex_cost_footer(name):
    footer = []
    footer.append(r"\bottomrule")
    footer.append(r"\end{tabular}")
    footer.append(r"\end{scriptsize}")
    footer.append(r"\caption{%s}" % name)
    return footer


def latex_costs(LWE, N, skip=None, small=False, secret_bounds=None):

    ret = []
    for i, n in enumerate(N):
        line = ["%4d"%n]
        lwe = LWE(n)
        cur = estimate_lwe(skip=skip, small=small, secret_bounds=secret_bounds, **unpack_lwe_dict(lwe))
        line.extend(latex_cost_row(cur))
        line = " & ".join(line) + "\\\\"
        ret.append(line)
        if get_verbose() >= 1:
            print

    header = latex_cost_header(cur)
    if small:
        name = "%s with $\s[(i)] \sample \{%d,%d\}$"%(LWE.__name__, secret_bounds[0], secret_bounds[1])
    else:
        name = LWE.__name__
    footer = latex_cost_footer(name)

    ret = header + ret + footer

    return "\n".join(ret)


def fhe_params(L, n):
    # Homomorphic Evaluation of the AES Circuit talks about σ^2 as variance so σ is stddev not width
    # parameter
    stddev = RR(3.2)
    xi = ZZ(8)
    q = ZZ(2)**(16.5*L + 5.4) * xi**(2*L-3) * n**L
    alpha = sigmaf(stddev)/q
    return n, alpha, q


def latex_fhe_costs(N, l, secret_bounds, skip=None):
    ret = []
    for n in N:
        line = ["%6d"%n]
        params = fhe_params(l, n)
        cur = estimate_lwe(*params, skip=skip, small=True, secret_bounds=secret_bounds)
        line.extend(latex_cost_row(cur))
        line = " & ".join(line) + "\\\\"
        ret.append(line)
        if get_verbose() >= 1:
            print

    header = latex_cost_header(cur)

    name = "FHE with $L=%d$ with $\s[(i)] \sample \{%d,%d\}$"%(l, secret_bounds[0], secret_bounds[1])
    footer = latex_cost_footer(name)

    ret = header + ret + footer
    return "\n".join(ret)


def make_all_tables():
    N = (64, 128, 256, 512, 1024)
    print latex_costs(Regev, N, skip=["arora-gb"])
    print
    print latex_costs(Regev, N, small=True, secret_bounds=(0, 1), skip=["arora-gb"])
    print
    print latex_costs(LindnerPeikert, N)
    print
    print latex_costs(LindnerPeikert, N, small=True, secret_bounds=(0, 1), skip=["arora-gb"])

    print latex_fhe_costs([2**i for i in range(6, 12)], l=2,  skip="Arora-GB", secret_bounds=(0, 1))
    print
    print latex_fhe_costs([2**i for i in range(6, 15)], l=10, skip="Arora-GB", secret_bounds=(0, 1))
    print


def make_all_plots():
    v = get_verbose()
    set_verbose(1)
    N = range(64, 400, 16)
    plot_costs(Regev, N, skip=["arora-gb", "mitm"])
    plot_costs(LindnerPeikert, N, skip=["arora-gb", "mitm"])

    plot_costs(Regev, N, small=True, secret_bounds=(0, 1), skip=["arora-gb"])
    plot_costs(LindnerPeikert, N, small=True, secret_bounds=(0, 1), skip=["arora-gb"])
    set_verbose(v)


class SimpleLWE(LWE):
    """
    LWE parameters with `σ=\sqrt{n/2π}` and `q = n^2`.
    """
    def __init__(self, n):
        """
        LWE parameters with `σ=\sqrt{n/2π}` and `q = n^2`.

        :param n: security parameter n (= dimension)

        """

        from sage.stats.distributions.discrete_gaussian_integer import DiscreteGaussianDistributionIntegerSampler
        from sage.rings.arith import next_prime

        q = ZZ(next_prime(n**2))
        s = sqrt(n)
        D = DiscreteGaussianDistributionIntegerSampler(s/sqrt(2*pi.n()), q)
        LWE.__init__(self, n=n, q=q, D=D, secret_dist=(0, 1))

# -----------------------------------------------------------------------------
# File:    python/zinull/rng.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All simulated data.
# Purpose: The random numbers of R (>= 3.6.0, default kinds): set.seed() with
#          its initial scrambling, the Mersenne-Twister, unif_rand(), and the
#          samplers built on it, ported operation by operation: runif(),
#          rnorm() by inversion, exponential draws (exp_rand), rbinom()
#          (inversion and BTPE), rpois() (Ahrens and Dieter PD) and rgeom().
#          A seed gives the same draws as in R, so the same simulated data.
#          One stream per process, as in R.
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull import rng
#          rng.set_seed(889192); rng.runif(3); rng.rbinom(1, 250, 0.2)
# Source:  R 4.6.1: src/main/RNG.c, src/nmath/snorm.c, qnorm.c, sexp.c,
#          rbinom.c, rpois.c, rgeom.c, runif.c; src/library/stats/src/random.c
#          (recycling of the parameters, missing values).
# Notes:   R's samplers keep the constants of the last parameters between
#          calls; these are functions of the parameters alone, so they are
#          computed at every call here. muladd() marks the products that R's
#          C code adds with one rounding on arm64 (see rmath.py).
# -----------------------------------------------------------------------------
import math

import numpy as np
from numba import njit

from .rmath import DBL_MAX, M_1_SQRT_2PI, M_LN2, _pow_di, muladd

N = 624  # words of the Mersenne-Twister state
M = 397

# The stream: N state words, then the position of the next word.
_state = np.zeros(N + 1, dtype=np.uint32)
_state[N] = N + 1



@njit(cache=True)
def _seed(state, seed):
    """Fill the state from a seed as `set.seed()` does.

    Parameters: state (uint32 array of length N + 1), seed (integer).
    Returns: None; `state` is overwritten.
    """
    s = np.uint32(seed)
    for _ in range(50):  # initial scrambling
        s = np.uint32(69069) * s + np.uint32(1)
    for j in range(N + 1):
        s = np.uint32(69069) * s + np.uint32(1)
        if j > 0:
            state[j - 1] = s
    state[N] = N  # the first word is replaced by the position


@njit(cache=True)
def _unif(state):
    """One uniform number in (0, 1): `unif_rand()`.

    Parameters: state (the stream). Returns: float.
    """
    mti = state[N]
    if mti >= N:
        for kk in range(N):
            y = (state[kk] & np.uint32(0x80000000)) | (state[(kk + 1) % N] & np.uint32(0x7FFFFFFF))
            v = state[(kk + M) % N] ^ (y >> np.uint32(1))
            if y & np.uint32(1):
                v ^= np.uint32(0x9908B0DF)
            state[kk] = v
        mti = 0
    y = state[mti]
    state[N] = mti + 1
    y ^= y >> np.uint32(11)
    y ^= (y << np.uint32(7)) & np.uint32(0x9D2C5680)
    y ^= (y << np.uint32(15)) & np.uint32(0xEFC60000)
    y ^= y >> np.uint32(18)
    x = float(y) * 2.3283064365386963e-10
    if x <= 0.0:
        return 0.5 * 2.328306437080797e-10
    if 1.0 - x <= 0.0:
        return 1.0 - 0.5 * 2.328306437080797e-10
    return x


@njit(cache=True)
def _horner(r, coefficients):
    """Polynomial in r, highest power first, evaluated as qnorm() does.

    Parameters: r, coefficients (array).
    Returns: float.
    """
    value = coefficients[0]
    for c in coefficients[1:]:
        value = muladd(value, r, c)
    return value


_CENTRE_NUM = np.array([2509.0809287301226727, 33430.575583588128105, 67265.770927008700853,
                        45921.953931549871457, 13731.693765509461125, 1971.5909503065514427,
                        133.14166789178437745, 3.387132872796366608])
_CENTRE_DEN = np.array([5226.495278852854561, 28729.085735721942674, 39307.89580009271061,
                        21213.794301586595867, 5394.1960214247511077, 687.1870074920579083,
                        42.313330701600911252, 1.0])
_MIDDLE_NUM = np.array([7.7454501427834140764e-4, .0227238449892691845833, .24178072517745061177,
                        1.27045825245236838258, 3.64784832476320460504, 5.7694972214606914055,
                        4.6303378461565452959, 1.42343711074968357734])
_MIDDLE_DEN = np.array([1.05075007164441684324e-9, 5.475938084995344946e-4,
                        .0151986665636164571966, .14810397642748007459, .68976733498510000455,
                        1.6763848301838038494, 2.05319162663775882187, 1.0])
_TAIL_NUM = np.array([2.01033439929228813265e-7, 2.71155556874348757815e-5,
                      .0012426609473880784386, .026532189526576123093, .29656057182850489123,
                      1.7848265399172913358, 5.4637849111641143699, 6.6579046435011037772])
_TAIL_DEN = np.array([2.04426310338993978564e-15, 1.4215117583164458887e-7,
                      1.8463183175100546818e-5, 7.868691311456132591e-4,
                      .0148753612908506148525, .13692988092273580531, .59983220655588793769,
                      1.0])


@njit(cache=True)
def _qnorm(p):
    """Standard normal quantile of a probability in (0, 1): `qnorm()`.

    Parameters: p. Returns: float.
    """
    q = p - 0.5
    if abs(q) <= 0.425:
        r = muladd(-q, q, 0.180625)
        return q * _horner(r, _CENTRE_NUM) / _horner(r, _CENTRE_DEN)
    r = np.sqrt(-np.log(1.0 - p if q > 0 else p))
    if r <= 5.0:
        r += -1.6
        value = _horner(r, _MIDDLE_NUM) / _horner(r, _MIDDLE_DEN)
    else:
        r += -5.0
        value = _horner(r, _TAIL_NUM) / _horner(r, _TAIL_DEN)
    return -value if q < 0.0 else value


# ----- Normal and exponential draws -----

@njit(cache=True)
def _norm_rand(state):
    """One standard normal number by inversion: `norm_rand()`.

    Two uniform numbers per value. Parameters: state. Returns: float.
    """
    big = 134217728.0  # 2^27
    u = _unif(state)
    u = int(big * u) + _unif(state)
    return _qnorm(u / big)


_EXP_Q = np.array([0.6931471805599453, 0.9333736875190459, 0.9888777961838675,
                   0.9984959252914960, 0.9998292811061389, 0.9999833164100727,
                   0.9999985691438767, 0.9999998906925558, 0.9999999924734159,
                   0.9999999995283275, 0.9999999999728814, 0.9999999999985598,
                   0.9999999999999289, 0.9999999999999968, 0.9999999999999999,
                   1.0000000000000000])


@njit(cache=True)
def _exp_rand(state):
    """One standard exponential number: `exp_rand()` (Ahrens and Dieter 1972).

    Parameters: state. Returns: float.
    """
    q = _EXP_Q
    a = 0.0
    u = _unif(state)
    while u <= 0.0 or u >= 1.0:
        u = _unif(state)
    while True:
        u += u
        if u > 1.0:
            break
        a += q[0]
    u -= 1.0
    if u <= q[0]:
        return a + u
    i = 0
    ustar = _unif(state)
    umin = ustar
    while True:
        ustar = _unif(state)
        if umin > ustar:
            umin = ustar
        i += 1
        if not u > q[i]:
            break
    return muladd(umin, q[0], a)


# ----- Binomial -----

@njit(cache=True, error_model="numpy")
def _stirling_tail(x2, x1):
    return (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0


@njit(cache=True, error_model="numpy")
def _rbinom(state, nin, pp):
    """One binomial number: `rbinom()` of the C library (inversion when
    n p < 30, else BTPE of Kachitvichyanukul and Schmeiser 1988).

    Parameters: state, nin (number of trials), pp (probability).
    Returns: float, NaN for a parameter that is missing or out of range.
    """
    if math.isnan(nin) or math.isinf(nin):
        return np.nan
    r = np.rint(nin)
    if r != nin:
        return np.nan
    if math.isnan(pp) or math.isinf(pp) or r < 0 or pp < 0.0 or pp > 1.0:
        return np.nan
    if r == 0 or pp == 0.0:
        return 0.0
    if pp == 1.0:
        return r
    if r >= 2147483647:
        return np.nan  # R draws by qbinom() here; not ported
    n = int(r)

    p = min(pp, 1.0 - pp)
    q = 1.0 - p
    np_ = n * p
    r = p / q
    g = r * (n + 1)
    ix = 0

    if np_ < 30.0:  # inversion
        qn = _pow_di(q, n)
        while True:
            ix = 0
            f = qn
            u = _unif(state)
            found = False
            while True:
                if u < f:
                    found = True
                    break
                if ix > 110:
                    break
                u -= f
                ix += 1
                f *= (g / ix - r)
            if found:
                break
    else:  # BTPE
        ffm = np_ + p
        m = int(ffm)
        fm = float(m)
        npq = np_ * q
        p1 = int(muladd(2.195, math.sqrt(npq), -(4.6 * q))) + 0.5
        xm = fm + 0.5
        xl = xm - p1
        xr = xm + p1
        c = 0.134 + 20.5 / (15.3 + fm)
        al = (ffm - xl) / muladd(-xl, p, ffm)
        xll = al * muladd(0.5, al, 1.0)
        al = (xr - ffm) / (xr * q)
        xlr = al * muladd(0.5, al, 1.0)
        p2 = p1 * (1.0 + c + c)
        p3 = p2 + c / xll
        p4 = p3 + c / xlr
        while True:
            u = _unif(state) * p4
            v = _unif(state)
            if u <= p1:  # triangular region
                ix = int(muladd(-p1, v, xm) + u)
                break
            if u <= p2:  # parallelogram region
                x = xl + (u - p1) / c
                v = muladd(v, c, 1.0) - abs(xm - x) / p1
                if v > 1.0 or v <= 0.0:
                    continue
                ix = int(x)
            elif u > p3:  # right tail
                ix = int(xr - math.log(v) / xlr)
                if ix > n:
                    continue
                v = v * (u - p3) * xlr
            else:  # left tail
                ix = int(xl + math.log(v) / xll)
                if ix < 0:
                    continue
                v = v * (u - p2) * xll
            k = abs(ix - m)
            if k <= 20 or k >= npq / 2 - 1:  # explicit evaluation
                f = 1.0
                if m < ix:
                    for i in range(m + 1, ix + 1):
                        f *= (g / i - r)
                elif m != ix:
                    for i in range(ix + 1, m + 1):
                        f /= (g / i - r)
                if v <= f:
                    break
            else:  # squeezing with bounds on log(f(x))
                amaxp = (k / npq) * (muladd(float(k), k / 3.0 + 0.625, 0.1666666666666) / npq + 0.5)
                ynorm = -k * k / (2.0 * npq)
                alv = math.log(v)
                if alv < ynorm - amaxp:
                    break
                if alv <= ynorm + amaxp:
                    x1 = float(ix + 1)
                    f1 = fm + 1.0
                    z = n + 1 - fm
                    w = n - ix + 1.0
                    z2 = z * z
                    x2 = x1 * x1
                    f2 = f1 * f1
                    w2 = w * w
                    bound = muladd(xm, math.log(f1 / x1), (n - m + 0.5) * math.log(z / w))
                    bound = muladd(float(ix - m), math.log(w * p / (x1 * q)), bound)
                    bound = bound + _stirling_tail(f2, f1) + _stirling_tail(z2, z) + \
                        _stirling_tail(x2, x1) + _stirling_tail(w2, w)
                    if alv <= bound:
                        break
    if pp > 0.5:
        ix = n - ix
    return float(ix)


# ----- Poisson and geometric -----

_FACT = np.array([1.0, 1.0, 2.0, 6.0, 24.0, 120.0, 720.0, 5040.0, 40320.0, 362880.0])


@njit(cache=True, error_model="numpy")
def _rpois(state, mu):
    """One Poisson number: `rpois()` of the C library (Ahrens and Dieter
    1982: table look-up below a mean of 10, normal deviates from 10 on).

    Parameters: state, mu (mean). Returns: float, NaN for a mean that is
    missing, infinite or negative.
    """
    if math.isnan(mu) or math.isinf(mu) or mu < 0:
        return np.nan
    if mu <= 0.0:
        return 0.0

    if mu < 10.0:  # inversion over a table of cumulative probabilities
        m = max(1, int(mu))
        table_end = 0  # the table pp[] is filled up to this index
        pp = np.empty(36)
        p0 = math.exp(-mu)
        p = p0
        q = p0
        while True:
            u = _unif(state)
            if u <= p0:
                return 0.0
            if table_end != 0:
                start = 1 if u <= 0.458 else min(table_end, m)
                for k in range(start, table_end + 1):
                    if u <= pp[k]:
                        return float(k)
                if table_end == 35:
                    continue
            table_end += 1
            for k in range(table_end, 36):
                p *= mu / k
                q += p
                pp[k] = q
                if u <= q:
                    return float(k)
            table_end = 35

    s = math.sqrt(mu)
    d = 6.0 * mu * mu
    big_l = math.floor(mu - 1.1484)
    g = muladd(s, _norm_rand(state), mu)
    pois = -1.0
    u = 0.0
    fk = 0.0
    difmuk = 0.0
    if g >= 0.0:
        pois = math.floor(g)
        if pois >= big_l:
            return pois
        fk = pois
        difmuk = mu - fk
        u = _unif(state)
        if d * u >= difmuk * difmuk * difmuk:
            return pois

    omega = M_1_SQRT_2PI / s
    b1 = 0.0416666666666666667 / mu
    b2 = 0.3 * b1 * b1
    c3 = 0.1428571428571428571 * b1 * b2
    c2 = muladd(-15.0, c3, b2)
    c1 = muladd(45.0, c3, muladd(-6.0, b2, b1))
    c0 = muladd(-15.0, c3, muladd(3.0, b2, 1.0 - b1))
    c = 0.1069 / mu

    quotient_step = g >= 0.0  # first pass with the normal deviate (step Q)
    E = 0.0
    while True:
        if not quotient_step:
            E = _exp_rand(state)
            u = muladd(2.0, _unif(state), -1.0)
            t = 1.8 + (abs(E) if u >= 0 else -abs(E))
            if not t > -0.6744:
                continue
            pois = math.floor(muladd(s, t, mu))
            fk = pois
            difmuk = mu - fk
        # 'subroutine' F: px, py, fx, fy
        if pois < 10:
            px = -mu
            py = math.pow(mu, pois) / _FACT[int(pois)]
        else:
            delta = 0.0833333333333333333 / fk
            delta = delta * muladd(-(4.8 * delta), delta, 1.0)
            v = difmuk / fk
            if abs(v) <= 0.25:
                poly = muladd(0.1250060, v, -0.1384794)
                poly = muladd(poly, v, 0.1421878)
                poly = muladd(poly, v, -0.1661269)
                poly = muladd(poly, v, 0.2000118)
                poly = muladd(poly, v, -0.2500068)
                poly = muladd(poly, v, 0.3333333)
                poly = muladd(poly, v, -0.5)
                px = muladd(fk * v * v, poly, -delta)
            else:
                px = muladd(fk, math.log(1.0 + v), -difmuk) - delta
            py = M_1_SQRT_2PI / math.sqrt(fk)
        x = (0.5 - difmuk) / s
        x *= x
        fx = -0.5 * x
        fy = omega * muladd(muladd(muladd(c3, x, c2), x, c1), x, c0)
        if not quotient_step:
            if c * abs(u) <= muladd(py, math.exp(px + E), -(fy * math.exp(fx + E))):
                break
        else:
            quotient_step = False
            if muladd(-u, fy, fy) <= py * math.exp(px - fx):
                break
    return pois


@njit(cache=True, error_model="numpy")
def _rgeom(state, p):
    """One geometric number: `rgeom()` of the C library.

    Parameters: state, p (success probability). Returns: float.
    """
    if math.isnan(p) or math.isinf(p) or p <= 0 or p > 1:
        return np.nan
    return _rpois(state, _exp_rand(state) * ((1 - p) / p))


# ----- Vectors of draws, parameters recycled as in R -----

@njit(cache=True)
def _runif(state, n):
    out = np.empty(n)
    for i in range(n):
        out[i] = _unif(state)
    return out


@njit(cache=True)
def _rnorm(state, n, mean, sd):
    out = np.empty(n)
    for i in range(n):
        out[i] = muladd(sd, _norm_rand(state), mean)
    return out


@njit(cache=True)
def _rexp_std(state, n):
    out = np.empty(n)
    for i in range(n):
        out[i] = _exp_rand(state)
    return out


@njit(cache=True)
def _rbinom_n(state, n, size, prob):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rbinom(state, size[i % len(size)], prob[i % len(prob)])
    return out


@njit(cache=True)
def _rpois_n(state, n, mean):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rpois(state, mean[i % len(mean)])
    return out


@njit(cache=True)
def _rgeom_n(state, n, prob):
    out = np.empty(n)
    for i in range(n):
        out[i] = _rgeom(state, prob[i % len(prob)])
    return out


def _vector(x):
    return np.atleast_1d(np.asarray(x, dtype=np.float64))


def set_seed(seed):
    """`set.seed(seed)`: start the stream of this process at a seed.

    Parameters: seed (integer). Returns: None.
    """
    _seed(_state, int(seed) & 0xFFFFFFFF)


def runif(n):
    """`runif(n)`: n uniform numbers in (0, 1). Returns: float array."""
    return _runif(_state, n)


def rnorm(n, mean=0.0, sd=1.0):
    """`rnorm(n, mean, sd)`: n normal numbers. Returns: float array."""
    return _rnorm(_state, n, float(mean), float(sd))


def rexp(n):
    """`rexp(n)`: n standard exponential numbers. Returns: float array."""
    return _rexp_std(_state, n)


def rbinom(n, size, prob):
    """`rbinom(n, size, prob)`: n binomial numbers; size and prob are recycled.

    Returns: float array, NaN (R: NA) where a parameter is missing or out
    of range.
    """
    return _rbinom_n(_state, n, _vector(size), _vector(prob))


def rpois(n, mean):
    """`rpois(n, mean)`: n Poisson numbers. Returns: float array."""
    return _rpois_n(_state, n, _vector(mean))


def rgeom(n, prob):
    """`rgeom(n, prob)`: n geometric numbers. Returns: float array."""
    return _rgeom_n(_state, n, _vector(prob))

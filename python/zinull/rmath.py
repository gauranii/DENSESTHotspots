# -----------------------------------------------------------------------------
# File:    python/zinull/rmath.py
# Paper:   Gauran, Park, Lim, Park, Zylstra, Peterson, Kann, Spouge (2018).
#          Empirical Null Estimation Using Zero-Inflated Discrete Mixture
#          Distributions and its Application to Protein Domain Data.
#          Biometrics 74:458-471. All tables and figures.
# Purpose: The functions of R's math library that the code in R/ calls,
#          ported operation by operation so that they return the same doubles
#          as R: lgamma() and lfactorial(), the power operator `^`, and
#          log(), exp(), sqrt() as R answers them for any double (no Python
#          exception).
# Inputs:  None.
# Outputs: None.
# Usage:   from zinull.rmath import lgamma, r_pow
#          lgamma(3.5); r_pow(10.0, -10.0)
# Source:  R 4.6.1: src/nmath/lgamma.c (lgammafn), gamma.c (gammafn),
#          lgammacor.c, chebyshev.c (chebyshev_eval); src/main/arithmetic.c
#          (R_POW, R_pow, R_pow_di).
# Notes:   R's C code is compiled for arm64 with a product that feeds a sum
#          or a difference in the same expression done in one rounding (fused
#          multiply-add). Every such place is written here as muladd(); the
#          other products and sums are rounded one at a time, as in R.
# -----------------------------------------------------------------------------
import math
import platform

import numpy as np
from llvmlite import ir
from numba import njit
from numba.core import types
from numba.extending import intrinsic

# R is compiled with fused multiply-add on arm64 only.
FUSED = platform.machine().lower() in ("arm64", "aarch64")

if FUSED:

    @intrinsic
    def muladd(typingctx, a, b, c):
        """a * b + c with a single rounding (fused multiply-add).

        Parameters: a, b, c (floats). Returns: float.
        """
        sig = types.float64(types.float64, types.float64, types.float64)

        def codegen(context, builder, signature, args):
            double = ir.DoubleType()
            fnty = ir.FunctionType(double, [double, double, double])
            fn = builder.module.declare_intrinsic("llvm.fma", [double], fnty)
            return builder.call(fn, args)

        return sig, codegen

else:

    @njit(cache=True, error_model="numpy")
    def muladd(a, b, c):
        """a * b + c with two roundings, for machines where R does not fuse.

        Parameters: a, b, c (floats). Returns: float.
        """
        return a * b + c


DBL_EPSILON = 2.220446049250313e-16
DBL_MAX = 1.7976931348623157e308
M_LN2 = 0.693147180559945309417232121458
M_LN_SQRT_2PI = 0.918938533204672741780329736406
M_1_SQRT_2PI = 0.398942280401432677939946059934
M_LOG10_2 = 0.301029995663981195213738894724
LGAMMA_XMAX = 2.5327372760800758e305
GAMMA_XSML = 2.2474362225598545e-308
LGAMMACOR_XBIG = 94906265.62425156

_GAMCS = np.array(
    [
        +.8571195590989331421920062399942e-2,
        +.4415381324841006757191315771652e-2,
        +.5685043681599363378632664588789e-1,
        -.4219835396418560501012500186624e-2,
        +.1326808181212460220584006796352e-2,
        -.1893024529798880432523947023886e-3,
        +.3606925327441245256578082217225e-4,
        -.6056761904460864218485548290365e-5,
        +.1055829546302283344731823509093e-5,
        -.1811967365542384048291855891166e-6,
        +.3117724964715322277790254593169e-7,
        -.5354219639019687140874081024347e-8,
        +.9193275519859588946887786825940e-9,
        -.1577941280288339761767423273953e-9,
        +.2707980622934954543266540433089e-10,
        -.4646818653825730144081661058933e-11,
        +.7973350192007419656460767175359e-12,
        -.1368078209830916025799499172309e-12,
        +.2347319486563800657233471771688e-13,
        -.4027432614949066932766570534699e-14,
        +.6910051747372100912138336975257e-15,
        -.1185584500221992907052387126192e-15,
    ]
)
_ALGMCS = np.array(
    [
        +.1666389480451863247205729650822e+0,
        -.1384948176067563840732986059135e-4,
        +.9810825646924729426157171547487e-8,
        -.1809129475572494194263306266719e-10,
        +.6221098041892605227126015543416e-13,
    ]
)


# ----- log(Gamma(x)) -----

@njit(cache=True, error_model="numpy")
def _chebyshev_eval(x, a, n):
    """Chebyshev series with the first n coefficients of a: chebyshev_eval().

    Parameters: x (in [-1.1, 1.1]), a (array), n (number of terms).
    Returns: float.
    """
    if x < -1.1 or x > 1.1:
        return np.nan
    twox = x * 2
    b2 = 0.0
    b1 = 0.0
    b0 = 0.0
    for i in range(1, n + 1):
        b2 = b1
        b1 = b0
        b0 = muladd(twox, b1, -b2) + a[n - i]
    return (b0 - b2) * 0.5


@njit(cache=True, error_model="numpy")
def _gammafn_small(x):
    """Gamma(x) for 0 < x <= 10: that range of gammafn().

    Parameters: x (float in (0, 10]). Returns: float.
    """
    n = int(x)
    y = x - n  # in [0, 1)
    n -= 1
    value = _chebyshev_eval(muladd(y, 2.0, -1.0), _GAMCS, 22) + 0.9375
    if n == 0:
        return value
    if n < 0:  # 0 < x < 1
        if y < GAMMA_XSML:
            return np.inf
        return value / x
    for i in range(1, n + 1):
        value *= y + i
    return value


@njit(cache=True, error_model="numpy")
def _lgammacor(x):
    """The correction of Stirling's formula for x >= 10: lgammacor().

    Parameters: x (float, >= 10). Returns: float.
    """
    if x < 10:
        return np.nan
    if x < LGAMMACOR_XBIG:
        tmp = 10 / x
        return _chebyshev_eval(muladd(tmp * tmp, 2.0, -1.0), _ALGMCS, 5) / x
    return 1 / (x * 12)


@njit(cache=True, error_model="numpy")
def _lgammafn(x):
    """log(Gamma(x)) for x >= 0 as R computes it: lgammafn().

    Parameters: x (float, >= 0; this port is not for negative x).
    Returns: float.
    """
    if math.isnan(x):
        return x
    if x <= 0:
        return np.inf
    if x < 1e-306:
        return -math.log(x)
    if x <= 10:
        return math.log(abs(_gammafn_small(x)))
    if x > LGAMMA_XMAX:
        return np.inf
    if x > 1e17:
        return x * (math.log(x) - 1.0)
    if x > 4934720.0:
        return muladd(x - 0.5, math.log(x), M_LN_SQRT_2PI) - x
    return muladd(x - 0.5, math.log(x), M_LN_SQRT_2PI) - x + _lgammacor(x)


@njit(cache=True, error_model="numpy")
def _lgamma_array(x):
    out = np.empty(len(x))
    for i in range(len(x)):
        out[i] = _lgammafn(x[i])
    return out


def lgamma(x):
    """`lgamma(x)` for x >= 0. Parameters: x (float or array).
    Returns: float, or float array."""
    if np.ndim(x) == 0:
        return float(_lgammafn(float(x)))
    return _lgamma_array(np.asarray(x, dtype=np.float64))


def lfactorial(x):
    """`lfactorial(x)`, lgamma(x + 1). Parameters: x (float or array).
    Returns: float, or float array."""
    if np.ndim(x) == 0:
        return float(_lgammafn(float(x) + 1.0))
    return _lgamma_array(np.asarray(x, dtype=np.float64) + 1.0)


# ----- log() and exp() of any double -----

@njit(cache=True, error_model="numpy")
def _log(x):
    return math.log(x) if x >= 0.0 else np.nan


@njit(cache=True, error_model="numpy")
def _exp(x):
    return math.exp(x)


@njit(cache=True, error_model="numpy")
def _log_array(x):
    out = np.empty(len(x))
    for i in range(len(x)):
        out[i] = _log(x[i])
    return out


@njit(cache=True, error_model="numpy")
def _exp_array(x):
    out = np.empty(len(x))
    for i in range(len(x)):
        out[i] = math.exp(x[i])
    return out


def log(x):
    """`log(x)` of the C library: -Inf at zero, NaN below it.

    Parameters: x (float or array). Returns: float, or float array.
    """
    if np.ndim(x) == 0:
        return float(_log(float(x)))
    return _log_array(np.asarray(x, dtype=np.float64))


def exp(x):
    """`exp(x)` of the C library: Inf when it overflows.

    Parameters: x (float or array). Returns: float, or float array.
    """
    if np.ndim(x) == 0:
        return float(_exp(float(x)))
    return _exp_array(np.asarray(x, dtype=np.float64))


# ----- Powers -----

@njit(cache=True, error_model="numpy")
def _pow_di(x, n):
    """x^n for an integer n by repeated squaring: R_pow_di()."""
    xn = 1.0
    if math.isnan(x):
        return x
    if n != 0:
        if math.isinf(x):
            return math.pow(x, float(n))
        negative = n < 0
        if negative:
            n = -n
        while True:
            if n & 1:
                xn *= x
            n >>= 1
            if n:
                x *= x
            else:
                break
        if negative:
            xn = 1.0 / xn
    return xn


def pow_di(x, n):
    """x^n for an integer n as R's C code computes it (R_pow_di).

    Parameters: x (float), n (integer). Returns: float.
    """
    return float(_pow_di(float(x), int(n)))


@njit(cache=True, error_model="numpy")
def _r_pow(x, y):
    """x^y as R's `^` computes it for doubles: R_POW() and R_pow()."""
    if y == 2.0:
        return x * x
    if x == 1.0 or y == 0.0:
        return 1.0
    if x == 0.0:
        if y > 0.0:
            return 0.0
        if y < 0.0:
            return np.inf
        return y
    return math.pow(x, y)


def r_pow(x, y):
    """`x^y` of R for two doubles: x * x for y = 2, 1 for x = 1 or y = 0,
    else pow() of the C library.

    Parameters: x, y (floats). Returns: float.
    """
    return float(_r_pow(float(x), float(y)))


# ----- Square root -----

@njit(cache=True, error_model="numpy")
def _sqrt(x):
    return math.sqrt(x) if x >= 0.0 else np.nan


def sqrt(x):
    """`sqrt(x)`: NaN below zero. Parameters: x (float). Returns: float."""
    return float(_sqrt(float(x)))

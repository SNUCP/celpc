from math import sqrt, pi, log2, log, ceil
import sympy

"""
Not Compressed
N: 19 21 23 25
n: 12 13 14 15
"""

N = 2**25
n = 2**15

"""
r = 32
b = 156
"""

r = 16
b = 63388

d = 2**11
ell = r * n / d
m = N / n

D = 24

print("N: %d n: %d m: %d ell: %d" % (N, n, m, ell))

delta = 1.005
eta = sqrt(log(2 * (2**30) * (2**128)) / pi)
kp = ceil(128 / log2(2 * d))

s1 = sqrt(3) * eta
s2 = sqrt(3) * sqrt(kp) * eta
s3 = sqrt(3) * (b * r / 2) * eta
print("s1: %d s2: %d s3: %d" % (s1, s2, s3))

sig1 = 2 * sqrt(3) * eta
sig2 = 2 * sqrt(3) * sqrt(kp) * eta
sig3 = 2 * sqrt(3) * (b * r / 2) * eta
print("sig1: %d sig2: %d sig3: %d" % (sig1, sig2, sig3))

mu = ceil(2048 / d)
nu = ceil(4096 / d)

q = 2**112 + 1


beta_max = 2 ** (2 * sqrt(mu * d * log2(q) * log2(delta)))

beta_open = ((m + 1) * sig1 + sqrt(m + 2) * sig2) ** 2 * nu * d
beta_open += (
    ((m + 1) * sig1 + sqrt(m + 2) * sig2 + (m + 1) * pow(2, D - 1)) ** 2 * mu * d
)
beta_open += (b + 1) ** 2 * ((m + 1) * s1 + sqrt(m + 2) * s2) ** 2 * ell * d
beta_open = beta_open**0.5

beta_eval = ((m + 1) * (b + 1) * (r / 2) * sig1 + sqrt(m + 2) * sig3) ** 2 * nu * d
beta_eval += (
    (
        (m + 1) * (b + 1) * (r / 2) * sig1
        + sqrt(m + 2) * sig3
        + (m + 1) * (b + 1) * (r / 2) * pow(2, D - 1)
    )
    ** 2
    * mu
    * d
)
beta_eval += (
    (b + 1) ** 2 * ((m + 1) * (b + 1) * (r / 2) * s1 + sqrt(m + 2) * s2) ** 2 * ell * d
)
beta_eval = beta_eval**0.5

beta_pc = beta_eval + (b + 1) * (m + 1) * (d * r / 2) * beta_open
beta_pc = 4 * beta_pc

print("log2(beta_open): ", log2(beta_open))
print("log2(beta_eval): ", log2(beta_eval))
print("log2(beta_pc): ", log2(beta_pc)),
print("log2(beta_max): ", log2(beta_max))


commit_size = (log2(q) - D) * (m + kp + 2) * mu * d

open_size = (
    kp
    * d
    * mu
    * log2(5 * (m + 1) * sig1 + 5 * sqrt(m + 2) * sig2 + (m + 1) * pow(2, D - 1))
)
open_size += kp * d * nu * log2(5 * (m + 1) * sig1 + 5 * sqrt(m + 2) * sig2)
open_size += kp * d * ell * log2(5 * b * (m + 1) * s1 + 5 * b * sqrt(m + 2) * s2)

eval_size = (
    d
    * mu
    * log2(
        5 * (m + 1) * b * (r / 2) * sig1
        + 5 * sqrt(m + 2) * sig3
        + (m + 1) * b * (r / 2) * pow(2, D - 1)
    )
)
eval_size += d * nu * log2(5 * (m + 1) * b * (r / 2) * sig1 + 5 * sqrt(m + 2) * sig3)
eval_size += (
    d * ell * log2(5 * (m + 1) * (b**2) * (r / 2) * s1 + 5 * b * sqrt(m + 2) * s3)
)


total_size = (commit_size + open_size + eval_size) / (2**23)
print("commit size: ", commit_size / (2**23), " MB")
print("open size: ", open_size / (2**23), " MB")
print("eval size: ", eval_size / (2**23), " MB")
print("total size: ", total_size, " MB")

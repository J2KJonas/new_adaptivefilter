# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 2021

@author: Peter Kjær Fisker (pfiske15@student.aau.dk)
         with heavy inspiration from examples by Jan Østergaard (jo@es.aau.dk)
"""

from scipy.signal import lfilter, firwin
from scipy.linalg import hankel
import numpy as np
from numpy.random import randn
import matplotlib.pyplot as plt

# Simulation setup

# Signal parameters and statistics
M = 100                                             # length of the filter
wo = firwin(M, 0.2, window='hamming')               # optimal filter coefficients
sigv = 0.01                                         # observation noise variance
sigu = 1                                            # variance of the white input signal
Ru = sigu * np.eye(M)                               # correlation matrix
rud = np.matmul(Ru, wo)                             # cross-correlation vector
sigd = sigv + np.matmul(np.matmul(wo.T, Ru), wo)    # power of the desired signal
Jwo = sigv                                          # minimum of the cost function

# Simulation parameters
N = 10000                                           # number of iterations

# Parameters of the adaptive filter
w = np.zeros((M, N + 1))                            # estimated filter coefficients
e = np.zeros(N)                                     # errors
y = np.zeros(N)                                     # output of the filter
Leff = 500                                          # effective window length
lam = 1 - 1 / Leff                                  # forgetting factor
delt = sigu                                         # variance of initial correlation matrix
P = np.eye(M) / delt                                # P(0)

# Generate the signals
u = np.sqrt(sigu) * randn(N + M - 1)                # white gaussian noise
U = np.flipud(hankel(u[0:M], u[M - 1:N + M]))       # matrix of input vectors
z = lfilter(wo, 1.0, u)                             # create noiseless desired signal
d = z[M - 1:N + M] + np.sqrt(sigv) * randn(N)       # desired vector

# ------------------------------------------
# Run the adaptive filter
for n in range(N):
    z = np.matmul(P, U[:, n])                       # temporary variable
    k = z / (lam + np.dot(U[:, n].T, z))            # gain vector
    y[n] = np.dot(U[:, n].T, w[:, n])               # output of the filter
    e[n] = d[n] - y[n]                              # a priori error
    w[:, n + 1] = w[:, n] + k * e[n]                # update the filter coefficients
    P = (P - np.outer(k, z.T))/lam                  # update inverse correlation matrix


plt.stem(wo, label="Actual", basefmt='b')
plt.stem(w[:,-1], markerfmt='C1o', label="Estimated", basefmt='b')
plt.legend()
plt.title(f'Forgetting factor: {lam}, Filter order: {M}')
plt.show()
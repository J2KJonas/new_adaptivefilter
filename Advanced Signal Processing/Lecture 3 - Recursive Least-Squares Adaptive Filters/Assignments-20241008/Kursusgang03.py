# -*- coding: utf-8 -*-
"""
Created on Fri Sep 10 2021

@author: Peter Kjær Fisker (pfiske15@student.aau.dk)
         with heavy inspiration from examples by Jan Østergaard (jo@es.aau.dk)
"""

import numpy as np
import soundfile as sf
from time import time

# Load Data Files
x, fs1 = sf.read("Signals_MM3/Music.wav")
d, fs2 = sf.read("Signals_MM3/Noise.wav")
u, Fs = sf.read("Signals_MM3/Noisy_Music.wav")

N = len(u)

# Filter order
M = 12

# LMS
t_start = time()
mu = 0.01
w = np.random.randn(M)
padded_u = np.pad(u, (M - 1, 0), 'constant')
y = np.zeros(N)

for n in range(0, N):
    u_vect = padded_u[n:n + M]
    e = d[n] - np.dot(w.T, u_vect)
    w = w + mu * e * u_vect
    y[n] = np.dot(w.T, u_vect)

filtered_signal_LMS = u - y
time_LMS = time() - t_start

SNR_LMS = 10 * np.log10((np.sum(x ** 2)) / (np.sum((filtered_signal_LMS - x) ** 2)))

# -----------------------------------------------------------------
# NLMS
t_start = time()
mu = 1
w = np.random.randn(M)
padded_u = np.pad(u, (M - 1, 0), 'constant')
y = np.zeros((N,))
Eps = 0.0001

for n in range(N):
    u_vect = padded_u[n:n + M]
    mu1 = mu / (Eps + np.linalg.norm(u_vect, 2) ** 2)
    e = d[n] - np.dot(w.T, u_vect)
    w = w + mu1 * e * u_vect
    y[n] = np.dot(w.T, u_vect)

filtered_signal_NLMS = u - y

time_NLMS = time() - t_start

SNR_NLMS = 10 * np.log10((np.sum(x ** 2)) / (np.sum((filtered_signal_NLMS - x) ** 2)))

# ------------------------------------------------
# RLS
t_start = time()
lamda_1 = 1 - 1 / (0.1 * M)
delta = 0.01

P = 1 / delta * np.eye(M)
w = np.random.randn(M)
padded_u = np.pad(u, (M - 1, 0), 'constant')
y = np.zeros((N,))

for n in range(N):
    u_vect = padded_u[n:n + M]
    PI = np.matmul(P, u_vect)
    gain_k = PI / (lamda_1 + np.dot(u_vect.T, PI))
    prior_error = d[n] - np.dot(w.T, u_vect)
    w = w + prior_error * gain_k

    # Because soundfile uses (n,) as shape for the nx1 vector we do some reshaping-hax
    temp1 = np.matmul(u_vect.T, P).reshape((1, M))
    temp2 = np.matmul(gain_k.reshape((M, 1)), temp1)

    P = P / lamda_1 - temp2 / lamda_1
    y[n] = np.dot(w.T, u_vect)

filtered_signal_RLS = u - y

time_RLS = time() - t_start
SNR_RLS = 10 * np.log10((np.sum(x ** 2)) / (np.sum((filtered_signal_RLS - x) ** 2)))

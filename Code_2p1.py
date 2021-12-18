# -*- coding: utf-8 -*-
"""
Created on Sun Aug  4 17:26:20 2019

@author: msahr
"""

import numpy as np

Niter=50
h = np.zeros((Niter, Niter))
A = [[3,2,1],[1,3,2],[2,1,3]]
b = [7, 5, 12]
x0 = [1, 0, 0]

r = b - np.asarray(np.matmul(A, x0)).reshape(-1)
x = []
v = [0 for i in range(Niter)]

x.append(r)

v[0] = r / np.linalg.norm(r)

for i in range(Niter):
    w = np.asarray(np.matmul(A, v[i])).reshape(-1)

    for j in range(i):
        h[j, i] = np.matmul(v[j], w)
        w = w - h[j, i] * v[j]
    if i < Niter-1 :
        h[i + 1, i] = np.linalg.norm(w)
        if (h[i + 1, i] != 0 and i != Niter - 1):
            v[i + 1] = w / h[i + 1, i]

    b = np.zeros(Niter)
    b[0] = np.linalg.norm(r)

    ym = np.linalg.lstsq(h, b,rcond=None)[0]
    x.append(np.dot(np.transpose(v), ym) + x0)


print(x)

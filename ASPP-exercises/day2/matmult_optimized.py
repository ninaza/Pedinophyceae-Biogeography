#!/usr/bin/env python3

import random

N = 250

# Create NxN matrix X using list comprehension
X = [[random.randint(0, 100) for _ in range(N)] for _ in range(N)]

# Create Nx(N+1) matrix Y using list comprehension
Y = [[random.randint(0, 100) for _ in range(N + 1)] for _ in range(N)]

# Initialize result matrix using list comprehension
result = [[0] * (N + 1) for _ in range(N)]

# Perform matrix multiplication
for i in range(N):
    for j in range(N + 1):
        for k in range(N):
            result[i][j] += X[i][k] * Y[k][j]

# Print result
for r in result:
    print(r)


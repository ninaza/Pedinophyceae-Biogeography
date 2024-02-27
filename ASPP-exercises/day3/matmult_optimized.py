import numpy as np

N = 250

# Create random matrices using NumPy
X = np.random.randint(0, 100, size=(N, N))
Y = np.random.randint(0, 100, size=(N, N + 1))

# Perform matrix multiplication
result = np.matmul(X, Y)

for r in result:
    print(r)
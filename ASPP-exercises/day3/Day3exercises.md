# Day3 exercises

## 2. Numpy

### a. 

`nv = np.zeros(10)
nv[4] = 1
print(nv)
[0. 0. 0. 0. 1. 0. 0. 0. 0. 0.]`


### b.

`v = np.arange(10, 50)
print(v)
[10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33
 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49]`

 ### c.

`rv = v[::-1]
print(rv)
[49 48 47 46 45 44 43 42 41 40 39 38 37 36 35 34 33 32 31 30 29 28 27 26
 25 24 23 22 21 20 19 18 17 16 15 14 13 12 11 10]`

 ### d.

 `m = np.arange(9).reshape(3,3)
 print(m)
[[0 1 2]
 [3 4 5]
 [6 7 8]]`

 ### e.

`a = np.array([1,2,0,0,4,0])
non_zero_indices = np.nonzero(a)
print(non_zero_indices)
(array([0, 1, 4]),)`

### f.

`v = np.random.random((30))
v.mean()`

### g.

`a2d = np.ones(25).reshape(5,5)
a2d[1:-1, 1:-1] = 0
print(a2d)
[[1. 1. 1. 1. 1.]
 [1. 0. 0. 0. 1.]
 [1. 0. 0. 0. 1.]
 [1. 0. 0. 0. 1.]
 [1. 1. 1. 1. 1.]]`

### h.

`m = np.zeros((8, 8)
m[1::2, ::2] = 1
m[::2, 1::2] = 1
print(m)
[[0. 1. 0. 1. 0. 1. 0. 1.]
 [1. 0. 1. 0. 1. 0. 1. 0.]
 [0. 1. 0. 1. 0. 1. 0. 1.]
 [1. 0. 1. 0. 1. 0. 1. 0.]
 [0. 1. 0. 1. 0. 1. 0. 1.]
 [1. 0. 1. 0. 1. 0. 1. 0.]
 [0. 1. 0. 1. 0. 1. 0. 1.]
 [1. 0. 1. 0. 1. 0. 1. 0.]]`

 ### i.

`m = np.array([[0, 1], [1, 0]])
m_8x8 = np.tile(m, (4, 4))
print(m_8x8)
[[0 1 0 1 0 1 0 1]
 [1 0 1 0 1 0 1 0]
 [0 1 0 1 0 1 0 1]
 [1 0 1 0 1 0 1 0]
 [0 1 0 1 0 1 0 1]
 [1 0 1 0 1 0 1 0]
 [0 1 0 1 0 1 0 1]
 [1 0 1 0 1 0 1 0]]`

### j.

`Z = np.arange(11)
i = np.array([3,4,5,6,7,8])
Z[i] *= -1
print(Z)`

### k.

`Z = np.random.random(10)
Z = np.sort(Z)
print(Z)`

### l.

`A = np.random.randint(0,2,5)
B = np.random.randint(0,2,5)
equal = np.array_equal(A, B)
print(equal)`

### m.

`Z = np.arange(10, dtype=np.int32)
print(Z)
[0 1 2 3 4 5 6 7 8 9]
print(Z.dtype)
int32
Z **= 2
print(Z)
[ 0  1  4  9 16 25 36 49 64 81]
print(Z.dtype)
int32`

### n.

`A = np.arange(9).reshape(3,3)
B = A + 1
C = np.dot(A,B)
D = np.diag(C)
print(D)
[ 18  66 132]
`
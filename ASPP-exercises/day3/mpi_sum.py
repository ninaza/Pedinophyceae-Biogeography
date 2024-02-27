#!/usr/bin/env python3

import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

data = np.random.randint(0, 100, size=1)

total_sum = comm.reduce(data, op=MPI.SUM, root=0)

if rank == 0:
    print("Total sum:", total_sum)
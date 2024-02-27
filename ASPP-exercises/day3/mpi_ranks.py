#!/usr/bin/env python3

from mpi4py import MPI

comm = MPI.COMM_WORLD

rank = comm.Get_rank()

print("This is the process with rank ", rank)
      
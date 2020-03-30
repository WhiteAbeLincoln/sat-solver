#!/usr/bin/env bash

# Written by Abe White and Jonathan Petersen
# Parallel Computing - CS 5500 / 6500

NumProcessors=4

mpirun -np $NumProcessors Solver.out

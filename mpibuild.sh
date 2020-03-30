#!/usr/bin/env bash

# Written by Abe White and Jonathan Petersen
# Parallel Computing - CS 5500 / 6500

# Configure Build
SOURCE_FILES=(
    src/main.cpp
    src/mpi/Communicator.cpp
)

# Run Build
mpic++ -std=c++17 -o Solver.out -I src ${SOURCE_FILES[*]}

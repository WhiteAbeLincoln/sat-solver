// Written by Abe White and Jonathan Petersen
// Parallel Computing - CS 5500 / 6500

#include <iostream>

#include "mpi/Communicator.hpp"

/**
 * Program Entry Point.
 * @param argc Number of arguments.
 * @param argv Array of arguments.
 * @returns 0 on successful execution, error code otherwise.
 */
int main(int argc, char* argv[])
{
  auto c = mpi::Communicator(argc, argv);

  std::cout << "Hello World from processor " << c.m_myRank << "!" << std::endl;

  return 0;
}

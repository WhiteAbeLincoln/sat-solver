// Written by Abe White and Jonathan Petersen
// Parallel Computing - CS 5500 / 6500

#include "mpi/Communicator.hpp"

using mpi::Communicator;

const MPI_Comm Communicator::COMMUNICATOR = MPI_COMM_WORLD;

Communicator::Communicator(int argc, char* argv[]) : m_myRank(0), m_size(0)
{
  // Initialize MPI
  MPI_Init(&argc, &argv);

  // Setup My Actor Info
  MPI_Comm_rank(COMMUNICATOR, &m_myRank);
  MPI_Comm_size(COMMUNICATOR, &m_size);

  // Setup Destinations
  m_nextHighest = (m_myRank + 1) % m_size;
  m_nextLowest = (m_myRank - 1) % m_size;
}

Communicator::~Communicator()
{
  // Cleanup MPI
  MPI_Finalize();
}

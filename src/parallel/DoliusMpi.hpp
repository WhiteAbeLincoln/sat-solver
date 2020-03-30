// Written by Abe White and Jonathan Petersen
// Parallel Computing - CS 5500 / 6500

#ifndef DOLIUS_MPI_HPP
#define DOLIUS_MPI_HPP

#include <fstream>
#include <string>
#include <vector>

#include "mpi/Communicator.hpp"

namespace parallel
{
enum class SolutionType
{
  UNKNOWN,
  UNSAT,
  SAT
};

class Dolius
{
public:
  using Clause = std::vector<int>;
  using GuidingPath = std::vector<Clause>;
  using LearntClauseIterator = Clause::iterator;
  using Solver = int; // TODO: Replace with our solver type once we have one.

  // Initialization
  Dolius(int numVariables, int numClauses, int argc, char* argv[]);
  ~Dolius();

  Dolius(const Dolius& other) = delete;
  Dolius(Dolius&& other) = default;

  Dolius& operator=(const Dolius& other) = delete;
  Dolius& operator=(Dolius&& other) = default;

  void setCNFFile(std::string filename);

  // Search Controls
  void run();
  void stop();

  // Clause Database Modification
  void addLearntClause(const Clause& clause);
  void addClause(const Clause& clause);

  // Iterators
  LearntClauseIterator getClauseIterator();

  // Guiding Path
  bool createGuidingPath(GuidingPath& a, GuidingPath& b);
  void addToGuidingPath(const Clause& clause);

  // SAT Information
  bool solutionFound() const;
  bool isSolutionFoundSAT() const;
  int getNbVar() const;
  int getSolutionLiteral(int var) const;
  int getNumLeartClauses() const;

private:
  // Solver m_solver;
  bool m_active;
  GuidingPath m_guidingPath;
  int m_numVariables;
  mpi::Communicator m_comm;
  SolutionType m_status;
  std::ifstream m_cnfFile;
};

} // namespace parallel

#endif

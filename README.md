# sat-solver
Serial and Parallel implementations of an SAT solver for CS5500

Test with the runtest.sh file. Expects a list of test files in DIMACS format, with optional parameters -t to specify
the command under test, -c to specify a canonical implementaion to compare with, and -s or -u to specify that the
result should be SAT or UNSAT respectively.

The command under test should exit with code 10 when SAT or code 20 when UNSAT.

Currently not all branching rule heuristic functions can handle the whole test suite. I don't know if this is a bug
or just inherent to the heuristic. Consequently, it might make more sense to use a portfolio parallel solution rather
than the divide-and-conquer.

# sat-solver
Serial and Parallel implementations of an SAT solver for CS5500

Test with the runtest.sh file. Expects a list of test files in DIMACS format, with optional parameters -t to specify
the command under test, -c to specify a canonical implementaion to compare with, and -s or -u to specify that the
result should be SAT or UNSAT respectively. Prints a tab separated list of the failed files and exit codes at the end of 
execution.

Both commands should exit with code 10 when SAT or code 20 when UNSAT.


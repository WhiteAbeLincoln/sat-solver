Parallel SAT Solver
===================
_By Abe White and Jonathan Petersen_

This project solves the [Boolean Satisfiability Problem](https://en.wikipedia.org/wiki/Boolean_satisfiability_problem) in a parallel manner using the [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface) library.

Project Layout
--------------
- `data` - Input files containing problems go here.
- `src` - Source code goes here.
- `.vscode` - Contains VS Code specific configuration. Can be ignored on other IDEs.

- `mpibuild.sh` and `mpirun.sh` contain the build and run scripts for the project.
- `.clang-format` can be used with [clang-format](https://clang.llvm.org/docs/ClangFormat.html) to automatically format the code.

# Overview

Generates the Doxygen documentation at the following URL:
-	https://synthetic-methods.github.io/xtal--math

Invokes `doxygen/index.sh` via `.github/workflows/pages.yml`:
-	Generating `exclude/xtal--math/all.hh` via `clang -E -C -P`.
-	Generating `exclude/xtal--math/html` via `doxygen`.

# Note

The customized CSS/HTML wrangling is ad-hoc,
pending a more comprehensive adadptation of Doxygen's templates.

# To Do...

Fold this branch into `main` once the kinks have been resolved.

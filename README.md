pbmpi
=====

phylobayes mpi

Performance work
----------------

Speed optimizations that follow exactly the same Markov chain as the original code (same data, model, seed and number of MPI ranks). Start with `docs/PERFORMANCE.md`.

- `docs/PERFORMANCE.md`: results, history, scaling, limits, how to reproduce
- `docs/code-optimizations.md`: each code change and why the results are unchanged
- `docs/build-and-compiler-options.md`: compiler flags (`make ARCHFLAGS=-mavx2`) and their effect
- `docs/pinning-and-cpu-affinity.md`: run one MPI rank per physical core (about 1.5x faster)
- `docs/validation.md`: how identity with the original is tested, and what is not tested

# Build and compiler options

> Part of the performance documentation; see the scope note at the top of `PERFORMANCE.md`.

This note lists every build option that matters for speed or for results, what was measured, and what is recommended.

## 1. Summary

- The default build (`make` in `sources/`) uses `-Wall -O3 -std=c++11`. It runs on any x86-64 CPU and on any compiler that can build the original code.
- For the best speed on a CPU with AVX2, build with `make ARCHFLAGS=-mavx2`. This changes only the machine code, not the arithmetic. On `master` and on the optimized code, chains from the `-mavx2` build are identical to chains from the default build.
- **Never** add `-ffast-math`, and never add `-mfma` or `-march=native` without `-ffp-contract=off`. They can change results.
- The measured gain of `-mavx2` is small and noisy: about 1.06x on the current kernel (unpinned runs). It was larger (about 1.1x to 1.3x) before the kernel was rewritten. See section 4.
- Other flags tested (`-mno-avx256-split-unaligned-load`, `-mno-avx256-split-unaligned-store`, `-mtune=znver1`) gave no gain that could be told apart from noise. Those tests were unpinned. A pinned re-test has not been done.
- `-DPB_MPI_STATS` adds per-message timing on the slaves. It is for diagnosis only.
- `-DPB_NO_VECTOR_KERNEL` switches off the hand-written vector kernel and uses plain loops.

## 2. How to build

```
cd sources
make                                  # default flags, builds ../data/pb_mpi and the other tools
make ARCHFLAGS=-mavx2                 # same, plus AVX2 machine code
make ARCHFLAGS="-mavx2 -DPB_MPI_STATS"   # AVX2 and per-message statistics
```

`ARCHFLAGS` is appended to `CPPFLAGS` (`CPPFLAGS= -Wall -O3 -std=c++11 $(ARCHFLAGS)`). The compiler wrapper is `mpic++`. To override the whole flag set, pass `CPPFLAGS=...` on the command line.

`scripts/install_build.sh REF [ARCHFLAGS] [SUFFIX]` builds a git ref into `builds/<shortsha>-<suffix>/` and writes `BUILD_INFO` there. `BUILD_INFO` records the commit, flags, compiler, the exact MPI module name, the number of `%ymm` (AVX) instructions found in the binary, the build host and the regression runs that the commit passed. Use it for production builds, so that every chain can be traced to a known binary.

Environment used for every measurement in this repository: GCC 8.5.0 (Red Hat 8.5.0-30), Open MPI 4.1.2 (`openmpi/4.1.2_slurm-24.11.1_mpi1-compat`), AMD EPYC 7713 (Zen 3).

## 3. Options and what they do

| Option | Changes machine code | Can change results | Measured effect | Recommendation |
|---|---|---|---|---|
| `-O3` (default) | yes | no (no value-changing optimizations without `-ffast-math`) | baseline | keep |
| `-std=c++11` (default) | no | no | none | keep |
| `ARCHFLAGS=-mavx2` | yes: 4-wide vector instructions | no (see below) | about 1.06x to 1.28x, noisy | use on AVX2 CPUs |
| `-mavx2 -mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store` | yes | no | not separable from noise (unpinned) | do not adopt yet |
| `-mtune=znver1` | yes (scheduling only) | no | no gain (unpinned) | do not use |
| `-mtune=znver2` | n/a | n/a | does not build with GCC 8.5 (`bad value`, the highest value is `znver1`) | needs GCC 9 or later |
| `-mfma`, `-march=native` | yes | **yes** unless `-ffp-contract=off` | not measured | avoid |
| `-ffast-math` | yes | **yes** | not measured | never |
| `-g -fno-omit-frame-pointer` | no change to the arithmetic | no | used by `bench.sbatch` so `perf` can unwind stacks | profiling only |
| `-DPB_MPI_STATS` | adds timing calls on slaves | no | overhead not measured | diagnosis only |
| `-DPB_NO_VECTOR_KERNEL` | switches the kernel back to plain loops | no | not measured | for debugging, or odd compilers |

### 3.1 Why `-mavx2` does not change the result

The main kernel (`Propagate.cpp`, see `code-optimizations.md`) computes matrix-by-vector products. For each output element, it adds the terms in a fixed order and uses separate multiply and add instructions. The compiler may process 4 output elements at once with AVX2, but each element still receives the same sequence of operations. No sum is split or reordered. IEEE-754 multiplication and addition give the same value in a vector register as in a scalar register. So the result is bit-identical.

Evidence, all from `compare_runs.sh` (identical `.treelist`, `.chain`, `.trace` except timing columns):

- `master`, default flags against `master` with `-mavx2`: identical (38 taxa x 59,994 sites, 64 ranks, 8 cycles).
- 312c2bd, default flags against 312c2bd with `-mavx2`: identical (40x5000, 8 ranks, 20 cycles).
- `master`, default flags against 312c2bd with `-mavx2` and pinned ranks: identical (same data, 20 cycles).
- Three earlier pairs of commits (with and without `-mavx2`) also gave identical chains.

### 3.2 Why FMA is dangerous

A fused multiply-add computes `a*b + c` with one rounding instead of two. It gives a slightly different value. The vector kernel writes `a += x * y`. If the compiler is allowed to contract this into an FMA, the sum changes in the last bits, and a long chain would then differ from `master`.

- `-mavx2` alone does not enable FMA instructions. The compiled kernel contains no FMA instruction (checked by disassembly of an earlier build, which contained zero `vfmadd` instructions).
- `-mfma` and `-march=native` on Zen 3 do enable FMA. GCC's default for contraction in GNU C++ is `-ffp-contract=fast`. This was not verified for this compiler version in strict ISO mode. Do not rely on the default: when you use `-mfma` or `-march=native`, add `-ffp-contract=off` explicitly.
- The `Makefile` comment says the same.
- Whether a build with `-mfma -ffp-contract=off` is still identical was not tested.

### 3.3 The "split unaligned" flags

GCC 8.5 splits each unaligned 256-bit vector load or store into two 128-bit halves by default for generic tuning. The flags `-mno-avx256-split-unaligned-load` and `-mno-avx256-split-unaligned-store` turn that off. Disassembly of an AVX2 build shows many `vinsertf128` and `vextractf128` instructions, which is the splitting.

Measurements on 8 unpinned ranks (CAT-GTR, 40x5000, 3 runs each, wall seconds):

| Code | `-mavx2` | `-mavx2` plus the two flags |
|---|---|---|
| commit 51fdf5f (before the kernel rewrite) | 278, 266, 245 | 227, 245, 192 |
| commit 312c2bd (current kernel) | 151, 161, 163 | 145, 136, 162 |

The means favour the two flags (221 s against 263 s, and 148 s against 158 s), but the spread inside each group (up to 60 s) is as large as the difference. The runs were unpinned, and unpinned runs vary a lot (see `pinning-and-cpu-affinity.md`). **No conclusion can be drawn.** The flags do not change results (identical chains on the earlier commit). A pinned test with at least 3 runs per variant is the open item before adopting them.

### 3.4 `-mtune`

`-mtune=znver2` is rejected by GCC 8.5. `-mtune=znver1` built and gave identical chains, with no speed gain (unpinned: mean 276 s against 263 s for plain `-mavx2` on the older code).

### 3.5 `-DPB_MPI_STATS`

With this macro, each slave records, for each message type, the number of messages, the time spent executing them, and the idle time after them. It writes `mpistats.<rank>.txt` in the current directory when the program ends normally. The last line (`TOTAL`) gives the totals per rank. Use it to find load imbalance (for example, see `pinning-and-cpu-affinity.md`). The macro has no effect on the chain. Its cost in run time was not measured separately.

### 3.6 `-DPB_NO_VECTOR_KERNEL`

Switches `Propagate.cpp` to the simple column-by-column loops. Use it to test whether a problem comes from the vector kernel, or to build on a compiler that does not support GNU vector extensions (this happens automatically for compilers that do not define `__GNUC__`).

## 4. Measured effect of `-mavx2`

All runs: CAT-GTR, 40x5000, 8 ranks, unpinned, 20 cycles, wall seconds. Three runs each. Noise between runs of the same build is 6% to 21% in this data set, so only large differences mean something.

| Commit | Default flags (mean) | `-mavx2` (mean) | Ratio |
|---|---|---|---|
| 31a9b87 | 265.5 | 229.9 | 1.15x |
| 312c2bd | 170.3 | 160.6 | 1.06x |
| 14297a9 (single run each) | 335.0 | 302.6 | 1.11x |
| a183303 (single run each) | 303.6 | 237.4 | 1.28x |

The effect shrinks as the code gets faster. In the current kernel, the vector work is written with explicit vector types (see `code-optimizations.md`). Without `-mavx2` the compiler uses 2-wide SSE2 vectors. With `-mavx2` it uses 4-wide AVX vectors. These numbers have not been repeated with pinned ranks.

## 5. Portability

The vector kernel in `Propagate.cpp` is selected at compile time:

| Condition | Kernel |
|---|---|
| `__GNUC__` defined (GCC, clang) and `-DPB_NO_VECTOR_KERNEL` not given | GNU vector extensions (`__attribute__((vector_size(...)))`) |
| `__AVX__` defined (`-mavx` or `-mavx2`) | 4 doubles per vector, 5 accumulators for a 20-state block |
| `__AVX__` not defined | 2 doubles per vector, 10 accumulators |
| any other compiler, or `-DPB_NO_VECTOR_KERNEL` | plain loops, same order of operations |

Notes:

- Vector extensions are a GCC and clang feature. Compilers that do not define `__GNUC__` use the plain loops. The plain loops give the same results.
- The code uses compound-literal syntax `(pbvec){a,b,...}`, which is a GNU extension for C++. GCC and clang accept it. Other compilers that define `__GNUC__` without supporting it (this is unusual) would fail to build; use `-DPB_NO_VECTOR_KERNEL` for them.
- On non-x86 CPUs (for example ARM), the 16-byte vector type maps to the CPU's native vector unit or is split into scalars by the compiler. This was **not tested**. Results should still be identical as long as the compiler does not fuse multiply and add. On ARM, `-ffp-contract=off` may be needed because FMA is part of the base instruction set there. This is untested.
- Built and tested only with GCC 8.5.0. Other GCC versions and clang were not tested. After changing the compiler, run `scripts/bench/regression.sh` (see `validation.md`).
- An AVX2 binary stops with an illegal-instruction error on a CPU without AVX2. This is standard behavior and was not tested here. On a mixed cluster, build one binary per CPU class. `install_build.sh` records the number of `%ymm` instructions in `BUILD_INFO` so a build can be checked.

## 6. Recommendations

1. For production on AVX2 nodes: `make ARCHFLAGS=-mavx2`, or `scripts/install_build.sh <ref>` (it defaults to `-mavx2`).
2. For a portable binary: `make` with no `ARCHFLAGS`.
3. Always run pinned (see `pinning-and-cpu-affinity.md`). Pinning gives about 1.5x, more than any flag tested.
4. Do not use `-ffast-math`. Do not use `-mfma` or `-march=native` without `-ffp-contract=off`, and re-run the regression if you do.
5. After any change of compiler, compiler version, CPU type or flags, run the regression (`validation.md`).

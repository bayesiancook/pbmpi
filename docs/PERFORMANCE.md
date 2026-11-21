# Performance work on PhyloBayes-MPI: results and how to reproduce them

> **Scope.** These documents describe the performance work as a whole. It is proposed as three pull requests (pruning and SPR updates; likelihood kernels; caches for the mixture moves), each on top of the stepping-stone fixes of pull request 39. Commit names, branch names other than `master` (such as `perf-opt`), and the scripts under `scripts/` refer to the development repository `hyphaltip/pbmpi` (branch `perf-opt`); the scripts are not part of the pull requests.

This is the entry point for the performance documentation of the `perf-opt` branch. It states what was changed, what speed-up was measured, how reliable the numbers are, and where to read more.

| Document | Content |
|---|---|
| `PERFORMANCE.md` (this file) | results, history, scaling, limits |
| `code-optimizations.md` | each code change: what, why, why the result is unchanged, risks |
| `build-and-compiler-options.md` | compiler flags, what they do, what was measured |
| `pinning-and-cpu-affinity.md` | why one rank per physical core matters; how to do it |
| `validation.md` | the identical-results rule, tests, results, gaps |

## 1. Summary

- The optimized code (`perf-opt`; the speed changes end at commit 312c2bd, and the later commits are scripts, documentation and the two corrections in pull request 39) follows **exactly the same Markov chain as `master`** for the same data, model, seed and number of ranks. Tested on 14 model configurations (for the two stepping-stone ones, apart from two output columns that `master` prints as `-nan`; see pull request 39). `validation.md` gives the evidence and the gaps.
- Speed-up in wall time, measured against `master` in the same batch of jobs (8 ranks, pinned, default flags, two runs each): **4.2x to 5.0x** for CAT-GTR, up to **5.5x** for LG, **2.1x to 3.4x** for CAT-Poisson, **1.1x to 2.2x** for codon models. No configuration is slower than `master`.
- Speed-up at production size: 38 taxa x 59,994 amino-acid sites, CAT-GTR, 64 ranks, pinned: **16.35 s per cycle** against **71.28 s** for `master` (4.36x).
- **Pinning each rank to its own physical core is worth about 1.5x by itself**, more than any compiler flag tested. It changes no result. See `pinning-and-cpu-affinity.md`.
- The optimizations do not change the arithmetic. They remove repeated work, avoid allocation, cache values that were recomputed, and let the compiler use vector registers.
- All numbers come from one cluster (AMD EPYC 7713), one compiler (GCC 8.5.0) and one MPI (Open MPI 4.1.2). The model table has two pinned runs per cell; the history and scaling tables are single runs.

## 2. What was done, in one paragraph per area

1. **Less repeated work in the pruning (likelihood) code.** The slaves ran two full likelihood passes where one was enough. Several moves recomputed the whole tree after changing one branch. Conditional likelihoods are now recomputed only where needed.
2. **Faster matrix-vector kernel.** The kernel that propagates conditional likelihoods through a branch (`MatrixSubstitutionProcess::PropagateImpl`) is the largest single cost. It now uses contiguous transposed matrices, shared exponentials, a leaf shortcut and a register-blocked vector loop. Each output element still receives its terms in the original order.
3. **Cached logarithms in the mixture moves.** Allocation and Dirichlet-weight moves recomputed `log` and `lgamma` values that depend only on the component. They are now computed once per component.
4. **Memory layout.** Each conditional-likelihood vector is one contiguous block, and a small scratch buffer replaces a large allocation per call.
5. **Optional diagnostics.** `-DPB_MPI_STATS` records per-message execution and wait time on every slave.
6. **Tooling.** Slurm scripts for benchmarking, fixed-seed regression against `master`, a regression for `readpb_mpi`, a pinned-install script and a rank-scaling benchmark.

`code-optimizations.md` describes each change in detail.

## 3. Results

### 3.1 All models, `perf-opt` (4d7ef1e) against `master`, pinned

Conditions: 8 MPI ranks, each pinned to its own physical core (`--hint=nomultithread`, `mpirun --bind-to core --map-by core`), default compiler flags (`-Wall -O3 -std=c++11`, no AVX2), simulated alignments (40 taxa x 5000 amino-acid sites; 20 codon taxa x 500 codons for the `mutsel` models), `-rnd 1`. Reference: `master` at 3c77656. Two complete passes on 2026-11-20 (folders `bench/pin1_*` and `bench/pin2_*`), at most 6 jobs at a time on a node. Wall time includes MPI start-up and program set-up.

| Model flags | Cycles | Chain identical | master (s), run 1 / run 2 | optimized (s), run 1 / run 2 | Ratio, run 1 / run 2 |
|---|---|---|---|---|---|
| `-cat -gtr` | 20 | yes | 515.7 / 549.8 | 122.4 / 110.4 | 4.21x / 4.98x |
| `-cat -gtr -spr 2 -nni 2` | 10 | yes | 144.3 / 188.6 | 40.7 / 47.7 | 3.55x / 3.95x |
| `-ncat 10 -gtr` | 15 | yes | 240.1 / 256.9 | 47.8 / 52.8 | 5.02x / 4.87x |
| `-ncat 1 -lg` | 15 | yes | 245.9 / 228.2 | 46.5 / 41.3 | 5.29x / 5.53x |
| `-cat -poisson` | 60 | yes | 287.4 / 250.0 | 124.9 / 117.5 | 2.30x / 2.13x |
| `-cat -poisson -dc` | 30 | yes | 156.6 / 132.2 | 45.6 / 52.9 | 3.43x / 2.50x |
| `-ncat 10 -poisson` | 60 | yes | 219.0 / 289.4 | 111.1 / 106.2 | 1.97x / 2.73x |
| `-cat -poisson -sis 1000 2 3` | 1000 | yes, apart from the two corrected columns | 157.0 / 161.2 | 89.2 / 101.5 | 1.76x / 1.59x |
| `-cat -gtr -sis 1000 2 3` (20 taxa x 1000 sites; reference is fixed `master`) | 1000 | yes, apart from the two corrected columns | 164.9 / 187.3 | 160.7 / 169.8 | 1.03x / 1.10x |
| `-mutselc -cat` | 8 | yes | 219.7 / 203.4 | 123.4 / 93.6 | 1.78x / 2.17x |
| `-mutselc -ncat 3` | 8 | yes | 44.2 / 31.6 | 35.4 / 24.6 | 1.25x / 1.28x |
| `-mutsel -cat` | 8 | yes | 59.9 / 55.1 | 49.8 / 49.1 | 1.20x / 1.12x |
| `-mutsel -cat -freeomega` | 8 | yes | 61.5 / 50.2 | 48.1 / 38.1 | 1.28x / 1.32x |
| `-mutsel -ncat 3` | 8 | yes | 44.3 / 32.2 | 37.3 / 24.9 | 1.19x / 1.29x |

Read the table with care:

- Two runs per cell. Timings of one cell differ by up to about 30% between the two passes, because other jobs share the node.
- The speed-up is largest for CAT-GTR, the finite GTR mixture and LG (4x to 5.5x), 2x to 3.4x for the Poisson models, and 1.1x to 2.2x for the codon (`mutsel`) models. The likely reason is that the optimized code paths (the amino-acid pruning kernel and the mixture moves) are a smaller share of the time in the codon models. The codon models were not profiled, so this is an inference.
- Stepping-stone with CAT-GTR shows no gain: the test is one stepping point and is almost entirely per-site stepping code, which was not optimized.
- An earlier unpinned table (8 logical CPUs on 4 cores) gave similar ratios with wall times about 1.4x longer for both builds.

### 3.2 Production-size data set, pinned

Data: 38 taxa, 59,994 amino-acid sites; model CAT-GTR; a chain resumed from one saved state with the same seed; AVX2 builds; 64 ranks, pinned to cores.

| Build | Cycles run | Seconds per cycle (mean) |
|---|---|---|
| `master` (3c77656) | 8 | 71.28 |
| `master` without AVX2 flags | 8 | 69.67 (median 63.67) |
| 312c2bd | 20 | 16.35 (median 16.39) |

Ratio 71.28 / 16.35 = 4.36x. The first 8 cycles of `.trace` are identical between `master` and 312c2bd (see `validation.md`, section 3.2). One run per row.

### 3.3 History: where the speed came from

Cumulative, one run each unless stated, CAT-GTR, 40x5000, 8 ranks, unpinned, default flags, 20 cycles, wall seconds. Each row includes all earlier commits.

| Run | Commit | What the commit adds | Wall (s) |
|---|---|---|---|
| base_a, base_b | 3c77656 (`master`) | none | 954, 900 |
| opt1 | d2456fc | Unfold: one likelihood pass instead of two | 648 |
| opt2 | 9cb60c2 | Propagate: small scratch, shared exponentials, leaf shortcut | 683 |
| opt3 | f56c1f7 | cache `log(profile)` per component | 533 |
| opt4 | 288cc6e | drop redundant full updates in GibbsSPR (276960e) | 509 |
| opt5 | 14297a9 | column-wise matrix-vector products, hoisted virtual calls | 335 |
| opt6 | 81ad2c2 | one contiguous block per conditional-likelihood vector | 344 |
| opt7 | a183303 | partial post-order updates in GibbsSPR | 304 |
| `cmp_base0` (3 runs) | 31a9b87 | adds e62498b (cached logs in the weight move) | 309, 255, 233 |
| `cmp_mine0` (3 runs) | 312c2bd | register-blocked Propagate kernel, leaf cache, hash-free exp cache | 179, 170, 162 |

Notes:

- These are cumulative times. **No single change was measured on its own.** The commit that adds the most visible speed in this series is 14297a9 (508 s to 335 s), but one run per row cannot separate the effect of a change from run-to-run variation (6% to 21% for identical builds). Changes of less than about 20% between adjacent rows are not meaningful.
- opt2 is slower than opt1 (683 s against 648 s) although it removes work. This is within the noise, and it shows how coarse the series is.
- The last two rows come from a different batch and different node load from the first rows.

### 3.4 Where the time goes now

A `perf` profile of a slave rank on the 40x5000 CAT-GTR run at 312c2bd gave these shares of self time (not re-measured for this document):

| Function | Share |
|---|---|
| `MatrixSubstitutionProcess::PropagateImpl` | 31% |
| MPI polling (idle waiting) | 9% |
| `SubstitutionProcess::Multiply` | 8% |
| QR and diagonalization (`LinAlg`) | 7% |
| `LogStatProb` | 4% |

The master spends about 91% of its time in MPI polling (from an earlier profile). It is mostly idle while the slaves work. A profile of an earlier commit (81ad2c2, unpinned) showed `PropagateImpl` at 48%, so the kernel work reduced its share.

### 3.5 Scaling with the number of ranks

Data: 38 taxa x 59,994 sites, CAT-GTR, pinned, AVX2 build of 4d7ef1e (`builds/4d7ef1e-avx2`). Each row is one run of 20 cycles, resumed from the same saved state (2026-11-20). Speed-up and efficiency are relative to 16 ranks. One rank is the master and does little work. The last column is the earlier measurement with build 31a9b87, which lacks the kernel of 312c2bd.

| Ranks | Seconds per cycle | Speed-up vs 16 | Efficiency | Core-seconds per cycle | Earlier build 31a9b87 (s per cycle) |
|---|---|---|---|---|---|
| 16 | 41.19 | 1.00x | 100% | 659 | 53.96 |
| 24 | 29.17 | 1.41x | 94% | 700 | 38.15 |
| 32 | 24.28 | 1.70x | 85% | 777 | 32.92 |
| 48 | 20.14 | 2.05x | 68% | 967 | 25.35 |
| 64 | 18.01 | 2.29x | 57% | 1153 | 20.75 |
| 96 | not re-measured | | | | 14.83 |
| 128 | not re-measured | | | | 14.15 |

What this means:

- The current build is 1.15x to 1.36x faster than 31a9b87 at every rank count measured, and the gain is largest at low rank counts.
- Efficiency falls faster than before: 57% at 64 ranks against 65% with the earlier build. The kernel work removed computation, and the communication and waiting time did not shrink with it.
- A cycle costs about 1.75x more core time at 64 ranks than at 16 ranks.
- Fewer ranks is the better choice for total throughput. More ranks is the better choice for the wall time of one chain.
- 96 and 128 ranks were not re-measured, because they did not fit under the group CPU limit while production chains were running. With the earlier build, going from 96 to 128 ranks reduced the time per cycle by only 5%.
- Each row is a single run, and the rows ran on different nodes (16 to 48 ranks on r07, 64 ranks on r10). An earlier 64-rank run of 312c2bd on r42 gave 16.35 s per cycle, so node-to-node differences of about 10% are present.
- The 20 cycles of each run have the same `.trace` (apart from the timing columns) as the run of the earlier build at the same rank count.

### 3.6 What did not help (not merged)

| Idea | Result | Where |
|---|---|---|
| Send the elementary branch-length commands in batches (182d90b) | chains identical, **no speed-up** at 8 and 64 ranks; reverted in 7683baa | branch `mpi-batching` |
| Ask for the number of states once per site in the pruning helpers (5941ab9) | about 8% on shared cores, **nothing measurable on pinned cores** (95.0 s against 95.9 s), validated only for CAT-GTR | branch `hoist-nstate` |
| A first register-accumulator prototype (22ccd4d) | 1.2x over the previous code, superseded by the kernel of 312c2bd | not kept |
| `-mtune=znver1`, `-mno-avx256-split-unaligned-*` | no gain separable from noise (unpinned) | `build-and-compiler-options.md` |

## 4. Reliability of the numbers

- **Run-to-run spread.** Identical builds on identical settings differ by 6% to 21% (for example 954 s and 900 s for `master`; 288 s and 262 s for one optimized build). Single-run differences below that range are not evidence.
- **Pinned and unpinned runs are not comparable.** Unpinned 8-rank runs use 4 physical cores and are about 1.5x slower. The early history in section 3.3 is unpinned; the tables in sections 3.1, 3.2 and 3.5 are pinned. Use the unpinned numbers to see the size of the effect, not to rank small variants.
- **One data set family.** The timing data come from simulated alignments of 40 taxa x 5000 sites, one real data set of 38 taxa x 59,994 sites, and short codon runs. Other shapes (many taxa, few sites, very large alignments) were not measured.
- **Ratios against `master` depend on `master` timings that vary** from batch to batch.
- **Not measured:** peak memory (the optimized code adds small per-matrix caches and a scratch buffer, and stores each conditional-likelihood vector in one block; the change in memory was not measured), start-up cost, and behavior with a different MPI library.

## 5. How to reproduce

All scripts run from the repository root on a Slurm cluster. Output goes to `bench/`, which is git-ignored.

```
# 1. Build and time a fixed-seed chain. Pinned, whole cores:
MPIRUN_OPTS="--bind-to core --map-by core" \
  sbatch --hint=nomultithread scripts/bench/bench.sbatch <ref> <id> 20

# 2. Compare two runs: chain identity and speed-up
scripts/bench/compare_runs.sh bench/<idA> bench/<idB>

# 3. Regression against master for 14 model configurations
PREFIX=<name> scripts/bench/regression.sh submit <ref> master
PREFIX=<name> scripts/bench/regression.sh check

# 4. Regression for readpb_mpi
sbatch scripts/bench/readpb_regression.sbatch <ref> master

# 5. Install a pinned build with a record of how it was made
sbatch scripts/install_build.sh <ref>

# 6. Rank scaling from a saved chain state
sbatch -n <ranks> --mem <mem> scripts/bench/scaling.sbatch <build dir> <id> <ncycles> <.param> <alignment>
```

`scripts/bench/make_testdata.sbatch` makes the simulated alignments. Edit the `-C milan` and partition lines in the scripts for another cluster.

## 6. Open items

1. Re-measure the rank-scaling table at 96 and 128 ranks with the current build (16 to 64 ranks are done).
2. (Done: `-mutsel -ncat 3` is not slower than `master`, section 3.1.)
3. Pinned test of the compiler flags in `build-and-compiler-options.md`.
4. `pb_mpi -cat -gtr -sis` aborted with a NaN in `master` itself. This is fixed (uninitialized empirical-prior arrays of the relative rates; see `validation.md`, section 3.4). The `logZ` and `effsize` columns of `.stepping` printed `-nan` for large steps, in `master` as well; this separate numerical problem is fixed in db61cc5 (`validation.md`, section 3.4). One open item: stepping-stone CAT-GTR gave a second result in 4 of 52 otherwise identical runs, not reproduced since (`validation.md`, section 3.5).
5. Remaining hot spots: `Multiply`, QR and diagonalization, `LogStatProb`, and idle time waiting for the master.

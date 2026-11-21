# Validation: do the optimized programs give the same results as the original?

> Part of the performance documentation; see the scope note at the top of `PERFORMANCE.md`.

This note states the rule, shows how it is tested, lists the results, and lists what is **not** tested.

## 1. The rule

Every optimization must give the same result as `master` for the same data, model, seed (`-rnd`) and number of MPI ranks. "The same result" means the same Markov chain:

- the same sampled trees (`.treelist`),
- the same saved parameter states (`.chain`),
- the same `.trace`, except the timing columns,
- the same stepping-stone output (`.stepping`) where it applies.

Small rounding differences would be allowed only if they were identified, explained and shown to be rounding-level. **No such difference was found.** All comparisons below are exactly identical, byte for byte.

Anything that changes the arithmetic is a separate class. It must be shown as an explicit accuracy-versus-speed trade-off and must be off by default. This branch contains no such change.

## 2. How a comparison works

`scripts/bench/bench.sbatch` builds a git ref and runs a fixed-seed chain (`-rnd 1`). `scripts/bench/compare_runs.sh <dirA> <dirB>` then compares two such runs:

| File | Comparison |
|---|---|
| `chain.treelist` | `cmp`, whole file |
| `chain.chain` | `cmp`, whole file (every saved parameter state, written with 15 significant digits) |
| `chain.stepping` | `cmp`, whole file (when present) |
| `chain.trace` | `diff` of all columns except column 2 (`time`) and column 3 (`topo`) |

Column 3 of the trace is the percentage of time spent in one timer, not a topology count (see `Trace()` in `sources/RASCATGTRSBDPGammaPhyloProcess.h`). Both ignored columns are timing columns. The tree itself is covered by `.treelist`.

The script exits with status 1 if any file differs or nothing could be compared.

`scripts/bench/regression.sh` runs 14 model configurations with a `master` build and a new build, then calls `compare_runs.sh` for each pair:

```
PREFIX=<name> scripts/bench/regression.sh submit <new-ref> master   # submits 28 jobs
PREFIX=<name> scripts/bench/regression.sh check                     # after the jobs finish
```

`PREFIX` names the output folders. Use a different prefix for each session so that two runs do not overwrite each other.

`scripts/bench/readpb_regression.sbatch <new-ref> master [model]` makes a short chain with the old build, then runs each `readpb_mpi` analysis with both builds on copies of that chain and compares every output file. `readpb_mpi` has no seed option and seeds from the clock. The script therefore preloads a small shim (`fixusec.so`) that fixes the microsecond field of `gettimeofday`, so that stochastic analyses can be compared.

## 3. Results for `perf-opt` at 4d7ef1e

The speed changes end at commit 312c2bd. The later commits add scripts, documentation and the two corrections of pull request 39. Sections 3.1 to 3.3 were measured on 312c2bd or earlier and name their commit.

Test conditions: 8 MPI ranks, each pinned to its own physical core (`--hint=nomultithread`, `mpirun --bind-to core --map-by core`), default compiler flags (`-Wall -O3 -std=c++11`, no AVX2), simulated alignments (40 taxa x 5000 amino-acid sites; 20 codon taxa x 500 codons for the `mutsel` models), `-rnd 1`. Reference: `master` at 3c77656. Two complete passes on 2026-11-20 (folders `bench/pin1_*` and `bench/pin2_*`), at most 6 jobs at a time on a node.

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

Remarks:

- **All 12 configurations without `-sis`** give identical `.treelist`, `.chain` and `.trace` (apart from the timing columns) in both passes.
- **The two stepping-stone configurations** give a `.stepping` file that is identical to the reference in every column except `logZ` and `effsize`. The reference prints `-nan` there, and the optimized branch prints finite values because of the correction db61cc5 (section 3.4). `compare_runs.sh` reports this case separately.
- **Stepping-stone with CAT-GTR.** `master` itself aborts with `in Plink::SetRelativeTime: setting to nan`. The reference for this row is therefore `master` plus the 8-line correction 9f2a6fc (branch `fix-sis-gtr-nan`). This test is almost entirely per-site stepping code, which the optimizations do not touch, so there is no gain. See also section 3.5.
- **`-mutsel -ncat 3`** is not slower than `master`. An earlier single unpinned run suggested 0.83x; it was not confirmed here, nor in 3 pinned runs of 30 cycles (`master` 147, 171, 154 s; 312c2bd 134, 137, 103 s).
- Timings still vary between the two passes by up to about 30% for one cell, because other jobs share the node. The ratios are good to a few tenths.

### 3.1 `readpb_mpi`

`readpb_regression.sbatch` was run on 2026-11-20 on the final sources (312c2bd plus the two stepping corrections and the four small fixes of `code-optimizations.md`, section 12) against `master`, on the small alignment `sim20x1000`. All 30 analyses gave identical output files:

| Model | Analyses | Result |
|---|---|---|
| `-cat -gtr` | 14: default, `-sitelogl`, `-anc`, `-r`, `-ss`, `-ppred`, `-allppred`, `-comp`, `-div`, `-map`, `-jointcv`, `-sitecv`, `-rr`, `-siteprofilesuffstat` | identical |
| `-cat -poisson` | 12: the same list without `-rr` and `-siteprofilesuffstat` | identical |
| `-ncat 10 -gtr` | 4: default, `-sitelogl`, `-ppred`, `-r` | identical |

An earlier run of the same job on a candidate before 312c2bd (commit 22ccd4d) had given the same result. The job was not run on 312c2bd itself.

### 3.2 A production-size check

A chain from a real data set (38 taxa, 59,994 amino-acid sites, CAT-GTR) was resumed from one saved state for several cycles with `master` and with 312c2bd, both built with `-mavx2`, 64 ranks, pinned. The first 8 cycles of `.trace` (all columns except the two timing columns) are identical. Master ran 8 cycles; 312c2bd ran 20. Only the first 8 were compared. `.treelist` of the 8 cycles was not compared in this check.

### 3.3 Compiler flags do not change the chain

`master` built with default flags and `master` built with `-mavx2` (64 ranks, same saved state, 8 cycles): `.trace` (all columns except the two timing columns) and `.treelist` are identical.

For 312c2bd (CAT-GTR, 40x5000, 8 ranks, 20 cycles, compared with `compare_runs.sh`):

- default flags against `-mavx2`: `.treelist`, `.chain` and `.trace` identical;
- `master` with default flags, unpinned, against 312c2bd with `-mavx2`, pinned: all three identical. This single comparison changes the compiler flag, the placement of the ranks and the code at the same time.

For earlier candidates (commit 51fdf5f and the later-superseded commit 22ccd4d), runs with `-mavx2 -mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store` and with `-mavx2 -mtune=znver1` also gave chains identical to the default-flag chain. Those flag variants were not re-run on 312c2bd.

### 3.4 Stepping-stone with CAT-GTR (`-sis`): the abort in `master` and its fix

**Cause.** `GTRProfileProcess::Create` allocated `emprralpha` and `emprrbeta` (the empirical-prior parameters of the relative rates) and never gave them values. In the runs examined they were all zero. They are read when `profilefrac != 1.0`, as `a = f + (1-f)*emprralpha` and `b = f + (1-f)*emprrbeta`. `SteppingRun` calls `GlobalSetEmpiricalFrac(frac)` at every stepping point, even when no empirical prior was loaded, and `frac` is 0 for the first point. `profilefrac` then becomes 0, `a` and `b` become 0, and `Gamma(a + count, b + beta)` gives NaN relative rates. The substitution matrix becomes NaN, path sampling stores NaN in `Plink::rel_time`, and the next write to that object aborts, because `Plink::SetRelativeTime` tests its old value. Models without relative rates to estimate (CAT-Poisson) were not affected, because their empirical-prior arrays default to 1.0.

**Evidence** (instrumented scratch build of `master`, 20 taxa x 1000 sites, 8 ranks): all master-side parameters were finite before every broadcast to the slaves, up to the call where `profilefrac` first printed 0. From that call on, all 190 relative rates were NaN, and the first NaN stationary-matrix row followed immediately. The `emprralpha` and `emprrbeta` values printed as 0.

**Fix.** Initialize both arrays to 1.0 (Gamma(1,1), the ordinary prior), as the Dirichlet weights, kappa and branch lengths already do. 8 lines in `sources/GTRProfileProcess.cpp`.

**Validation.**

- The failing case exits 0 with the fix and exited 1 without it.
- Chains that do not use `-sis` never read these arrays (`profilefrac` is 1.0). Fixed `master` against unfixed `master`: `catgtr` (20 cycles) and `ncat10gtr` (15 cycles) identical in `.treelist`, `.chain` and `.trace`; `stepping_catpoisson` identical in `.stepping`.
- Fixed `perf-opt` against fixed `master`: `catgtr`, `ncat10gtr`, `stepping_catpoisson` and `stepping_catgtr` all identical (`.stepping` for the stepping models). The `stepping_catgtr` run is one stepping point, 0 to 1000 sites, on `sim20x1000`.

**The `logZ` and `effsize` columns (fixed in db61cc5).** These two columns of `.stepping` printed `-nan` when the first step had a very negative log-likelihood, in `master` too (also for CAT-Poisson). In the accumulator of `SteppingRun` (`sources/Model.h`), the first rescaling was `0 * exp(+large)`. The fix skips the rescaling while the running sums are still empty. It does not change the chain, and it is not a speed change.

- Stepping CAT-Poisson with and without the fix, on `master` and on `perf-opt`: the `.stepping` files differ only in columns 3 (`logZ`) and 8 (`effsize`), which go from `-nan` to finite values.
- A model test of the accumulator with 20,000 random inputs: the old logic gave NaN or infinity in 38% of cases; the new logic was always finite, equal to the old wherever the old was finite, and within 1.3e-15 of a reference calculation.
- Before the fix, both sides printed `-nan` in the same columns, so an identical `.stepping` file showed that `meanlogp` and `varlogp` agreed but said nothing about `logZ`.

### 3.5 Open: stepping-stone CAT-GTR gave a second result in 4 of 52 runs (not reproduced)

The test is `-cat -gtr -sis 1000 2 3 -rnd 1` on `sim20x1000`, 8 ranks, one stepping point. With the same data, seed, rank count and compiler flags, the `meanlogp` column was -41136.2 in 48 runs and -36484.3 in 4 runs.

- The 4 runs were all on 2026-11-19 between 22:01 and 22:15, on nodes r41 and r42: one of `perf-opt` 1077165 and three of db61cc5. Other runs on the same nodes in the same hour gave -41136.2.
- It is not a build difference: the db61cc5 binary that gave -36484.3 three times gave -41136.2 in 7 of 7 runs the next day. It is not the `logZ` fix either, because both code versions show both values.
- Attempts to reproduce it on 2026-11-20 all gave -41136.2: 7 repeats of that binary (node r07), 8 runs through `bench.sbatch` (r30), and 8 runs fixed to r41 and r42 (4 of fixed `master`, 4 of `perf-opt`).
- Fixed `master` (branch `fix-sis-gtr-nan`) gave -41136.2 in all 8 of its runs, so whether `master` can show the second value is not known.
- No cause was found. Filling the heap with 0x00 or 0xFF, and `MALLOC_PERTURB_`, did not change the result. valgrind memcheck on the first 100 sites reported no uninitialised reads in `pb_mpi`.
- One suspect, not shown to matter: the master adds up the replies of the slaves in arrival order (`MPI_ANY_SOURCE` in `PhyloProcess::GlobalComputeNodeLikelihood` and similar loops), so a sum can depend on timing at rounding level.
- Stepping CAT-Poisson and all models without `-sis` have never shown run-to-run variation.

To test again later: repeat the command above many times and compare `chain.stepping`; if a second value appears, add per-cycle output to find the first cycle that differs. Until then, treat a single differing stepping CAT-GTR run as possibly this effect, and repeat it before concluding that a code change altered the chain.

### 3.6 The branches of the upstream pull requests

For the upstream repository the speed changes are split into three branches, each on top of `fix-stepping-stone` (`master` plus the two stepping corrections, pull request 39):

| Branch | Content (sections of `code-optimizations.md`) |
|---|---|
| `perf-pruning` | 1, 4, 7, 8: fewer conditional-likelihood updates in Unfold and GibbsSPR; optional message statistics |
| `perf-kernels` | 2, 5, 6, 10: the Propagate and pruning kernels, contiguous vectors |
| `perf-mixture` | 3, 9: cached logarithms in the allocation and weight moves |

Each branch, and the merge of the three (`perf-combined`), was built alone and compared with `fix-stepping-stone` on the 14 models of the table above (2026-11-20, 8 ranks, pinned). All four are identical to the reference on all 14 models, including all columns of `.stepping` for the two stepping-stone models. All five builds have no compiler warning.

- The three branches do not depend on each other. One point needs care when they are combined: `perf-kernels` adds `PropagateFromLeaf`, and the partial post-order update of `perf-pruning` can use it for leaves. `perf-pruning` alone uses plain `Propagate` there. `perf-combined` has one extra commit that switches to the leaf shortcut. The result is the same with or without it.
- The sources of `perf-combined` are the same as those of `perf-opt`.

Wall time in seconds of these runs (one run per cell; the `perf-pruning` column was run in a later batch than the others):

| Model | base | `perf-pruning` | `perf-kernels` | `perf-mixture` | `perf-combined` |
|---|---|---|---|---|---|
| `catgtr` | 525 | 308 | 250 | 296 | 87 |
| `catgtr_nni` | 185 | 110 | 78 | 152 | 38 |
| `ncat10gtr` | 169 | 137 | 54 | 178 | 38 |
| `lg` | 259 | 108 | 50 | 179 | 51 |
| `catpoisson` | 194 | 222 | 177 | 157 | 110 |
| `catpoisson_dc` | 148 | 130 | 116 | 90 | 54 |
| `ncat10poisson` | 153 | 118 | 96 | 214 | 108 |
| `stepping_catpoisson` | 155 | 169 | 151 | 155 | 108 |
| `stepping_catgtr` | 169 | 191 | 188 | 199 | 184 |
| `mutselc_cat` | 202 | 259 | 198 | 113 | 126 |
| `mutselc_ncat3` | 45 | 38 | 32 | 44 | 30 |
| `mutsel_cat` | 55 | 59 | 48 | 47 | 52 |
| `mutsel_freeomega` | 68 | 51 | 53 | 52 | 45 |
| `mutsel_ncat3` | 43 | 36 | 29 | 34 | 27 |

These single runs are a rough guide only. Other jobs share the nodes, and the same build varied by up to 45% between batches (see below). Use section 3 of `PERFORMANCE.md` for the speed of the complete work. What the table does show: the gains of the three branches are separate and add up (CAT-GTR: 525 s for the base, 250 to 308 s for each branch alone, 87 s for all three).

Three cells where a branch looked slower than the base were repeated 3 times, base and branch in the same batch:

| Model | Branch | base (s) | branch (s) |
|---|---|---|---|
| `catpoisson` | `perf-pruning` | 297, 265, 281 | 230, 229, 234 |
| `mutselc_cat` | `perf-pruning` | 253, 250, 253 | 241, 243, 242 |
| `ncat10poisson` | `perf-mixture` | 215, 216, 148 | 213, 213, 150 |

So `perf-pruning` is not slower on these models, and `perf-mixture` makes no difference for the finite Poisson mixture. The slower single cells above were noise.

## 4. What is not validated

Be explicit about these in any pull request:

1. **Different numbers of MPI ranks.** All comparisons use the same rank count for both builds. The chain depends on the rank count (see `pinning-and-cpu-affinity.md`), so a chain from 8 ranks is not expected to equal a chain from 64 ranks, for `master` or for the optimized code. This was not tested by running different rank counts.
2. **Stepping-stone with CAT-GTR** (`-sis`) is validated only for one stepping point on a 20-taxon, 1000-site alignment (section 3.4), and 4 of 52 runs of that test gave a different result for an unknown reason (section 3.5). Several stepping points, the 40x5000 alignment and the `-empstepping` option were not tested with CAT-GTR.
3. **Models and options outside the table.** The 14 configurations do not use, for example, a fixed topology (`-T`, `-Tbl`), a non-default `-dgam`, or `-catfix`. Restarting a stopped chain from its `.param` file (`pb_mpi <name>`) is covered only by the single production-size check in section 3.2. Any option not named in the table should be treated as untested.
4. **Large or long runs.** The longest identity comparison is 20 cycles of CAT-GTR on 40x5000 (about 8 saved cycles on the 38x59,994 data set). Rounding differences that appear only after many thousands of cycles would not be seen. The arithmetic is argued to be unchanged in `code-optimizations.md`, so none is expected.
5. **Other compilers, other CPUs.** All tests used GCC 8.5.0 and Open MPI 4.1.2 on AMD EPYC 7713 (Zen 3). A different compiler or CPU can change floating-point results through instruction contraction or different library code. Re-run the regression after any such change.
6. **Statistical validation.** Not needed while chains are identical, and not done.

## 5. Checklist for a new change

1. Build the candidate and `master`; run `PREFIX=<name> scripts/bench/regression.sh submit <ref> master`, then `check`. All models must report identical.
2. Run `scripts/bench/readpb_regression.sbatch <ref> master`.
3. Add a row for any new model or option you touched.
4. Time with pinned runs and at least 3 repeats (see `pinning-and-cpu-affinity.md`).
5. If a result is not identical, stop. Find the first differing cycle, explain it, and classify the change under rule 2 of section 1.

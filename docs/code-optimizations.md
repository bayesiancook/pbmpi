# Code optimizations: what changed and why the results are unchanged

> Part of the performance documentation; see the scope note at the top of `PERFORMANCE.md`.

This note describes every code change on `perf-opt` that affects speed, in the order of the commits. Compared are `master` (3c77656) and `perf-opt` at 312c2bd, the last commit that changes speed. The later commits are scripts, documentation and two corrections of errors in `master` (9f2a6fc and db61cc5, see pull request 39).

How to read the evidence tags:

- **[code]**: read in the diff or in the source at the cited line.
- **[bench]**: read in a benchmark output file.
- **[inference]**: follows from the code, but no run proves it.
- **[rechecked]**: the item was checked a second time, independently of the first reading.

Line numbers refer to the sources at HEAD, unless marked `master:`.

## 0. Principle

Every change obeys one rule. **Same arithmetic, same order of operations, same sequence of random numbers.** The changes remove work, they do not approximate it. Six kinds of change are used:

1. Skip a computation whose result is thrown away or recomputed without change (sections 1, 4, 7).
2. Compute a value once and reuse it when the inputs are bit-identical (sections 2, 3, 9, 10).
3. Replace many small allocations by one large allocation (sections 2, 6).
4. Change the memory layout so that the compiler can use vector instructions, while each output keeps its own sequence of operations (sections 5, 10).
5. Avoid repeated virtual calls in loop conditions (section 5).
6. Diagnostics that are off by default (section 8).

Two things are **not** done: no sum is reordered, and no multiply-add is fused. If a future change needs either, it is an accuracy-versus-speed trade-off and must be opt-in (see `validation.md`).

## 1. Unfold computes the conditional likelihoods once (d2456fc)

**What was slow.** `GlobalUnfold` sent the message `UNFOLD`. The slave ran a complete two-pass update of the conditional likelihoods inside `Unfold()`. The master then sent `UPDATE`, which ran the same update again (master:PhyloProcess.cpp:1306-1315, master:MatrixPhyloProcess.cpp:50). The CAT models do this 5 times per cycle [code].

**What changed.**

- Two new flags with default initializers: `unfoldupdate = true` on the slave and `postorderunfold = false` on the master (PhyloProcess.h:623, 625) [code].
- In the slave handler for `UNFOLD`, `unfoldupdate` is false around the call (PhyloProcess.cpp:1821-1827). The three `Unfold()` bodies skip their update when it is false (PhyloProcess.cpp:44, MatrixPhyloProcess.cpp:50, GeneralPathSuffStatMatrixPhyloProcess.cpp:48) [code].
- A new message, `UPDATE_POSTORDER`, makes the slave run `PostOrderUpdateConditionalLikelihoods()` (PhyloProcess.cpp:241-245, 1838). That runs the post-order pass, `MultiplyByStationaries` and `ComputeLikelihood`, but no pre-order pass [code].
- In the three CAT model loops (RASCATGTRSBDPGammaPhyloProcess.h:323, RASCATSBDPGammaPhyloProcess.h:286, RASCATFiniteGammaPhyloProcess.h:307), the unfolds for `rep < 4` use `GlobalUnfoldBeforeCollapse()`. The next step after them is a collapse, and the collapse reads only post-order vectors [code].

**Why the result is unchanged.**

- A full update reads only the data, the branch lengths, the parameters and vectors that it writes itself in the same pass (PhyloProcess.cpp:170-178, 333-386). The second update on `master` recomputed the same values [code, inference].
- `Collapse` draws allocations and samples node states. `SampleNodeStates` reads the post-order vectors. It reads a pre-order vector only after the parent has written it (PhyloProcess.cpp:433-459). Stale pre-order values are therefore overwritten before use [code].
- After the post-order update, the master sends the same likelihood request. The root likelihood reads only post-order vectors (PhyloProcess.cpp:290-331) [code].
- No random-number call is added, removed or moved [code].

**What a reviewer should check.**

- The shortcut is valid only because a collapse follows directly. A future change that inserts another move between the unfold and the collapse would read pre-order vectors that were not updated. Add a comment or an assertion if you extend these loops.
- The skipped vectors are initialized (1.0, offset 0; SubstitutionProcess.cpp:112-116), not left undefined [code].
- Direct slave calls to `Unfold()` (PoissonPhyloProcess.cpp:348 and the gene models) keep `unfoldupdate == true` and behave as on `master` [code].
- `PostOrderUpdateConditionalLikelihoods` is not virtual. `MatrixPhyloProcess::UpdateConditionalLikelihoods` has the same body today. A future override with a different body would be bypassed [code, inference].

**Portability.** Default member initializers need C++11. The Makefile already uses `-std=c++11`.

**Measured.** Not measured separately. Cumulative run opt1 (this commit and all earlier ones): 647.9 s, against 954.0 s and 900.2 s for `master` [bench].

## 2. Matrix `Propagate`: scratch buffer, shared exponentials, leaf shortcut (9cb60c2; final form in 312c2bd)

**What was slow** (master:Propagate.cpp:82-232) [code]:

- Each call allocated `new double[Nsite*Nrate*nstate]` and freed it again (master:89, 231).
- It computed `exp(length*eigenval[k])` for every site, rate and state (master:148).
- For a leaf with an observed state, it ran the full product `P^-1 . up`, although `up` is a unit vector.

**What changed.**

- `Propagate` and a new virtual `PropagateFromLeaf` (SubstitutionProcess.h:143; the default calls `Propagate`) both call `PropagateImpl(from, to, time, condalloc, leafstate)` (Propagate.cpp:82-88, 203) [code].
- `PostOrderPruning` calls `PropagateFromLeaf` when the child is a leaf (PhyloProcess.cpp:341-347) [code].
- Scratch members are in MatrixSubstitutionProcess.h:50-68: `propaux`, `propentries`, `propheads`, `expbuffer`, `leafindex`, `leafbuffer`, `leafclipped`, `propstamp`. They are cleared at the start of each call (Propagate.cpp:211-217). No cache survives from one call to the next [code].

**Why the result is unchanged.**

- Exponentials: `length = time*rate` and `ex[k] = exp(length*eigenval[k])` (Propagate.cpp:247-251) are the same expressions as on `master`. A cached value is reused only for the same matrix object and an identical `rate` value (Propagate.cpp:237) [code].
- Leaf shortcut: on `master` the sum for `aux[k]` runs over all `l` with `up` equal to a unit vector. Every term except `l = obs` is a product with zero. The new code computes `(0.0 + col[k]) * expdiag[k]` (Propagate.cpp:275). The `0.0 +` keeps the sign of zero identical to the full sum [code, inference].
- Missing data (state `-1`) takes the general path (Propagate.cpp:226, 262) [code].

**What a reviewer should check.**

- The order of the fatal-error checks changed. The new code tests `up[k] < 0` before the products (Propagate.cpp:318-328). `master` tested for NaN first (master:176-195). The leaf path has no `up < 0` check. Only the choice of the error message can differ [code].
- If the inverse eigenvector matrix contained `inf` or `NaN` in a column other than the observed state, `master` would produce NaN and exit. The leaf path would not see it [inference].
- `PropagateFromLeaf` must only be called right after `Initialize` with the same data. The two call sites do this (PhyloProcess.cpp:224, 343) [code].

**Measured.** Not measured separately. Cumulative run opt2: 682.7 s [bench].

## 3. Allocation moves: `log(profile)` cached per component (f56c1f7)

**What was slow.** `LogStatProb(site, cat)` computed `log(profile[cat][k])` for every site, component and state (master:PoissonMixtureProfileProcess.cpp:109, master:ExpoConjugateGTRMixtureProfileProcess.cpp:119). It is called inside the loops over sites and components in `SlaveMixMove` and `SlaveIncrementalFiniteMove` [code].

**What changed.**

- `MixtureProfileProcess.h:186-206` adds `logprofilecache`, `logprofilecachesize`, `PrepareLogProfileCache(ncat)`, `ReleaseLogProfileCache()` and `GetCachedLogProfile(cat)` (it returns 0 when `cat` is outside the cache) [code].
- Prepare and Release surround the site loops (PoissonSBDPProfileProcess.cpp:176/241, MatrixSBDPProfileProcess.cpp:261/325, PoissonFiniteProfileProcess.cpp:144/185, MatrixFiniteProfileProcess.cpp:94/135) [code].
- The two `LogStatProb` bodies use the cache when it exists and run the old loop otherwise (PoissonMixtureProfileProcess.cpp:108-119, ExpoConjugateGTRMixtureProfileProcess.cpp:118-128) [code].

**Why the result is unchanged.** `log(x)` of the same double gives the same double. The products and the order of additions over the states are unchanged. Components outside the cache fall back to the direct `log`. The random-number calls in the loops are untouched [code].

**What a reviewer should check.**

- Invalidation is manual. The cache is wrong if a profile changes between Prepare and Release. No profile write and no early `return` were found inside the bracketed loops [code].
- Release sets the size to 0 before the `MPI_Send` and before the profile move [code].
- `GeneralPathSuffStatMatrixMixtureProfileProcess::LogStatProb` (.cpp:123) does not use the cache. For those models Prepare adds work (`K0 * dim` logarithms per repetition) and gives no gain [code, inference].
- Prepare reads `profile[cat][k]` for every `cat < ncat`. This assumes that every component has `GetDim()` entries. Not verified for every model [inference].

**Measured.** Not measured separately. Cumulative run opt3: 533.0 s [bench].

## 4. GibbsSPR: drop the redundant full updates around the series (276960e)

**What was slow.** `GibbsSPR(int nrep)` called a full update of all conditional likelihoods before the loop and again after it (master:PhyloProcess.cpp:575-609). Every inner `GibbsSPR()` already updated after the re-rooting (master:1521) and at its end (master:697) [code].

**What changed** (PhyloProcess.cpp:650-666). Both outer updates are removed. For `nrep <= 0`, one update is done and 0 is returned. Dead commented code is removed [code, rechecked].

**Why the result is unchanged.** A full update does not depend on earlier vector contents (see section 1). The first inner update follows the re-rooting. The removed trailing update repeated the last inner one. The master-side random-number order is unchanged. The slaves receive fewer messages, in the same relative order [code, inference].

**Behavior difference.** For `nrep == 0`, `master` returned `0/0` (NaN) after two updates. The new code returns 0 after one update. `MoveTopo` adds this value to its success count (PhyloProcess.cpp:642). The change is not in the commit message. Calls with `nrep == 0` are not expected in normal runs [code, rechecked].

**Measured.** Not measured separately. Cumulative run opt4: 508.7 s [bench].

## 5. Pruning kernels: column-wise products and hoisted virtual calls (14297a9)

**What was slow.** `master` read `inveigenvect[k][l]` and `eigenvect[k][l]` row by row through `double**` (master:Propagate.cpp:129-170). The helpers in master:SubstitutionProcess.cpp called the virtual functions `GetNrate(i)` and `GetNstate(i)` in loop conditions [code].

**What changed.**

- `SubMatrix` stores a transposed, contiguous copy of each eigenvector matrix: `ut` and `invut`, with `ut[l*Nstate+k] = u[k][l]` (SubMatrix.h:194-195, allocated in SubMatrix.cpp:72-73, freed in 117-118). They are filled at the end of `Diagonalise()`, just before `diagflag = true` (SubMatrix.cpp:305-312). The accessors `GetTransposedEigenVect` and `GetTransposedInvEigenVect` (SubMatrix.h:120-127) diagonalize on demand [code].
- `Diagonalise()` is the only place that sets `diagflag = true`, so the copies cannot be stale relative to `u` and `invu` [code].
- In SubstitutionProcess.cpp about 40 loops use `for (int j=0, nrate=GetNrate(i); ...)` and `for (int k=0, nst=GetNstate(i); ...)` [code].
- `ARCHFLAGS` was added to the Makefile (empty by default). See `build-and-compiler-options.md`.

**Why the result is unchanged.** For each output `k`, the terms are added in the order `l = 0..n-1`, starting from `0.0`, with a separate multiplication and addition (the scalar `MatVecCols`, Propagate.cpp:188-199). This is the per-element order of the loops in `master`. The hoisting only caches loop bounds [code].

**What a reviewer should check.**

- Contraction into fused multiply-add changes the result. With `-mfma` or `-march=native`, or on AArch64 where FMA is part of the base instruction set, `master` and `perf-opt` could contract `a*b+c` differently. Identity is not established for such builds. `-mavx2` has no FMA, and the AVX2 builds gave identical chains [inference, bench].
- Each diagonalization makes two extra copies of an `n x n` matrix [code].
- `__restrict__` is a GCC and clang spelling. MSVC does not accept it. `master` already uses variable-length arrays (master:SubstitutionProcess.cpp:434, 671; Tree.cpp:912), which MSVC also rejects [code].

**Measured.** Not measured separately. Cumulative run opt5: 335.0 s (302.6 s with `-mavx2`) [bench].

## 6. One contiguous block per conditional-likelihood vector (81ad2c2)

**What was slow.** One `new double*[]` per site and one `new double[nstate+1]` per site and rate, for each link, at each Unfold. The matching deletes ran at each Collapse (master:SubstitutionProcess.cpp:85-116) [code].

**What changed** (SubstitutionProcess.cpp:90-130). One `double**` block (sum of the rate counts) and one `double` block (sum over sites of `Nrate*(Nstate+1)`) per vector. `condl[i]` is first set to 0 for every site. The delete frees `condl[sitemin][0]`, `condl[sitemin]` and `condl`. The access `condl[i][j][k]` is unchanged [code].

**Why the result is unchanged.** Only addresses change. Initial values are the same (1.0, offset 0) [code].

**What a reviewer should check.**

- The delete needs the same `sitemin` and `sitemax` as the create, and `GetNrate(sitemin) >= 1` [code].
- Code that frees or replaces a single `condl[i][j]` would now break. None was found. `SaveCondl` (NNI.cpp:359-370) and `sitecondlmap` (Stepping.cpp:152-173) allocate their own arrays [code].
- Totals use `size_t` [code].

**Measured.** Not measured separately. Cumulative run opt6: 344.3 s. It is slower than opt5 (335.0 s). The difference is inside the run-to-run noise [bench].

## 7. GibbsSPR: partial post-order updates (a183303)

**What was slow.** Each `GibbsSPR()` ran three full two-pass updates: after the re-rooting (master:1521), after the detach (master:640) and after the attach (master:697) [code].

**What changed.**

- A new message `UPDATE_PARTIAL`. `GlobalPartialUpdateConditionalLikelihoods(target)` broadcasts the signal, then one integer (the link index), then sends the usual root likelihood request (PhyloProcess.cpp:1540-1549). The slave reads the integer and calls `PartialUpdateConditionalLikelihoods(GetLinkForGibbs(n)->GetNode())` (PhyloProcess.cpp:1834-1837) [code].
- Slave side (PhyloProcess.cpp:181-238): `MarkAncestors` marks the target node and all its ancestors. `PartialPostOrderPruning` goes down the tree and propagates only into marked children. At each marked node it runs Reset, Multiply over all child links and Offset. Then `MultiplyByStationaries`, `ComputeLikelihood` and the **full** `PreOrderPruning` run. If the target is not found, a full update runs [code].
- Callers: the re-rooting for `rep > 0` (PhyloProcess.cpp:1602-1612), the step after the detach (686-709) and the step after the attach (765-768) [code].

**Why the result is unchanged** [inference from tracing Tree.cpp:721-800 and 884-903; no separate test beyond the chain comparisons]:

- Vectors are indexed by link, and link objects persist across `RootAt`, `Detach` and `Attach`.
- Re-rooting: only the branches on the path between the old and the new root change direction. Those nodes are the marked ones. The off-path child vectors were post-order results of the previous update, for the same link, branch and subtree.
- Detach: the pruned subtree keeps its vector (same subtree, same branch). Only the new parent and its ancestors need recomputation.
- Attach: the scan overwrote one vector (PhyloProcess.cpp:785). The node of `up` is marked, so it is recomputed.
- Every recomputed node runs the same statement sequence as the full `PostOrderPruning`. The pre-order pass is complete.
- The first repetition uses a full update, so the argument has a starting point.
- The random-number order on the master is unchanged [code].

**What a reviewer should check.**

- The precondition "everything was up to date" is the caller's duty. It holds in `GibbsSPR(int)` because the repetition with `rep == 0` does a full update after the re-rooting [code].
- The mark array is indexed by `Node::GetIndex()` with a range check that exits (PhyloProcess.cpp:204-207) [code].
- The comment above `GibbsSPR(int)` said that each call "ends with a complete update", while it now ends with a partial update. The comment is corrected (section 12, item 2) [code, rechecked].
- The two broadcasts and the reply have the same shape as existing messages. No deadlock path was found [code, inference].

**Measured.** Not measured separately. Cumulative run opt7: 303.6 s (237.4 s with `-mavx2`) [bench].

## 8. Optional per-message statistics (72311b3, `-DPB_MPI_STATS`)

Not an optimization. A second version of `WaitLoop` (PhyloProcess.cpp:1657-1713) is compiled only when `PB_MPI_STATS` is defined. It counts each message type and records execution time and idle time with `MPI_Wtime`. When the slave receives `KILL`, it writes `mpistats.<rank>.txt` in the current directory. The normal `WaitLoop` (1715-1722) is unchanged. It makes no random-number call and sees the same messages [code, rechecked].

Limits: arrays of size 128 with bounds checks; no file is written on abnormal exit; the table of names ends with an extra name `BATCH` that the enum does not have (see section 12).

## 9. `MoveDirWeights`: cached logarithms and log-gamma values (e62498b)

**What was slow.** `MixtureProfileProcess::MoveDirWeights` runs on the master while all slaves wait. It calls `LogStatPrior()` twice per (repetition, component). Each call evaluated, for every component, `dim` logarithms and `dim+1` calls of `logGamma` (master:MixtureProfileProcess.cpp:259-268) [code].

**What changed.**

- `MoveDirWeights` surrounds its loop with `PrepareLogProfileCache(GetNcomponent())` and `ReleaseLogProfileCache()` (MixtureProfileProcess.cpp:335, 353) [code].
- `LogStatPrior(cat)` uses the cached logarithms when present. It always uses the new `CachedLogGamma(slot, x)` (263-273). That function keeps the last `(x, value)` per slot in `loggammakey` and `loggammaval` (MixtureProfileProcess.h:208-210) and recomputes when `key != x` [code].

**Why the result is unchanged.** `Random::logGamma` is a pure function of its argument (Random.cpp:533-547). The cached value is the value for the identical double. The order of additions in `LogStatPrior(cat)` is unchanged. The random draws are unchanged [code].

**What a reviewer should check.**

- The key sentinel is `-1.0`. For `x == -1.0`, `master` exits inside `logGamma` ("only positive argument"). The new code returns the initial value 0.0. Dirichlet weights are positive in normal use [code, inference].
- The log-gamma cache is active for all callers of `LogStatPrior(cat)`, not only for this move [code].
- `PoissonMixtureProfileProcess::MoveDirWeights` (.cpp:170) is a separate non-virtual function and is unchanged. Which version a model calls depends on static name lookup. This was not traced for each model [code, inference].

**Measured.** Not measured separately.

## 10. `Propagate` kernel: vector accumulators, leaf cache, hash-free exponential cache (312c2bd)

**What was slow.** After sections 2 and 5: one hash lookup per site and rate; accumulators kept in memory by the auto-vectorized loops; a leaf result recomputed for every site that shares the same matrix, rate and observed state [code].

**What changed.**

- `MatVecBlock20` (Propagate.cpp:132-165) computes 20 outputs with 5 accumulators of 4 doubles (when `__AVX__` is defined) or 10 accumulators of 2 doubles. `MatVecCols` (167-184) runs blocks of 20, then a scalar tail. For 20 states there is one block. For 61 states (codons) there are 3 blocks and a scalar tail of 1. For 4 states only the scalar tail runs [code, rechecked].
- Exponential cache: `SubMatrix::propstamp` and `propslot` (SubMatrix.h:116-117) are scratch fields. If `matrix->propstamp != this->propstamp`, the matrix receives a new slot in `propheads` (Propagate.cpp:227-231). Each slot heads a linked list of `PropEntry{rate, expoff, next, leafbase}` (233-256). The process stamp is incremented at each call (211). This replaces the `std::unordered_map` of 9cb60c2 [code].
- Leaf cache: `leafindex[leafbase+obs]` points into `leafbuffer`. The first site computes the vector, clips negative values and stores it together with its clip count in `leafclipped`. Later sites copy the vector and add the same count to `infprobcount` (Propagate.cpp:262-314) [code].
- General path: the range checks are done without branches (`negative |=` and `bad |= !(down[k] >= 0.0)`). The old code runs only in the rare error case (316-369) [code].

**Why the result is unchanged.**

- Each vector lane is one output element. Each accumulator starts at 0.0 and adds `col[k]*x[l]` for `l = 0..n-1` in order. This is the same sequence per element as the scalar loop and as `master`. There is no horizontal reduction and no reassociation [code, rechecked].
- The multiplication and the addition are separate operations. Without FMA contraction, a vector multiply and add give the same value as the scalar ones (see `build-and-compiler-options.md`, section 3.2) [code, inference].
- A copied leaf vector equals a recomputed one, because the inputs (matrix, rate, time, observed state) are bit-equal. `infprobcount` receives the same total [code, inference].
- `down[nstate] = up[nstate]` is still done for every site (371) [code].
- Evidence by comparison of chains: 312c2bd matches `master` on 13 models (see `validation.md`).

**What a reviewer should check.**

- Stamp collision. If two `MatrixSubstitutionProcess` objects in one rank share a `SubMatrix` and reach equal stamps, `propslot` could index outside `propheads`. Whether the multi-gene models share matrices was not verified [inference].
- `int& idx` (line 270) stays valid because `leafindex` is not resized before `idx = newidx` (the comment at line 305 says so) [code].
- `expdiag` is taken after the resize of `expbuffer` [code].
- `aux` is a private buffer, so the `__restrict__` contracts hold for `(inv, up, aux)` and `(eig, aux, down)` [code].
- An observed state at or above `nstate` would index outside `leafindex`. `master` would also write out of range in `Initialize` [inference].
- Models with one matrix per site get no cache hits and pay the bookkeeping. This was not measured [inference].

**Portability.**

- The kernel uses GCC and clang vector extensions (`__attribute__((vector_size, aligned(8), may_alias))` and a compound-literal initializer, Propagate.cpp:120-130). The guard is `#if defined(__GNUC__) && !defined(PB_NO_VECTOR_KERNEL)` (line 118) [code, rechecked].
- A scalar fallback with the same order of operations exists (lines 186-199). It is used for compilers that do not define `__GNUC__`, or with `-DPB_NO_VECTOR_KERNEL` [code, rechecked].
- The code is not specific to x86. Without `__AVX__` it uses 16-byte vectors, which GCC maps to SSE2 or NEON [inference; not built on ARM].
- `aligned(8)` allows unaligned loads, so no alignment is required [code].
- No `std::unordered_map` is left in these files [code].

**Measured.** The code at 31a9b87 (the commit before this one; 31a9b87 changed only a script) and the code at 312c2bd were timed in one comparison, 3 unpinned runs each, default flags: 265.5 s against 170.3 s (1.56x). Individual runs: 309.1, 254.6, 232.7 s and 179.4, 169.8, 161.6 s. The two ranges do not overlap, so this commit does give a real gain on this model. This is the only commit with a separate measurement [bench].

## 11. Not part of this branch

- **182d90b** (batching of branch-length messages) was reverted by 7683baa. Chains were identical, but there was no speed-up at 8 and 64 ranks. The code is on branch `mpi-batching`. The net change in `perf-opt` is none [code].
- **5941ab9** (branch `hoist-nstate`) is not merged. See `PERFORMANCE.md`, section 3.6.
- Script and documentation commits are not described here: 515579a, 96a4fad, 288cc6e, fd3da8a, 38aa626, 377973f, 50612b8, 51fdf5f, 31a9b87, f88fb3b, and the documentation commits after it. The corrections 9f2a6fc and db61cc5 are described in pull request 39. Commit f8cca16 is mostly scripts but also adds default member initializers for the unfold flags in `PhyloProcess.h`.

## 12. Items to fix or decide before a pull request

Items 1, 2, 3 and 6 are fixed (2026-11-20). The four fixes change no behavior: the build has no `-Wreorder` warning any more, and all 14 models give the same chains as before (`validation.md`).

1. **Fixed: `-Wreorder` warning (87 warning lines in the build of 312c2bd, 0 in the `master` build).** In `MixtureProfileProcess.h:29` the constructor initialized `profile(0), logprofilecachesize(0)`, but `logprofilecachesize` is declared before `profile`. The behavior was correct. The initializer list now follows the declaration order. The pinned builds 312c2bd and 4d7ef1e were made before this fix [code, rechecked].
2. **Fixed: outdated comment** above `GibbsSPR(int)` (see section 7) [code, rechecked].
3. **Fixed: `pb_message_names`** ended with `,BATCH`, a leftover of the reverted batching work. The enum has no `BATCH`, and that index was never counted. The name is removed [code].
4. **`GibbsSPR(0)` return value** changed from NaN to 0 (section 4), and it is not in the commit message. Mention it in the pull request.
5. **`CachedLogGamma` sentinel** `-1.0` (section 9). Consider a separate validity flag.
6. **Fixed: `PB_NO_VECTOR_KERNEL`** was not documented in the Makefile. A Makefile comment now describes it, in addition to `build-and-compiler-options.md`.
7. **Non-virtual copies** of virtual function bodies (`PostOrderUpdateConditionalLikelihoods`, `PartialUpdateConditionalLikelihoods`) can drift from `UpdateConditionalLikelihoods` in a subclass (section 1).
8. **Public scratch fields** `SubMatrix::propstamp` and `propslot` are used by another class. They carry the shared-matrix hazard of section 10 [code, inference].
9. **Prepare without benefit.** The Matrix SBDP and Finite allocation moves call Prepare even when the model's `LogStatProb` does not read the cache (section 3) [code].
10. **`-mutsel -ncat 3` is not slower.** One unrepeated, unpinned run showed 56.1 s for `master` against 67.2 s. Repeated pinned with 30 cycles and 3 runs each, `master` took 147, 171 and 154 s and 312c2bd took 134, 137 and 103 s, with identical chains [bench]. The gain is small for this model. Section 10 lists a possible cause: per-site matrices get no cache hits and pay bookkeeping. This is a guess [inference].
11. **Formatting:** mixed tabs and spaces at RASCATFiniteGammaPhyloProcess.h:306-312 [code].
12. No leftover debug output was found in the diff [code].

## 13. Behavior differences other than speed

- No new command-line option, no changed default and no change to the normal output files [code].
- New compile-time options: `PB_MPI_STATS` (extra files `mpistats.<rank>.txt`), `PB_NO_VECTOR_KERNEL`, and the Makefile variable `ARCHFLAGS` (empty by default) [code].
- MPI protocol: the enum values 0 to 65 are unchanged. `UPDATE_POSTORDER` (66) and `UPDATE_PARTIAL` (67) are appended. `UPDATE_PARTIAL` is followed by an extra integer broadcast. The master and all slaves must run the same build. A `master` slave does not know the new messages, and a `perf-opt` slave does not update in `UNFOLD` and relies on the `UPDATE` that follows. This is the normal case, because all ranks of one run start the same binary [code, rechecked for the enum].
- The master sends fewer `UPDATE` and likelihood requests per cycle [code].
- `GibbsSPR(0)` returns 0 instead of NaN (section 4).
- The choice of fatal-error message in `PropagateImpl` can differ (section 2).
- `CachedLogGamma(-1.0)` returns 0 instead of exiting (section 9).

## 14. Files changed under `sources/`

| File | Change |
|---|---|
| `Parallel.h` | two messages appended to the enum |
| `PhyloProcess.h`, `PhyloProcess.cpp` | unfold flags; `UPDATE_POSTORDER` and `UPDATE_PARTIAL`; partial update functions; GibbsSPR; `PropagateFromLeaf` call; `PB_MPI_STATS` wait loop |
| `MatrixPhyloProcess.cpp`, `GeneralPathSuffStatMatrixPhyloProcess.cpp` | `if (unfoldupdate)` guard |
| `RASCATFiniteGammaPhyloProcess.h`, `RASCATSBDPGammaPhyloProcess.h`, `RASCATGTRSBDPGammaPhyloProcess.h` | `GlobalUnfoldBeforeCollapse` for `rep < 4` |
| `Propagate.cpp`, `MatrixSubstitutionProcess.h` | new `PropagateImpl` and its scratch members |
| `SubMatrix.h`, `SubMatrix.cpp` | transposed eigenvector copies; `propstamp` and `propslot` |
| `SubstitutionProcess.h` | virtual `PropagateFromLeaf` |
| `SubstitutionProcess.cpp` | contiguous vectors; hoisted `GetNrate` and `GetNstate` |
| `MixtureProfileProcess.h`, `MixtureProfileProcess.cpp` | log-profile cache; `CachedLogGamma`; `MoveDirWeights` |
| `PoissonMixtureProfileProcess.cpp`, `ExpoConjugateGTRMixtureProfileProcess.cpp` | `LogStatProb` uses the cache |
| `PoissonSBDPProfileProcess.cpp`, `PoissonFiniteProfileProcess.cpp`, `MatrixSBDPProfileProcess.cpp`, `MatrixFiniteProfileProcess.cpp` | Prepare and Release around the allocation loops |
| `Makefile` | `ARCHFLAGS` |

## 15. Code paths not touched

- `PoissonSubstitutionProcess::Propagate` (the CAT-Poisson kernel). CAT-Poisson still gets the unfold, GibbsSPR, vector layout and mixture-move changes, which explains its smaller gain.
- `SitePropagate`, `SimuPropagate` and the per-site stepping-stone code (Stepping.cpp).
- NNI (NNI.cpp), the branch-length move and its message protocol.
- `NonMPIGibbsSPR` and `OldMPIGibbsSPR` (gene models).
- Substitution mapping, sufficient statistics, and the moves of rates, relative rates, kappa and profile values.
- `PoissonMixtureProfileProcess::MoveDirWeights` and `GeneralPathSuffStatMatrixMixtureProfileProcess::LogStatProb`.
- All `Move()` loops other than the three CAT classes still use the full `GlobalUnfold` (AACodonMutSel, CodonMutSel, RASCATGTRFinite, RASCATGTRDP, RASCATGamma and others). They still get the slave-side single update, the contiguous vectors, the matrix kernel (if matrix-based) and the GibbsSPR changes [code].

## 16. Cumulative timing for orientation

Single runs, CAT-GTR, 40x5000, 8 ranks, unpinned, default flags, 20 cycles, wall seconds [bench]. Identical builds differ by 6% to 21%, so steps smaller than about 20% are not resolved. The full history, with later commits, is in `PERFORMANCE.md`, section 3.3.

| Run | Commit | Wall (s) | Relative to `master` mean (927.1 s) |
|---|---|---|---|
| base_a, base_b | `master` | 954.0, 900.2 | 1.00x |
| opt1 | d2456fc | 647.9 | 1.43x |
| opt2 | 9cb60c2 | 682.7 | 1.36x |
| opt3 | f56c1f7 | 533.0 | 1.74x |
| opt4 | 288cc6e | 508.7 | 1.82x |
| opt5 | 14297a9 | 335.0 | 2.77x |
| opt6 | 81ad2c2 | 344.3 | 2.69x |
| opt7 | a183303 | 303.6 | 3.05x |

For all runs listed, `chain.treelist` is byte-identical to the `master` run, and `chain.trace` is identical in all columns except the two timing columns [bench].

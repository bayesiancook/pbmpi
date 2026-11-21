# Pinning MPI ranks to physical cores

> Part of the performance documentation; see the scope note at the top of `PERFORMANCE.md`.

This note explains why `pb_mpi` must run with one rank per physical core, how to do it on Slurm with Open MPI, and what was measured.

## 1. Summary

- On a node with hyper-threading (SMT), a Slurm job that asks for N CPUs gets N *hardware threads*. N MPI ranks then share N/2 physical cores.
- `pb_mpi` is a synchronous master/slave program. At each step the master waits for the slowest slave. Ranks that share a core run slower, and every other rank waits for them.
- Measured: the same chain (same seed, same build, same data) took **147.8 s** unpinned and **95.0 s** pinned (mean of 3 runs, 8 ranks). That is **1.56x**.
- Pinning does not change the result. The two chains are identical in `.treelist`, `.chain` and `.trace`.
- Use `sbatch --hint=nomultithread` and `mpirun --bind-to core --map-by core`.
- Compare timings only between runs that use the same placement. Unpinned timings are noisy and about 1.5x slower.

## 2. Background: CPUs, cores and threads

The UCR HPCC "milan" nodes used for these measurements have:

| Item | Value (from `lscpu` on the node) |
|---|---|
| CPU | AMD EPYC 7713 64-Core Processor |
| Sockets | 2 |
| Cores per socket | 64 |
| Threads per core | 2 |
| Logical CPUs | 256 |

Slurm counts logical CPUs. By default, `sbatch -n 8` allocates 8 logical CPUs. Two logical CPUs can be the two threads of one physical core. In that case 8 ranks run on 4 physical cores. This was the situation for all early benchmark runs in this project (the "8-rank" timings in the history table of `PERFORMANCE.md`).

With `--hint=nomultithread`, Slurm gives each task one whole physical core and reserves the sibling thread. The job then shows twice as many CPUs as it has ranks. For example, a 32-rank job shows `NumCPUs=64`.

## 3. Why unpinned runs are slower

What the measurements show:

1. **The ranks do not run at equal speed when unpinned.** `-DPB_MPI_STATS` records, for each slave rank, the total time spent executing messages. Same chain, 8 ranks (1 master + 7 slaves), `-mavx2`, 20 cycles, commit 312c2bd:

   | Slave rank | Execution, unpinned (s) | Waiting, unpinned (s) | Execution, pinned (s) | Waiting, pinned (s) |
   |---|---|---|---|---|
   | 1 | 127.1 | 8.1 | 99.9 | 6.4 |
   | 2 | 104.2 | 31.0 | 99.5 | 6.8 |
   | 3 | 127.8 | 7.4 | 99.6 | 6.7 |
   | 4 | 126.7 | 8.5 | 99.2 | 7.1 |
   | 5 | 127.0 | 8.2 | 99.4 | 6.9 |
   | 6 | 112.3 | 22.9 | 99.4 | 6.9 |
   | 7 | 113.6 | 21.6 | 99.8 | 6.5 |

   (Values are the `TOTAL` row of `mpistats.<rank>.txt`: execution time, and time spent idle after a message.) Pinned, all seven ranks agree within 0.7 s in execution time. Unpinned, execution time ranges from 104.2 s to 127.8 s. The slowest unpinned ranks (1, 3, 4, 5) are about 1.28x slower than the pinned ranks.

2. **The fast ranks wait.** Idle time is much larger on the faster unpinned ranks (rank 2: 31.0 s, rank 6: 22.9 s, rank 7: 21.6 s) than on any pinned rank (6.4 to 7.1 s). The master needs an answer from every slave, so a fast rank waits for the slowest one.

3. **Wall time follows the slowest rank.** In these statistics runs, wall time was 136.4 s unpinned and 106.7 s pinned, a ratio of 1.28, the same as the ratio of the slowest unpinned rank to a pinned rank. (The statistics build is slower than the normal build because of the timing code, and it is one run per placement.)

What is *not* measured (inference):

- The most likely cause of the slow ranks is that two ranks share one physical core and compete for its execution units, vector units and caches. The data are consistent with this explanation: with 8 ranks on 4 cores, some ranks run alone on a core and some share one, and the measured execution times differ between ranks. The exact rank-to-core placement was not recorded, so the match is not proven. No hardware counters were collected.
- Idle ranks poll for messages (a busy wait). A profile of the unpinned run shows about 15% of slave cycles and about 91% of master cycles inside Open MPI polling. A polling rank on the sibling thread takes cycles from its neighbour. This effect was not measured separately.
- Whether the master (rank 0) needs a core of its own was not tested.

## 4. Measurements

Setup for all rows: CAT-GTR (`-cat -gtr`), simulated alignment 40 taxa x 5000 amino-acid sites, `-rnd 1`, 20 saved cycles, 8 ranks, AMD EPYC 7713 "milan" node, GCC 8.5.0, Open MPI 4.1.2, build flag `-mavx2`. Wall time includes MPI start-up. Each row is 3 independent jobs.

| Commit | Placement | Wall time (s) | Mean | Median |
|---|---|---|---|---|
| 312c2bd | unpinned (`sbatch -n 8`) | 151.1, 145.0, 147.1 | 147.8 | 147.1 |
| 312c2bd | pinned | 99.8, 94.0, 91.1 | 95.0 | 94.0 |
| 5941ab9 (branch `hoist-nstate`, not merged) | unpinned | 137.1, 134.6, 136.5 | 136.0 | 136.5 |
| 5941ab9 (branch `hoist-nstate`, not merged) | pinned | 87.5, 95.5, 104.6 | 95.9 | 95.5 |

Notes:

- Unpinned versus pinned for 312c2bd: 1.56x (mean) and 1.57x (median).
- Pinned runs still vary by up to 14% (87.5 s to 104.6 s for 5941ab9). Pinning removes the large systematic effect. It does not remove all noise.
- Identity check: `scripts/bench/compare_runs.sh` on the unpinned and pinned run of 312c2bd reports identical `.treelist`, `.chain` and `.trace`.
- On the pinned runs, the difference between 312c2bd and 5941ab9 is smaller than the run-to-run spread. This is why 5941ab9 was not merged.
- Earlier variant comparisons in this project (compiler flags, kernel variants) used unpinned runs. Their differences are within the noise of unpinned runs and must not be used to rank variants. See `build-and-compiler-options.md`.

## 5. How to pin

### 5.1 Slurm job header

```
#SBATCH -N 1
#SBATCH -n <number of MPI ranks>
#SBATCH --hint=nomultithread
```

`--hint=nomultithread` makes Slurm allocate whole physical cores. Check with `scontrol show job <id>`: `NumTasks` is the number of ranks, and `NumCPUs` is twice that on an SMT-2 node.

### 5.2 mpirun

For Open MPI:

```
mpirun --bind-to core --map-by core -np $SLURM_NTASKS pb_mpi ...
```

- `--map-by core` places rank i on core i.
- `--bind-to core` stops the operating system from moving a rank to another core.
- Add `--report-bindings` once to confirm the placement. Open MPI prints one line per rank.

For the benchmark script:

```
MPIRUN_OPTS="--bind-to core --map-by core" sbatch --hint=nomultithread scripts/bench/bench.sbatch <ref> <id>
```

`scripts/bench/scaling.sbatch` already uses both options.

### 5.3 Other MPI libraries

Only Open MPI 4.1.2 was used. Equivalent options exist in other libraries (for example `-bind-to core` in MPICH and `I_MPI_PIN` variables in Intel MPI), but they were not tested here. When running under Slurm without Open MPI, `srun --cpu-bind=cores` is the usual alternative. This is also untested here.

## 6. Pitfalls

- **Rank count and results.** All ranks start from the same seed, and the alignment sites are divided among the slaves by rank. A chain therefore depends on the number of ranks. Changing the rank count changes the chain. Pinning does not. (This statement comes from reading `PB.cpp` and `PhyloProcess.cpp`. A run with different rank counts was not made to confirm it.)
- **Memory.** `--hint=nomultithread` doubles the CPU count that Slurm bills, but not the memory request. Set `--mem` explicitly.
- **Billing and queue.** A pinned job uses twice the logical CPUs per rank. Cluster accounting and fair-share charge for them.
- **Shared nodes.** If another job runs on the same node, its threads can still compete for memory bandwidth and caches. Pinning does not give isolation. Use `--exclusive` for benchmark runs when the cluster allows it.
- **More ranks than physical cores.** A request with more ranks than physical cores cannot be pinned one rank per core. Do not exceed the physical core count (128 on a milan node).
- **Hybrid or threaded runs.** `pb_mpi` has no threads. These options assume one single-threaded process per rank.

## 7. Recommendations

1. Run every production chain and every benchmark pinned.
2. Record the placement (`--hint`, `MPIRUN_OPTS`) with each timing.
3. Repeat each benchmark at least 3 times. Report the median and the range.
4. Do not compare an unpinned timing with a pinned one.
5. Re-measure rank scaling (see `PERFORMANCE.md`) with the current build. Most scaling points were measured with an earlier build.

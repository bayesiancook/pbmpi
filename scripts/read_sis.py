import math
from basicstat import *

def readbf(runnames, min_ess = 10):

    logz_idx = 2
    meanlog_idx = 3
    varlog_idx = 5
    npoint_idx = 6
    ess_idx = 7

    runs = []
    for run in runnames:
        with open(run, 'r') as infile:
            x = [[float(y) for y in line.rstrip('\n').split()] for line in infile]
            runs.append(x)

    nrun = len(runs)
    m = len(runs[0])
    for i in range(1,nrun):
        if len(runs[i]) != m:
            print("error: run {0} does not have same length as run 0: {1} versus {2}".format(runnames[i], len(runs[i]), m))
            sys.exit(1)

    # within-run variances
    across_run_vars = [unbiased_var([x[i][logz_idx] for x in runs]) for i in range(m)]
    mean_across_run_var = mean(across_run_vars)
    mean_bias = -0.5*mean_across_run_var

    logzs = [mean([x[i][logz_idx] for i in range(m)]) for x in runs]
    mean_logz = mean(logzs)
    mean_debiased_logz = mean_logz - mean_bias

    ess = [mean([x[i][ess_idx] for x in runs]) for i in range(m)]
    # ess = [mean([x[i][ess_idx] for i in range(m)]) for x in runs]
    mean_ess = mean(ess)
    nmin_ess = sum([s<min_ess for s in ess])

    totvar = mean_across_run_var/m
    stdev = math.sqrt(totvar/nrun)
    dev = student95critval(nrun)
    min_score = mean_debiased_logz - dev*stdev
    max_score = mean_debiased_logz + dev*stdev

    return (mean_logz, mean_debiased_logz, min_score, max_score, mean_bias, stdev, mean_ess, nmin_ess, m)

if __name__ == "__main__":

    import sys

    if len(sys.argv) == 1:
        print("readbf.py run_list")
        sys.exit(1)

    runnames = sys.argv[1:]

    (mean_raw_score, mean_score, min_score, max_score, mean_bias, stdev, mean_ess, nmin_ess, m) = readbf(runnames)

    print("{0:^5s} {1:^10s} {2:^10s} {3:^10s} {4:^10s} {5:^10s} {6:^10s} {7:^10s} {8:^10s}".format("n", "score", "debiased", "CI95min", "CI95max", "bias", "stdev", "mean ess", "#ess<10"))
    print("{0:5d} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f} {6:10.4f} {7:10.4f} ({8}/{9})".format(m, mean_raw_score, mean_score, min_score, max_score, mean_bias, stdev, mean_ess, nmin_ess, m))


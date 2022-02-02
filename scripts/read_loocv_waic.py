import math
from basicstat import *
from studentC import student_table

def readwaic(runnames, min_ess = 10):

    logl_idx = 1
    var_idx = 2
    loocv_idx = 3
    essloocv_idx = 4
    logpostmeanl_idx = 5
    esspostmeanl_idx = 6

    runs = []
    for run in runnames:
        with open(run, 'r') as infile:
            header = infile.readline()
            x = [[float(y) for y in line.rstrip('\n').split()] for line in infile]
            runs.append(x)

    nrun = len(runs)

    m = len(runs[0])
    for i in range(1,nrun):
        if len(runs[i]) != m:
            print("error: run {0} does not have same length as run 0: {1} versus {2}".format(runnames[i], len(runs[i]), m))
            sys.exit(1)

    ####
    # LOO-CV 
    ####

    # averaging across sites for each run
    loocvs = [mean([x[i][loocv_idx] for i in range(m)]) for x in runs]
    # then taking the mean and the variance across runs
    loocv = mean(loocvs)
    var_loocv = unbiased_var(loocvs)

    # bias calculation based on the site-specific variances across runs
    across_run_vars_loocv = [unbiased_var([x[i][loocv_idx] for x in runs]) for i in range(m)]
    mean_across_run_var_loocv = mean(across_run_vars_loocv)
    # sitevar_loocv = mean_across_run_var_loocv
    bias_loocv = 0.5*mean_across_run_var_loocv
    debiased_loocv = loocv - bias_loocv

    # site-specific ess
    ess_loocv = [mean([x[i][essloocv_idx] for x in runs]) for i in range(m)]
    # then averaged across sites
    mean_ess_loocv = mean(ess_loocv)
    # and counting those with low ess
    fmin_ess_loocv = abs(mean([(ess_loocv[i]<min_ess)*x[i][loocv_idx] for i in range(m)]) / loocv)
    nmin_ess_loocv = sum([s<min_ess for s in ess_loocv])

    ####
    # waic (official version, based on Wanatabe)
    ####

    # averaging log<L>_post across sites, for each run
    logpostmeanls = [mean([x[i][logpostmeanl_idx] for i in range(m)]) for x in runs]
    
    # mean and variance of <logL>_post across sites, for each run
    postmean_logls = [mean([x[i][logl_idx] for i in range(m)]) for x in runs]
    postvar_logls = [mean([x[i][var_idx] for i in range(m)]) for x in runs]

    # wAIC scores across runs
    waics = [logpostmeanls[i] - postvar_logls[i] for i in range(nrun)]
    # then taking the mean and variance across runs
    waic = mean(waics)
    var_waic = unbiased_var(waics)

    # bias calculation based on the site-specific variances across runs
    across_run_vars_logpostmeanl = [unbiased_var([x[i][logpostmeanl_idx] for x in runs]) for i in range(m)]
    mean_across_run_var_logpostmeanl = mean(across_run_vars_logpostmeanl)
    bias_waic = -0.5*mean_across_run_var_logpostmeanl
    debiased_waic = waic - bias_waic

    # site-specific ess
    ess_waic = [mean([x[i][esspostmeanl_idx] for x in runs]) for i in range(m)]
    # then averaged across sites
    mean_ess_waic = mean(ess_waic)
    # and counting those with low ess
    sitewaics = [mean([x[i][logpostmeanl_idx] - x[i][var_idx] for x in runs]) for i in range(m)]
    fmin_ess_waic = abs(mean([(ess_waic[i]<min_ess)*sitewaics[i] for i in range(m)]) / waic)
    nmin_ess_waic = sum([s<min_ess for s in ess_waic])

    return (debiased_waic, bias_waic, var_waic, mean_ess_waic, nmin_ess_waic, fmin_ess_waic, debiased_loocv, bias_loocv, var_loocv, mean_ess_loocv, nmin_ess_loocv, fmin_ess_loocv, m, nrun)

if __name__ == "__main__":

    import sys

    if len(sys.argv) == 1:
        print()
        print("model_fit.py <run_list>")
        print("<run_list>: list of .sitelogl files for at least 2 independent runs under same model")
        print()
        sys.exit(1)

    runnames = sys.argv[1:]

    nrun = len(runnames)
    if (nrun < 2):
        print("at least 2 independent runs are needed")
        sys.exit(1)

    (waic, bias_waic, var_waic, mean_ess_waic, nmin_ess_waic, fmin_ess_waic, loo, bias_loo, var_loo, mean_ess_loocv, nmin_ess_loocv, fmin_ess_loocv, m, nrun) = readwaic(runnames)

    a = student_table[nrun-1]

    stdev_waic = math.sqrt(var_waic)
    min_waic = waic - a*stdev_waic
    max_waic = waic + a*stdev_waic

    stdev_loo = math.sqrt(var_loo)
    min_loo = loo - a*stdev_loo
    max_loo = loo + a*stdev_loo

    print("{0:10s} {1:^10s} {2:^10s} {3:^10s} {4:^10s} {5:^10s} {6:^10s} {7:^10s} {8:^10s}".format("", "debiased score", "bias", "stdev", "CI95min", "CI95max", "ess", "%(ess<10)", "f(ess<10)"))
    print("{0:10s} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f} {6:10.4f} {7:10.3f} {8:10.3f}".format("LOO-CV", loo, bias_loo, stdev_loo, min_loo, max_loo, mean_ess_loocv, nmin_ess_loocv/m, fmin_ess_loocv, m))
    print("{0:10s} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f} {6:10.4f} {7:10.3f} {8:10.3f}".format("waic", waic, bias_waic, stdev_waic, min_waic, max_waic, mean_ess_waic, nmin_ess_waic/m, fmin_ess_waic, m))


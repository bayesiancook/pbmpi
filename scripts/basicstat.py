
# make scripts more portable (no need to use numpy or scipy)

from studentC import student_table

def mean(x):
    return sum(x)/len(x)

def var(x):
    m1 = mean(x)
    m2 = mean([y*y for y in x])
    return m2 - m1*m1

def unbiased_var(x):
    m = len(x)
    return m/(m-1)*var(x)

def covar(x1, x2):
    m1 = mean(x1)
    m2 = mean(x2)
    m12 = mean([x1[i]*x2[i] for i in range(len(x1))])
    return m12 - m1*m2

def correl(x1, x2):
    m1 = mean(x1)
    m11 = mean([y*y for y in x1])
    m2 = mean(x2)
    m22 = mean([y*y for y in x2])
    m12 = mean([x1[i]*x2[i] for i in range(len(x1))])
    cov12 = m12 - m1*m2
    var1 = m11 - m1*m1
    var2 = m22 - m2*m2
    return cov12 / math.sqrt(var1*var2)

def student95critval(df=10):

    if (df < 2):
        print("error in student95critval: df should be at least 2")
        return 2

    return student_table[df-1]


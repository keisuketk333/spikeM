#!/bin/python

##############################################################
#    Authors: Yasuko Matsubara & Christos Faloutsos
#    Date: 04-26-2014
##############################################################

##################
# Given: T, (parameter set):
#    T,
#    N, beta0, slope,
#    nc, Sc, bgnoise,
#    Pp, Pa, Ps
#
# Output: dB, B, U
############
# generate the dB[T], B[T] and U[T] arrays
# to simulate the rise-and-fall phenomenon in a blog setting
#      dB[n]    count of people that started blogging @n
#      U[n]     count of un-informed people
#      B[n] = sum(dB[n]): total number of informed people
#                        (who immediately start blogging
############
# Our model:
#    - power law decay functions
#      decay(n,beta,slope)  beta * n**(slope)
#                           n=1,2, ... slope=e.g., -1.5
#    - periodicity (period() could return a sinusoid/clock function of time
#    - a constant source of information, Sb (like, say, a TV, as in Bass)
#####################

import sys
import math as M

try:
    from scipy.stats import poisson
except:
    print("can not find scipy")
    print("for mac users:")
    print("sudo port install py27-scipy")
    raise

try:
    import pylab as P

    has_pylab = True
except:
    print("sudo apt-get install pylab      matplotlib ")
    print("    continuing without plotting functionality")

    has_pylab = False

import array as A


############################
# power law decay, with slope (e.g.,-1.5)
############################
def decay_pl(n, beta, slope):
    if n < 1:
        return 0
    else:
        return beta * M.pow(n, slope)


##########################
# exogenous function
##########################
def exo(t, nc, Sc):
    if t == nc:
        return Sc
    else:
        return 0


##########################
# sinusoidal periodicity
##########################
def period(t, Pp, Pa, Ps):
    if Pa == 0:
        return 1
    val = M.sin((2 * M.pi / Pp) * (t + Ps))
    val = (val + 1) * 0.5
    val = (1 - (val * Pa))
    # val = 1 + Pa*M.cos( (2*M.pi/Pp)*(t-Ps) )
    return val


def spikeM(T, N, beta0, slope, nc, Sc, bgnoise, Pp, Pa, Ps):
    slope = -abs(slope)
    nc = round(nc)

    # --- for multi-wdsize (please ignore)
    wd = 1
    T = T * wd
    nc = int(nc) * wd
    bgnoise = bgnoise / wd
    Sc = Sc / wd
    # --- for multi-wdsize (please ignore)

    U = [N] * T  # uninfected
    dB = [0] * T  # infected (total)
    B = [0] * T  # blogged/infected

    # init
    B0 = 0
    U[0] = N - B0
    dB[0] = B0
    B[0] = B0

    for n in range(0, T - 1):
        dsum = 0  # ~ number of sneezes
        for i in range(nc, n + 1):  # i.e, from zero to n
            Si = exo(i, nc, Sc)
            dsum += (dB[i] + Si) * decay_pl(n + 1 - i, beta0, slope)
        P_n1 = period(n + 1, Pp, Pa, Ps)
        dB[n + 1] = P_n1 * (U[n] * dsum + bgnoise)

        # upper-bound b[n+1]:
        if dB[n + 1] > U[n]:
            dB[n + 1] = U[n]

        U[n + 1] = U[n] - dB[n + 1]
        B[n + 1] = B[n] + dB[n + 1]

        # print(B[n+1], U[n+1], N, abs(B[n+1] + U[n+1] - N))
        assert (abs(B[n + 1] + U[n + 1] - N) < 0.001)

    # --- for multi-wdsize (please ignore)
    b_tmp = list(range(T))
    for j in range(0, int(T / wd)):
        b_tmp[j] = 0
        for jj in range(0, wd):
            b_tmp[j] += dB[j * wd + jj]
    dB = b_tmp
    T = int(T / wd)
    # --- for multi-wdsize (please ignore)

    dB = dB[0:T]
    B = B[0:T]
    U = U[0:T]
    return dB, U


if __name__ == "__main__":

    T = 100  # number of time-ticks
    N = 2000  # initial population
    betaN = 1.0
    beta0 = betaN / N
    slope = -1.5
    nc = 0
    Sc = 1.0
    bgnoise = 0.0

    Pp = 10
    Pa = 0.0
    Ps = 0

    (dB, U) = spikeM(T, N, beta0, slope, nc, Sc, bgnoise, Pp, Pa, Ps)

    if has_pylab:
        tlist = list(range(0, len(dB)))
        xlist = dB
        ylist = U
        ptitle = " N=%d beta*N=%.2f T=%d steps " % (N, beta0 * N, T)
        P.title(ptitle)
        P.xlabel("time")
        P.ylabel("num bloggers @ time ")
        P.axis([1, T, 0, max(U) * 10])
        P.xscale('log')
        P.yscale('log')
        P.plot(tlist, xlist)
        P.plot(tlist, ylist)
        P.show()
    else:
        print("no pylab - do 'make demo' to see gnuplot demo")

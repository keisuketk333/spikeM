# -*- coding: utf-8 -*-

##################
# Given: dat:
# dat:   input sequence
# pfreq: period of the cycle (e.g., pfreq=24 hours)
#        (if you don't need freq, set pfreq=-1;)
# outfn: output file name
# iter:  # of max iteration
# wantPlot:  (=1)if you need GUI plot, (=0) else
##################

import numpy as np
import math
import matplotlib.pyplot as plt
import spikeM
try:
    import lmfit
except ModuleNotFoundError:
    print("can not find lmfit - please see http://lmfit.github.io/lmfit-py/")

XTL = 1.e-7  # どの程度探索を進めるか決定する閾値
FTL = 1.e-7
MAXFEV = 100


def main():
    # # demo: fitting sequence "./sequence.dat"
    fn = './sequence.dat'
    # duration of sequence
    T = 24 * 4
    # # of max iteration
    ITER = 20
    # daily periodicity (24hours)
    pfreq = 24

    with open(fn) as f:
        dat = [float(s.strip()) for s in f.readlines()]

    dat = dat[0:T]
    outfn = 'output'
    wantPlot = 0  # no

    print('===================================')
    print('DEMO - fitting sample sequence')
    print('-----------------------------------')
    print('- filename = ' + fn)
    print('- duration = ' + str(T))
    print('- max iteration = ' + str(ITER))
    print('===================================')
    print(' ')
    spikeMfit(dat, pfreq, outfn, ITER, wantPlot)


def spikeMfit(dat, pfreq, outfn, iter, wantPlot):
    # fitting spikeM
    RSE, params = LMfit(dat, pfreq, iter, wantPlot)
    # plot fitting results
    plotsRNF(dat, params, outfn)
    # save spikeM parameters
    with open(outfn + '.param', mode='w') as f:
        f.write(str(params))


def LMfit(dat, pfreq, iter, wantPlot):
    # # parameter settings
    # duration of sequence
    T = len(dat)
    # init parameters
    params = init_params(dat, pfreq)
    LBs, UBs = init_const(dat, pfreq, T)
    # default
    order = []
    order = order + [3, 4]  # nc, Sc
    order = order + [1, 0]  # BetaN, N_max
    if pfreq != -1:
        order = order + [8, 7]  # pshift, p_rate
    order = order + [5]  # background noise

    #  # --- start fitting --- # #
    RSEP = np.inf
    RSE = np.inf
    for i in range(1, iter+1):
        # ftype = 'lin'
        # # if you want "tail-part-sensitive" fitting, please try below
        # if i < iter / 2 or i == iter:
        #     ftype='log'
        # else:
        #     ftype='lin'
        if i < iter/2:
            ftype = 'log'
        else:
            ftype = 'lin'
        params0 = params
        print('iter=' + str(i) + '/' + str(iter))
        # # if you want to plot
        if wantPlot == 1:
            plotsRNF(dat, params, [])  # 未完成
        # # for each param
        for loc in range(len(order)):
            lo = order[loc]
            if lo == 5 and ftype == 'lin':
                ftype = 'log'
            # start lmfit
            # nb
            if lo == 3:
                params[lo] = FD_search(dat, params, lo, ftype)
            else:
                P = _createP(lo, params, dat, T)
                try:
                    lmsol = lmfit.Minimizer(F_RNF, P, fcn_args=(dat, params, lo, ftype))
                    res = lmsol.leastsq(xtol=XTL, ftol=FTL)
                    # res = lmsol.leastsq(xtol=XTL, ftol=FTL, max_nfev=MAXFEV)
                    params = _updateP(res.params, lo, params)
                except:
                    print('Debug:', lo, params[lo])
                    params[lo] = params0[lo]
                params = _const(params, LBs, UBs)

        # compute RSError
        RSE = printRNF(dat, params)
        TH = 0.0001
        if abs(RSEP-RSE) < TH or RSE < 0.0008:
            break
        RSEP = RSE
    if np.isnan(RSEP):
        params = np.zeros(9)
    return RSE, params


def init_params(dat, pfreq):
    # init params
    # RNF-base
    N_max = sum(dat)  # 0
    betaN = 1.0       # 1
    slope = -1.5      # 2
    # RNF-X
    nc = 0            # 3
    Sc = 0.1          # 4
    bgnoise = 0.01    # 5
    # RNF-P
    Pp = pfreq        # 6
    Pa = 0.1          # 7
    Ps = 1.0e-5       # 8
    #
    params = [N_max, betaN, slope, nc, Sc, bgnoise, Pp, Pa, Ps]
    return params


def _createP(lo, params, dat, T):
    P = lmfit.Parameters()
    namelist = ['N_max', 'betaN', 'slope', 'nc', 'Sc', 'bgnoise', 'Pp', 'Pa', 'Ps']
    # create params
    # P.add(namelist[lo], value=params[lo], vary=True)
    if lo == 0:
        P.add(namelist[lo], value=params[lo], vary=True, min=sum(dat))
    elif lo == 1:
        P.add(namelist[lo], value=params[lo], vary=True, min=0.01, max=2.0)
    elif lo == 3:
        P.add(namelist[lo], value=params[lo], vary=True, min=0.0, max=T/2)
    elif lo == 4:
        P.add(namelist[lo], value=params[lo], vary=True, min=0.0)
    elif lo == 5:
        P.add(namelist[lo], value=params[lo], vary=True, min=1.0e-8)  # min=0.0
    elif lo == 7:
        P.add(namelist[lo], value=params[lo], vary=True, min=0.05, max=1.0)
    elif lo == 8:
        P.add(namelist[lo], value=params[lo], vary=True, min=1.0e-8)
    return P


def _updateP(P, lo, params):
    namelist = ['N_max', 'betaN', 'slope', 'nc', 'Sc', 'bgnoise', 'Pp', 'Pa', 'Ps']
    params[lo] = P[namelist[lo]].value
    return params


def init_const(dat, pfreq, T):
    LB_base = np.array([sum(dat), 0.01, -1.5])
    UB_base = np.array([np.inf, 2.0, -1.5])
    LB_X    = np.array([0, 0.0, 0])
    UB_X    = np.array([T/2, np.inf, np.inf])
    LB_P    = np.array([pfreq, 0.05, 0])
    UB_P    = np.array([pfreq, 1, pfreq])
    LBs     = np.concatenate([LB_base, LB_X, LB_P])
    UBs     = np.concatenate([UB_base, UB_X, UB_P])
    return LBs, UBs


def _const(params, LB, UB):
    params = np.abs(params)
    # pshift
    params[8] = params[8] % params[6]
    # L & U bounding
    for i in range(len(params)):
        if params[i] < LB[i]:
            params[i] = LB[i]
        if params[i] > UB[i]:
            params[i] = UB[i]
    return params


def removeSparse(X, wd, th):
    # wd = int(Decimal(wd/2).quantize(Decimal('1'), rounding=ROUND_CEILING))
    wd = int(np.ceil(wd/2))
    n = len(X)
    for t in range(1, n+1):
        st = t - wd
        ed = t + wd
        if st < 1:
            st = 1
        if ed > n:
            ed = n
        Y = X
        counts = 0
        for i in range(ed-st+1):
            if Y[st+i-1] < th:
                counts += 1
        length = ed - st
        if counts > length/2:
            X[t-1] = 0
    return X


def FD_search(dat, params, loc, scale):
    th = 1.0
    # if starting point is too sparse, then, ignore the point
    spWD = 4
    dat = removeSparse(dat, spWD, th)
    loclist = []
    for i in range(len(dat)):
        if dat[i] > th:
            loclist.append(i)
    if not loclist:
        st = 1
    else:
        st = loclist[0]
    if st < 0:
        st = 0
    th = max(dat)
    loclist = []
    for i in range(len(dat)):
        if dat[i] == th:
            loclist.append(i)
    ed = loclist[0] + 1
    # # #
    idxlist = []
    for i in range(st, ed+1):
        idxlist.append(i)
    sselist = np.zeros(len(idxlist))
    for i in range(len(idxlist)):
        params[loc] = idxlist[i]
        sselist[i] = F_RNF(-1, dat, params, -1, scale)
    minlist = [i for i, v in enumerate(sselist) if v == min(sselist)]
    if not minlist:
        estimate = 0  # 1
    else:
        estimate = idxlist[minlist[0]]
    return estimate

##########################
# Rise and Fall fitting
##########################


def F_RNF(P, dat, params, loc, scale):
    if loc != -1:
        params = _updateP(P, loc, params)
    T = len(dat)
    b, u = spikeM.spikeM(
        T,
        params[0], params[1]/params[0], params[2], params[3],
        params[4], params[5], params[6], params[7], params[8])
    if scale == 'lin':
        pass
    elif scale == 'log':
        b = [math.log(b[i] + 1) for i in range(len(b))]
        dat = [math.log(dat[i] + 1) for i in range(len(dat))]
    elif scale == 'R5':
        b = [math.pow(b[i], 0.2) for i in range(len(b))]
        dat = [math.pow(dat[i], 0.2) for i in range(len(dat))]
    sse = np.sqrt(np.mean((np.array(b) - np.array(dat)) ** 2))
    # print('-loc:[{}]---------------------------'.format(loc))
    # print('  b:', b)
    # print('  u:', u)
    # print('sse:', sse)
    # print(np.array(b).shape, np.array(dat).shape, sse.shape)
    return sse

# # for visualization


def printRNF(dat, params):
    RSE_LIN = F_RNF(-1, dat, params, -1, 'lin')
    # # output parameters
    print('===================================')
    print('N       = ', str(params[0]))
    print('beta*N  = ', str(params[1]))
    print('slope   = ', str(params[2]))
    print('nc      = ', str(params[3]))
    print('Sc      = ', str(params[4]))
    print('bgnoise = ', str(params[5]))
    print('pcycle (Pp, Pa, Ps) = ', str(params[6]), str(params[7]), str(params[8]))
    print('-----------------------------------')
    print('error (LIN)  = ', str(RSE_LIN))
    print('===================================')
    return RSE_LIN

# # for visualization (LOG & LIN scale)


def plotsRNF(dat, params, outfn):
    T = len(dat)
    b, u = spikeM.spikeM(
        T,
        params[0], params[1]/params[0], params[2], params[3],
        params[4], params[5], params[6], params[7], params[8])
    RMSE = F_RNF(-1, dat, params, -1, 'lin')
    if not outfn:
        print('No output file')
    else:
        # # --- linear plot --- # #
        n = [i for i in range(1, len(b)+1)]
        plt.scatter(n, dat, marker='o', color='black', facecolor='None', label='Original')
        plt.plot(n, b, color='red', label=r'$\Delta$B(n)')
        plt.legend()
        plt.xlim(min(n), max(n))
        plt.xlabel('Time (n)')
        plt.ylabel('Value (lin-lin)')
        plt.title('N =' + '{:.0f}, '.format(params[0]) + r'$\beta$*N =' + '{:.2f}, '.format(params[1]) + 'RMSE = {:.2f}'.format(RMSE))
        print('save as:' + outfn)
        plt.savefig(outfn + 'LIN(py).png')
        plt.close()

        # # --- log plot --- # #
        b = [b[i] + 1 for i in range(len(b))]
        u = [u[i] + 1 for i in range(len(u))]
        dat = [dat[i] + 1 for i in range(len(dat))]
        plt.scatter(n, dat, marker='o', color='black', facecolor='None', label='Original')
        plt.plot(n, b, color='red', label=r'$\Delta$B(n)')
        plt.plot(n, u, color='lime', label='U(n)', linestyle='--')
        plt.legend()
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(params[3], len(b))
        plt.xlabel('Time (n)')
        plt.ylabel('Value (log-log)')
        plt.savefig(outfn + 'LOG(py).png')


if __name__ == "__main__":
    main()

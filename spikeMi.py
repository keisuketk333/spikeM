#!/usr/bin/env python
#!/usr/bin/python

"""
   Interactive main routine, for the SpikeM (rise and fall) equation.
   Authors: Yasuko Matsubara and Christos Faloutsos
   Date: Feb. 2012
"""

try:
    from pylab import *
    from matplotlib.widgets import Slider, Button
    # from scipy.stats import poisson
except:
    print("can not find pylab - get it from http://matplotlib.sourceforge.net/")
    print("for mac users: e.g., ")
    print("sudo port install python27")
    print("sudo port select --set python python27")
    print("sudo port install py27-matplotlib")
    print("sudo port install py27-scipy")
    print("# for GUI")
    print("mkdir ~/.matplotlib")
    print("echo backend: MacOSX > ~/.matplotlib/matplotlibrc")
    raise
try:
    import spikeM 
except:
    print("can not find spikeM.py")
    raise

import sys
import getopt
import math as M
import pylab
#
try:
    import pylab as P
    has_pylab = True
except:
    print("sudo apt-get install pylab      matplotlib ")
    print("    continuing without plotting functionality")
    has_pylab = False

import array as A

LOGLOGP = 0
ax = subplot(111)
subplots_adjust(bottom=0.35)

def RNFplot(t, yvals, uyvals, dat, plotLog):
    if plotLog == 1:
        loglogPlot(t, yvals, uyvals, dat)
    if plotLog == -1:
        semilogPlot(t, yvals, uyvals, dat)
    if plotLog == 0:
        linlinPlot(t, yvals, uyvals, dat)
    xlabel("time")
    ylabel("# of informed at time n")
    
def loglogPlot(t, yvals, uyvals, dat):
    loglog(t, yvals, '-', lw=2, color='red')
    loglog(t, uyvals, lw=2, color='blue')
    loglog(dat, '-', color='green')
    axis([0, T, 0, N_max0*10])
 
def semilogPlot(t, yvals, uyvals, dat):
    semilogy(t, yvals, '-', lw=2, color='red')
    semilogy(t, uyvals, lw=2, color='blue')
    semilogy(dat, '-', color='green')
    axis([0, T, 0, N_max0*10])
 
def linlinPlot(t, yvals, uyvals, dat):
    plot(t, yvals, '-', lw=2, color='red')
    plot(dat, '-', color='green')
    # axis([0, T, 0, max(yvals)+50])

def RSErr(T, x, y):
   if len(y) == 0:
       return 0
   val = 0.0
   for i in range(0, T):
       val += (x[i] - y[i])**2
   val = ((val/T)**0.5)
   return val

# data setting
# load data
if len(sys.argv) == 1:
    dat = []
else:
    dat = pylab.loadtxt(sys.argv[1])

if len(dat) == 0:
    # basic parameters
    T = 120
    N_max0 = 2000
    betaN0 = 1.0
    my_beta0 = betaN0/N_max0
    nc0=0
    Sc0=1
    bgn0=0
    slope0 = 1.5
    pfreq0 = 24 
    prate0 = 0.4
    pshift0 = 0
else:
    # data specific parameters
    T = len(dat)  # 24*7 # 120
    N_max0   = 7000
    betaN0   = 0.750
    my_beta0 =betaN0/N_max0
    slope0   = 1.5
    nc0      =11
    Sc0      =300
    bgn0     =0.1
    pfreq0   = 24
    prate0   = 0.47
    pshift0  = 14
    if len(dat) < T:
        T = len(dat)

# spikeM model (s: infected, us: uninfected)
(s, us) = spikeM.spikeM(
    T, N_max0, my_beta0, -slope0,
    nc0, Sc0, bgn0,
    pfreq0, prate0, pshift0)

# create the horizontal axis
t = list(range(len(s)))
print(("duration = ", len(t)))  # len(s)
t = [float(item) for item in t]
RNFplot(t, s, us, dat, LOGLOGP)
err = RSErr(T, s, dat)
ptitle = "SpikeM beta=%.2e N=%2d (RMSE=%.2f)" % (my_beta0, N_max0, err)
title(ptitle)

# GUI: slide bars
axcolor = 'lightgoldenrodyellow'
hvcolor = 'darkgoldenrod'
axbgn   = axes([0.15, 0.05, 0.25, 0.03], facecolor=axcolor)
axSc    = axes([0.15, 0.10, 0.25, 0.03], facecolor=axcolor)
axnc    = axes([0.15, 0.15, 0.25, 0.03], facecolor=axcolor)
axbeta  = axes([0.15, 0.20, 0.25, 0.03], facecolor=axcolor)
aN_max  = axes([0.15, 0.25, 0.25, 0.03], facecolor=axcolor)

apshift = axes([0.6, 0.20, 0.25, 0.03],  facecolor=axcolor)
aprate  = axes([0.6, 0.25, 0.25, 0.03],  facecolor=axcolor)
#
sbgn    = Slider(axbgn,  'epsilon', 0,  0.1 ,valinit=bgn0)
sSc     = Slider(axSc,    'S_b',      0,  300,valinit=Sc0)
snc     = Slider(axnc,    'n_b',  0,  T,   valinit=nc0)
sN_max  = Slider(aN_max,  'N',  1000,  8000, valinit=N_max0)
sbetaN  = Slider(axbeta, 'beta*N', 0.1,  1.2, valinit=my_beta0*N_max0)
sprate  = Slider(aprate,  'P_a',   0,  1,   valinit=prate0)
spshift = Slider(apshift,  'P_s',  1,  T/2,  valinit=pshift0)

def update(val):
    bgn     = sbgn.val
    Sc      = sSc.val
    nc      = snc.val
    slope   = slope0
    N_max   = sN_max.val
    betaN   = sbetaN.val
    my_beta = betaN/N_max
    pfreq   = pfreq0
    prate   = sprate.val
    pshift  = spshift.val

    # re-do the plot
    ax = subplot(111)
    # clear old plot
    cla()

    (yvals,uyvals) = spikeM.spikeM(
        T,
        N_max, my_beta, -slope,
        nc, Sc, bgn, 
        pfreq, prate, pshift
    )
    t = list(range(len(yvals)))
    t = [float(item) for item in t]
    RNFplot(t, yvals, uyvals, dat, LOGLOGP)
    err = RSErr(T, yvals, dat)
    ptitle= "SpikeM beta=%.2e N_max=%2d (RMSE=%.2f)" %\
        (my_beta, N_max, err)
    title(ptitle)
    draw()

# GUI: buttons
Logax     = axes([0.6, 0.1, 0.1, 0.04])
buttonLog = Button(Logax, 'Log/Lin', color=axcolor, hovercolor=hvcolor)
resetax   = axes([0.7, 0.1, 0.1, 0.04])
button    = Button(resetax, 'Reset', color=axcolor, hovercolor=hvcolor)

def reset(event):
    sbgn.reset()
    sSc.reset()
    snc.reset()
    sbetaN.reset()
    sN_max.reset()
    sprate.reset()
    spshift.reset()

def logscale(event):
    global LOGLOGP
   
    if LOGLOGP != 1:
        LOGLOGP += 1
    else:
        LOGLOGP=-1  # -1 # 0

# event listeners
sbgn.on_changed(update)
sSc.on_changed(update)
snc.on_changed(update)
sbetaN.on_changed(update)
sN_max.on_changed(update)
sprate.on_changed(update)
spshift.on_changed(update)
button.on_clicked(reset)
buttonLog.on_clicked(logscale)
buttonLog.on_clicked(update)

show()



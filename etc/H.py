import pandas as pd
import seaborn as sns
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import tensorflow as tf
import numpy as np
import os
import random
from google.colab import output
from Lib.lattice import *
from math import *
from scipy import optimize, integrate, stats
%cd / content/drive/My\ Drive/Tc_estimation/trial
path = '38_WL/'


def p(s):
    return path+s


def rand(a, b):
    return (b - a) * np.random.rand() + a


def notify():
    output.eval_js('new Audio("https://smar-tone.com/data/su802.mp3").play()')


def shuffle_samples(*args):
    zipped = list(zip(*args))
    np.random.shuffle(zipped)
    shuffled = list(zip(*zipped))
    result = []
    for ar in shuffled:
        result.append(np.asarray(ar))
    return result


plt.rcParams.update({'font.size': 18})
plt.tight_layout()

cmap = plt.get_cmap("tab10")

Lattice = H
Tc_R = Lattice.Tc(0.0, 0.0)
print("Tc_R = ", Tc_R)
Tc_L = 0.05

minTc, maxTc = Tc_L, Tc_R
bins = 10
ck_init = [1.0, 1.0]
centering_rate, variance = 0.9, 0.8

xlim = [minTc-0.01, maxTc+0.01]  # for plot
min_f = 1e-8

binWidth = (maxTc-minTc)/bins
binCenters = np.arange(minTc+binWidth/2, maxTc, binWidth)

sigma = [[variance, 0], [0, variance]]


def my_exp(c, n):
    return 1 if c-n > 0 else np.exp(c-n)


def isFlat(hist):
    return all(hist > np.mean(hist)*0.8)


def k_to_binInd(k1, k2):
    Tc = Lattice.Tc(k1, k2)
    if Tc < minTc or maxTc < Tc:
        return -1
    elif maxTc == Tc:
        return bins-1
    else:
        return int((Tc-minTc)//binWidth+0.5)


def mu(k1, k2):
    return [k1*centering_rate, k2*centering_rate]


def RW(ck1, ck2):
    return np.random.multivariate_normal(mu(ck1, ck2), sigma)


def g_ratio(ck1, ck2, nk1, nk2):
    denom = stats.multivariate_normal.pdf((nk1, nk2), mu(ck1, ck2), sigma)
    numer = stats.multivariate_normal.pdf((ck1, ck2), mu(nk1, nk2), sigma)
    return numer/denom


hist = np.zeros(bins)
Tc_s = np.array([])
k1s, k2s = np.array([]), np.array([])
S_s = np.zeros(bins)
ck1, ck2 = ck_init
cInd = k_to_binInd(ck1, ck2)
step = 0
f = 1.0
hist_reset_step = np.array([])


while f > min_f:
    while True:
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd = k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd != -1 and np.random.rand() <= np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

        Js = Lattice.gen_Js(ck1, ck2)
        Tc = Lattice.Tc(ck1, ck2)

        hist[cInd] += 1.0
        S_s[cInd] += f
        Tc_s = np.append(Tc_s, Tc)
        k1s = np.append(k1s, ck1)
        k2s = np.append(k2s, ck2)

        step += 1
        if isFlat(hist):
            fig, ax = plt.subplots()
            ax.bar(binCenters, hist, width=binWidth)
            ax.set_xlabel(r"$T_c$")
            ax.set_xlim(xlim)
            ax.set_title("%.1e -> %.1e" % (f, f/2))
            ax.axvline(x=minTc, c='orange')
            ax.axvline(x=maxTc, c='orange')
            plt.show()
            hist_reset_step = np.append(hist_reset_step, step)
            break

    hist = np.zeros(bins)
    f /= 2.0

notify()

fig, ax = plt.subplots(figsize=(12, 6))
ax.plot(binCenters, S_s, c=cmap(0), marker='.')
ax.set_ylabel("ln (g)")
fig.show()

if False:
    np.save(p("%s_S.npy" % (Lattice.id)), S_s)

S_s = np.load(p("%s_S.npy" % (Lattice.id)))

steps = 10000

Tc_s = np.array([])
k1s, k2s = np.array([]), np.array([])
ck1, ck2 = ck_init
cInd = k_to_binInd(ck1, ck2)

for i in range(steps):
    while True:
        try:
            nk1, nk2 = RW(ck1, ck2)
            nInd = k_to_binInd(nk1, nk2)
        except:
            pass
        else:
            break
    if nInd != -1 and np.random.rand() <= np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
        cInd, ck1, ck2 = nInd, nk1, nk2
    Js = Lattice.gen_Js(ck1, ck2)
    Tc = Lattice.Tc(ck1, ck2)

    Tc_s = np.append(Tc_s, Tc)
    k1s = np.append(k1s, ck1)
    k2s = np.append(k2s, ck2)

l = 0
r = 10000
xticks = np.arange(l, r)

fig, axes = plt.subplots(2, 1, figsize=(16, 8))
fig.suptitle(r"Trajectory of $k$ and $T_c$")

axes[0].plot(xticks, k1s[l:r], label=r"$k_1$")
axes[0].plot(xticks, k2s[l:r], label=r"$k_2$")
axes[0].set_xlabel("step")
axes[0].set_ylabel(r"$k$")
axes[0].xaxis.set_visible(False)
axes[0].legend()

axes[1].plot(np.arange(l, r), Tc_s[l:r])
axes[1].set_xlabel("Step")
axes[1].set_ylabel(r"$T_c$")
axes[1].axhline(y=maxTc, c='orange')
axes[1].axhline(y=minTc, c='orange')
axes[1].set_ylim(xlim)

fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.canvas.draw()
axpos1 = axes[0].get_position()
axpos2 = axes[1].get_position()
axes[1].set_position([axpos2.x0, axpos2.y0, axpos1.width, axpos2.height])

fig.show()

interval = 1000
samples = 10000

centi = samples//100
deci = centi*10

Js_s = np.zeros((samples, Lattice.numJ))
y = np.zeros(samples)
ck1, ck2 = ck_init
cInd = k_to_binInd(ck1, ck2)

for i in range(samples):
    for _ in range(interval):
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd = k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd != -1 and np.random.rand() <= np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

    Js = Lattice.gen_Js(ck1, ck2)
    Tc = Lattice.Tc(ck1, ck2)
    Js_s[i] = Js
    y[i] = Tc

    if i % centi == centi-1:
        print(i//centi+1 if i % deci == deci-1 else ".", end='')
notify()

fig, ax = plt.subplots()
ax.hist(y, bins=20, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

Js_s, y = shuffle_samples(Js_s, y)

if False:
    np.savez(p("%s_Js_Tc.npz" % (Lattice.id)), Js_s, y)

Js_s, y = np.load(p("%s_Js_Tc.npz" % (Lattice.id))).values()

fig, ax = plt.subplots()
ax.hist(y, bins=30, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

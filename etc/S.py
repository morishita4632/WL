import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import lattice
from math import *
from scipy import optimize, integrate
path = '38_WL/'


def p(s):
    return path+s


def rand(a, b):
    return (b - a) * np.random.rand() + a


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

S.build()

# S
$J_1 = 1$は固定。$J_2 = 10 ^ k$とし、$k$を与える（その後規格化する）。

ランダムウォークでは、$k$に一様乱数（[-1, 1]など）を加える。

flat histgram の判定は「全てのビンが平均値の80 % 以上である」とした（原論文）。

print("max Tc = ", S.Tc(0.0))

minTc, maxTc = 0.1, S.Tc(0.0)
xlim = [0.0, 0.6]
bins = 30
RWstep = 1.0
min_f = 1e-8

binWidth = (maxTc-minTc)/bins
binCenters = np.arange(minTc+binWidth/2, maxTc, binWidth)


def my_exp(c, n):
    return 1 if c-n > 0 else np.exp(c-n)


def isFlat(hist):
    return all(hist > np.mean(hist)*0.8)


def k_to_binInd(k):
    Tc = S.Tc(k)
    if Tc < minTc or maxTc <= Tc:
        return -1
    else:
        return int((Tc-minTc)//binWidth+0.5)


hist = np.zeros(bins)
Tc_s = np.array([])
S_s = np.zeros(bins)
ck = 0.0
cInd = k_to_binInd(ck)
repeat = 0
f = 1.0

while f > min_f:
    while True:
        while True:
            try:
                dk = rand(-RWstep, RWstep)
                nk = ck+dk
                nInd = k_to_binInd(nk)
            except:
                pass
            else:
                break
        if nInd != -1 and np.random.rand() <= my_exp(S_s[cInd], S_s[nInd]):
            cInd = nInd
            ck = nk

        Js = S.gen_Js(ck)
        Tc = S.Tc(ck)

        hist[cInd] += 1.0
        S_s[cInd] += f
        Tc_s = np.append(Tc_s, Tc)

        repeat += 1
        if isFlat(hist):
            fig, ax = plt.subplots()
            ax.bar(binCenters, hist, width=binWidth)
            ax.set_xlabel(r"$T_c$")
            ax.set_xlim(xlim)
            ax.set_title("%.1e -> %.1e" % (f, f/2))
            ax.axvline(x=minTc, c='orange')
            ax.axvline(x=maxTc, c='orange')
            plt.show()
            break

    hist = np.zeros(bins)
    f /= 2.0

notify()

l = 000
r = l+1000
fig, ax = plt.subplots()
ax.set_title("Trajectory of " + r"$T_c$")
ax.plot(np.arange(l, r), Tc_s[l:r])
ax.set_xlabel("Step")
ax.set_ylabel("Tc")
ax.axhline(y=minTc, c='orange')
ax.axhline(y=maxTc, c='orange')
fig.show()

fig, ax = plt.subplots()
ax.plot(binCenters, S_s, marker='.')
ax.set_xlabel(r"$T_c$")
ax.set_xlim(xlim)
ax.set_ylabel("ln(g)")
fig.show()

if False:
    np.save(p("%s_S.npy" % (S.id)), S_s)

S_s = np.load(p("%s_S.npy" % (S.id)))

interval = 1000
samples = 10000

centi = samples//100
deci = centi*10

Js_s = np.zeros((samples, 2))
y = np.zeros(samples)
ck = 0.0
cInd = k_to_binInd(ck)

for i in range(samples):
    for _ in range(interval):
        while True:
            try:
                dk = rand(-RWstep, RWstep)
                nk = ck+dk
                nInd = k_to_binInd(nk)
            except:
                pass
            else:
                break
        if nInd != -1 and np.random.rand() <= my_exp(S_s[cInd], S_s[nInd]):
            ck = nk
            cInd = nInd

    Js = S.gen_Js(ck)
    Tc = S.Tc(ck)
    Js_s[i] = Js
    y[i] = Tc

    if i % centi == centi-1:
        print(i//centi+1 if i % deci == deci-1 else ".", end='')
notify()

fig, ax = plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

Js_s, y = shuffle_samples(Js_s, y)

if True:
    np.savez(p("%s_Js_Tc.npz" % (S.id)), Js_s, y)

Js_s, y = np.load(p("%s_Js_Tc.npz" % (S.id))).values()

fig, ax = plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

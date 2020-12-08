%cd /content/drive/My\ Drive/Tc_estimation/trial
path = '38_WL/'

import random, os, numpy as np, tensorflow as tf, matplotlib.pyplot as plt, matplotlib.ticker as ticker, seaborn as sns, pandas as pd
from scipy import optimize, integrate, stats
from math import *
from Lib.lattice import *
from google.colab import output

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

Lattice = T
Tc_R=Lattice.Tc(0.0, 0.0)
print("Tc_R = ", Tc_R)

Tc_M=S.Tc_Js(np.ones(2)/4)
print("Tc_M = ", Tc_M)

Tc_L=0.1
high_range_width=Tc_R-Tc_M
low_range_width=Tc_M-Tc_L

samples_all = 10000
samples_high =int(samples_all * high_range_width/(high_range_width+low_range_width))
samples_low = samples_all-samples_high
print("#samples...(%d, %d)"%(samples_high, samples_low))

# 上

minTc, maxTc = Tc_M, Tc_R
bins=10
RWstep=0.5
ck_init = [1.2, 1.2]

xlim=[minTc-0.01, maxTc+0.01] # for plot
min_f=1e-8

binWidth=(maxTc-minTc)/bins
binCenters=np.arange(minTc+binWidth/2, maxTc, binWidth)

# sigma = [[variance, 0], [0, variance]]

def my_exp(c, n):
    return 1 if c-n>0 else np.exp(c-n)

def isFlat(hist):
    return all(hist>np.mean(hist)*0.8)

def k_to_binInd(k1, k2):
    Tc=Lattice.Tc(k1, k2)
    if Tc<minTc or maxTc<Tc:
        return -1
    elif maxTc==Tc:
        return bins-1
    else:
        return int((Tc-minTc)//binWidth+0.5)

# def mu(k1, k2):
#     return [k1*centering_rate, k2*centering_rate]

def RW(ck1, ck2):
    # return np.random.multivariate_normal(mu(ck1, ck2), sigma)
    return [ck1+rand(-RWstep, RWstep), ck2+rand(-RWstep, RWstep)]

def g_ratio(ck1, ck2, nk1, nk2):
    return 1.0
    # denom=stats.multivariate_normal.pdf((nk1, nk2), mu(ck1, ck2), sigma)
    # numer=stats.multivariate_normal.pdf((ck1, ck2), mu(nk1, nk2), sigma)
    # return numer/denom

hist = np.zeros(bins)
Tc_s=np.array([])
k1s, k2s = np.array([]), np.array([])
S_s = np.zeros(bins)
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)
step=0
f=1.0
hist_reset_step=np.array([])



while f > min_f:
    while True:
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd=k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

        Js = Lattice.gen_Js(ck1, ck2)
        Tc=Lattice.Tc(ck1, ck2)

        hist[cInd]+=1.0
        S_s[cInd]+=f
        Tc_s=np.append(Tc_s, Tc)
        k1s=np.append(k1s, ck1)
        k2s=np.append(k2s, ck2)

        step+=1
        if isFlat(hist):
            fig, ax=plt.subplots()
            ax.bar(binCenters, hist, width=binWidth)
            ax.set_xlabel(r"$T_c$")
            ax.set_xlim(xlim)
            ax.set_title("%.1e -> %.1e" % (f, f/2) )
            ax.axvline(x=minTc, c='orange')
            ax.axvline(x=maxTc, c='orange')
            plt.show()
            hist_reset_step=np.append(hist_reset_step, step)
            break

    hist=np.zeros(bins)
    f/=2.0

notify()

fig, ax=plt.subplots(figsize=(12,6))
ax.plot(binCenters, S_s, c=cmap(0), marker='.')
ax.set_ylabel("ln (g)")
fig.show()

if False:
    np.save(p("%s_S_high.npy" % (Lattice.id)), S_s)

S_s = np.load(p("%s_S_high.npy" % (Lattice.id)))

steps=10000

Tc_s=np.array([])
k1s, k2s = np.array([]), np.array([])
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)

for i in range(steps):
    while True:
        try:
            nk1, nk2 = RW(ck1, ck2)
            nInd=k_to_binInd(nk1, nk2)
        except:
            pass
        else:
            break
    if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
        cInd, ck1, ck2 = nInd, nk1, nk2
    Js = Lattice.gen_Js(ck1, ck2)
    Tc=Lattice.Tc(ck1, ck2)

    Tc_s=np.append(Tc_s, Tc)
    k1s=np.append(k1s, ck1)
    k2s=np.append(k2s, ck2)

l=0
r=10000
xticks=np.arange(l, r)

fig,axes=plt.subplots(2, 1, figsize=(16,8))
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
axes[1].axhline(y=S_maxTc, c='orange')
axes[1].set_ylim(xlim)

fig.tight_layout(rect=[0,0,1,0.96])
fig.canvas.draw()
axpos1 = axes[0].get_position()
axpos2 = axes[1].get_position()
axes[1].set_position([axpos2.x0, axpos2.y0, axpos1.width, axpos2.height])

fig.show()

interval=1000
samples=samples_high

centi=samples//100
deci=centi*10

Js_s=np.zeros((samples, Lattice.numJ))
y=np.zeros(samples)
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)

for i in range(samples):
    for _ in range(interval):
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd=k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

    Js = Lattice.gen_Js(ck1, ck2)
    Tc=Lattice.Tc(ck1, ck2)
    Js_s[i] = Js
    y[i]=Tc

    if i%centi==centi-1:
        print(i//centi+1 if i%deci==deci-1 else ".", end='')
notify()

fig,ax=plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

Js_s, y = shuffle_samples(Js_s, y)

if False:
    np.savez(p("%s_Js_Tc_high.npz" % (Lattice.id)), Js_s, y)

Js_s, y=np.load(p("%s_Js_Tc_high.npz" % (Lattice.id))).values()

fig,ax=plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

# 下

minTc, maxTc = Tc_L, Tc_M
bins=10
ck_init = [1.0, 3.0]
centering_rate, variance = 0.9, 0.8

xlim=[minTc-0.01, maxTc+0.01] # for plot
min_f=1e-8

binWidth=(maxTc-minTc)/bins
binCenters=np.arange(minTc+binWidth/2, maxTc, binWidth)

sigma = [[variance, 0], [0, variance]]

def my_exp(c, n):
    return 1 if c-n>0 else np.exp(c-n)

def isFlat(hist):
    return all(hist>np.mean(hist)*0.8)

def k_to_binInd(k1, k2):
    Tc=Lattice.Tc(k1, k2)
    if Tc<minTc or maxTc<Tc:
        return -1
    elif maxTc==Tc:
        return bins-1
    else:
        return int((Tc-minTc)//binWidth+0.5)

def mu(k1, k2):
    return [k1*centering_rate, k2*centering_rate]

def RW(ck1, ck2):
    return np.random.multivariate_normal(mu(ck1, ck2), sigma)

def g_ratio(ck1, ck2, nk1, nk2):
    denom=stats.multivariate_normal.pdf((nk1, nk2), mu(ck1, ck2), sigma)
    numer=stats.multivariate_normal.pdf((ck1, ck2), mu(nk1, nk2), sigma)
    return numer/denom

hist = np.zeros(bins)
Tc_s=np.array([])
k1s, k2s = np.array([]), np.array([])
S_s = np.zeros(bins)
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)
step=0
f=1.0
hist_reset_step=np.array([])



while f > min_f:
    while True:
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd=k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

        Js = Lattice.gen_Js(ck1, ck2)
        Tc=Lattice.Tc(ck1, ck2)

        hist[cInd]+=1.0
        S_s[cInd]+=f
        Tc_s=np.append(Tc_s, Tc)
        k1s=np.append(k1s, ck1)
        k2s=np.append(k2s, ck2)

        step+=1
        if isFlat(hist):
            fig, ax=plt.subplots()
            ax.bar(binCenters, hist, width=binWidth)
            ax.set_xlabel(r"$T_c$")
            ax.set_xlim(xlim)
            ax.set_title("%.1e -> %.1e" % (f, f/2) )
            ax.axvline(x=minTc, c='orange')
            ax.axvline(x=maxTc, c='orange')
            plt.show()
            hist_reset_step=np.append(hist_reset_step, step)
            break

    hist=np.zeros(bins)
    f/=2.0

notify()

fig, ax=plt.subplots(figsize=(12,6))
ax.bar(binCenters, hist, width=binWidth)
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title("#bin = %d, RW params = (%.1f, %.1f)" % (bins, centering_rate, variance))
ax.tick_params(axis='y', colors=cmap(0))

ax2=ax.twinx()
ax2.plot(binCenters, S_s, c=cmap(1), marker='.')
ax2.tick_params(colors=cmap(1))
ax2.set_ylabel("ln (g)")
fig.show()

r=len(k2s)
try:
    l=hist_reset_step[-1].astype(int)
except:
    l=0

xticks=np.arange(l, r)

fig,axes=plt.subplots(2, 1, figsize=(16,8))
fig.suptitle(r"Trajectory of $k$ and $T_c$")

axes[0].plot(xticks, k1s[l:r], label=r"$k_1$")
axes[0].plot(xticks, k2s[l:r], label=r"$k_2$")
axes[0].set_xlabel("step")
axes[0].set_ylabel(r"$k$")
axes[0].xaxis.set_visible(False)
axes[0].set_xlim([l, r])
for x in hist_reset_step:
    axes[0].axvline(x=x, c='red')
axes[0].legend()

axes[1].plot(np.arange(l, r), Tc_s[l:r])
axes[1].set_xlabel("Step")
axes[1].set_ylabel(r"$T_c$")
axes[1].set_xlim([l, r])
axes[1].axhline(y=minTc, c='orange')
axes[1].axhline(y=maxTc, c='orange')
axes[1].axhline(y=maxTc-binWidth, c='orange')
for x in hist_reset_step:
    axes[1].axvline(x=x, c='red')

fig.tight_layout(rect=[0,0,1,0.96])
fig.canvas.draw()
axpos1 = axes[0].get_position()
axpos2 = axes[1].get_position()
axes[1].set_position([axpos2.x0, axpos2.y0, axpos1.width, axpos2.height])

fig.show()

fig, ax=plt.subplots(figsize=(12,6))
ax.plot(binCenters, S_s, c=cmap(0), marker='.')
ax.set_ylabel("ln (g)")
fig.show()

if True:
    np.save(p("%s_S_low.npy" % (Lattice.id)), S_s)

S_s = np.load(p("%s_S_low.npy" % (Lattice.id)))

steps=10000

Tc_s=np.array([])
k1s, k2s = np.array([]), np.array([])
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)

for i in range(steps):
    while True:
        try:
            nk1, nk2 = RW(ck1, ck2)
            nInd=k_to_binInd(nk1, nk2)
        except:
            pass
        else:
            break
    if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
        cInd, ck1, ck2 = nInd, nk1, nk2
    Js = Lattice.gen_Js(ck1, ck2)
    Tc=Lattice.Tc(ck1, ck2)

    Tc_s=np.append(Tc_s, Tc)
    k1s=np.append(k1s, ck1)
    k2s=np.append(k2s, ck2)

l=0
r=10000
xticks=np.arange(l, r)

fig,axes=plt.subplots(2, 1, figsize=(16,8))
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
axes[1].axhline(y=S_maxTc, c='orange')
axes[1].set_ylim(xlim)

fig.tight_layout(rect=[0,0,1,0.96])
fig.canvas.draw()
axpos1 = axes[0].get_position()
axpos2 = axes[1].get_position()
axes[1].set_position([axpos2.x0, axpos2.y0, axpos1.width, axpos2.height])

fig.show()

interval=1000
samples=samples_low

centi=samples//100
deci=centi*10

Js_s=np.zeros((samples, Lattice.numJ))
y=np.zeros(samples)
ck1, ck2 = ck_init
cInd=k_to_binInd(ck1, ck2)

for i in range(samples):
    for _ in range(interval):
        while True:
            try:
                nk1, nk2 = RW(ck1, ck2)
                nInd=k_to_binInd(nk1, nk2)
            except:
                pass
            else:
                break
        if nInd!=-1 and np.random.rand()<=np.exp(S_s[cInd]-S_s[nInd])*g_ratio(ck1, ck2, nk1, nk2):
            cInd, ck1, ck2 = nInd, nk1, nk2

    Js = Lattice.gen_Js(ck1, ck2)
    Tc=Lattice.Tc(ck1, ck2)
    Js_s[i] = Js
    y[i]=Tc

    if i%centi==centi-1:
        print(i//centi+1 if i%deci==deci-1 else ".", end='')
notify()

fig,ax=plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

Js_s, y = shuffle_samples(Js_s, y)

if False:
    np.savez(p("%s_Js_Tc_low.npz" % (Lattice.id)), Js_s, y)

Js_s, y=np.load(p("%s_Js_Tc_low.npz" % (Lattice.id))).values()

fig,ax=plt.subplots()
ax.hist(y, bins=bins, range=(minTc, maxTc))
ax.set_xlim(xlim)
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=minTc, c='orange')
ax.axvline(x=maxTc, c='orange')
fig.show()

# 全体

Js_s_low, y_low=np.load(p("%s_Js_Tc_low.npz" % (Lattice.id))).values()
Js_s_high, y_high=np.load(p("%s_Js_Tc_high.npz" % (Lattice.id))).values()
Js_s = np.vstack([Js_s_low, Js_s_high])
y = np.append(y_low, y_high)
Js_s, y = shuffle_samples(Js_s, y)

fig,ax=plt.subplots()
ax.hist(y, bins=50, range=(Tc_L, Tc_R))
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=Tc_L, c='orange')
ax.axvline(x=Tc_R, c='orange')
fig.show()

if False:
    np.savez(p("%s_Js_Tc.npz"% (Lattice.id)), Js_s, y)

Js_s, y = np.load(p("%s_Js_Tc.npz"% (Lattice.id))).values()

fig,ax=plt.subplots()
ax.hist(y, bins=50, range=(Tc_L, Tc_R))
ax.set_xlabel(r"$T_c$")
ax.set_title(r"$T_c$" + " distribution")
ax.axvline(x=Tc_L, c='orange')
ax.axvline(x=Tc_R, c='orange')
fig.show()
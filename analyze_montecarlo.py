import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from scipy.stats import norm

# 0 |   .   | 0.5   1    0
# 1 | .   . | 0.25  0.5  1 2
# 2 |. . . .| 0.125 0.25 3 4 5 6
# 3                      7

def show(fig):
    if False:
        fig.show()
    else:
        fig.suptitle(fig.get_label())
        fig.savefig("plots4/{}.png".format(fig.get_label().lower().replace(' ','_')))

def hier_ix(ix, shift_n1p1=True):
    level = int(np.log2(ix+1))
    offset = 0.5**(level + 1)
    stride = 0.5**level
    locix = ix - (2**level - 1)

    if shift_n1p1:
        return (offset + locix*stride)*2 - 1
    else:
        return offset + locix*stride

hier_x_centered = np.array([-1,1] + [hier_ix(ix,shift_n1p1=True) for ix in range(1500-2)])
s = 1000.0/6
hier_x_centered *= 3*s
hier_p = 1.0/np.sqrt(2*np.pi*s**2) * np.exp(-(hier_x_centered)**2 / (2 * s**2))

ix_u0 = 2
ix_v0 = 3
ix_u1 = 4
ix_v1 = 5
ix_u2 = 6
ix_v2 = 7
ix_vel0 = 8
ix_vel1 = 9
ix_vel2 = 10

all_dat = None

N = 2000
for i in range(N):
    if i % 100 == 0: print(i)
    dat = np.loadtxt('mc50/{}.dat'.format(i))
    if all_dat is None:
        all_dat = [dat]
    else:
        assert((dat[:,0] == all_dat[0][:,0]).all())
        assert((dat[:,1] == all_dat[0][:,1]).all())
        all_dat.append(dat)


all_dat = np.array(all_dat)

print("{} measurements x {} timesteps x {} channels".format(*all_dat.shape))

mean = np.mean(all_dat, axis=0)
sigma  = np.std(all_dat, axis=0)

last_time_step = all_dat[:,-1,:] # [measurement, ix]
# mean and simga for increasing N
running_mean = np.cumsum(last_time_step,axis=0) / np.arange(1,1+last_time_step.shape[0])[:,np.newaxis] # [measurements, ix]
running_sigma  = np.zeros_like(running_mean) # [measurements, ix]
for n in range(all_dat.shape[0]):
    diffs = last_time_step[:n,:] - running_mean[n,:]
    running_sigma[n,:] = np.sqrt(np.sum(diffs*diffs, axis=0) / (n - 1))

ix_show = ix_vel0

for ix_show in range(2,11):
    print("last mean: {:>.5e}, last sigma: {:>.5e}".format(mean[-1,ix_show], sigma[-1,ix_show]))

    fig = plt.figure("Monte Carlo over Time for {}".format(ix_show), figsize=(20,12))
    fig.clear()
    ax = fig.add_axes([0.1,0.1,0.65,0.65])
    ax.plot(mean[:,1], mean[:,ix_show], 'r', label=r"mean $\mu(t)$, $\mu(t)\pm\sigma(t)$")
    ax.fill_between(mean[:,1], mean[:,ix_show] + sigma[:,ix_show], mean[:,ix_show] - sigma[:,ix_show], color='r', alpha=.25)
    ax.legend(loc='best')
    ax.set_xlabel("t")
    axh = fig.add_axes([0.75+0.02,0.1,0.2,0.65])
    axh.yaxis.set_major_formatter(NullFormatter())
    axh.xaxis.set_major_formatter(NullFormatter())
    axh.hist(all_dat[:,-1,ix_show], color='g', bins=np.sqrt(N), alpha=0.5, normed=True, orientation="horizontal", label="histogram for t={:1.1f}s".format(all_dat[0,-1,1]))
    axv = fig.add_axes([0.1,0.75+0.02,0.65,0.17])
    axv.xaxis.set_major_formatter(NullFormatter())
    axv.plot(mean[:,1], sigma[:,ix_show]**2, color='r', label=r"variance $\sigma^2(t)$")
    axv.legend(loc='best')

    ylimits = (min(ax.get_ylim()[0], axh.get_ylim()[0]), max(ax.get_ylim()[1], axh.get_ylim()[1]))
    ax.set_ylim(ylimits)
    axh.set_ylim(ylimits)

    x = np.linspace(ylimits[0], ylimits[1], 1000)
    axh.plot(norm.pdf(x, mean[-1,ix_show], sigma[-1,ix_show]), x, 'k', linewidth=2, label=r"fit $\mu={:>.1E}, \sigma={:>.1E}$".format(mean[-1,ix_show], sigma[-1,ix_show]))
    axh.legend(loc='best')

    ax.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
    axv.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
    show(fig)

# Convergence of MC

fig = plt.figure("Monte Carlo Convergence - Mean")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.loglog(np.abs(running_mean[:,2:] - running_mean[-1,2:]), basex=2, basey=2)
ax.loglog(1/np.sqrt(np.arange(running_mean.shape[0])) / 2**6, basex=2, basey=2)
ax.set_xlabel("N")
show(fig)

fig = plt.figure("Monte Carlo Running Mean")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.plot(running_mean[:,2:])
ax.set_xlabel("N")
show(fig)

fig = plt.figure("Monte Carlo Convergence - Std")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.loglog(np.abs(running_sigma[:,2:] - running_sigma[-1,2:]), basex=2, basey=2)
ax.loglog(1/np.sqrt(np.arange(running_mean.shape[0])) / 2**6, basex=2, basey=2)
ax.set_xlabel("N")
show(fig)

fig = plt.figure("Monte Carlo Running Std")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.plot(running_sigma[:,2:])
ax.set_xlabel("N")
show(fig)

# Convergence of Trapezoidal Rule

all_tr_dat = None

N = 1500
for i in range(N):
    if i % 100 == 0: print(i)
    dat = np.loadtxt('tr50/{}.dat'.format(i))
    if all_tr_dat is None:
        all_tr_dat = [dat]
    else:
        assert((dat[:,0] == all_tr_dat[0][:,0]).all())
        assert((dat[:,1] == all_tr_dat[0][:,1]).all())
        all_tr_dat.append(dat)

all_tr_dat = np.array(all_tr_dat) # [measurement, time, ix]

trapezoidal = []
trapezoidal_sigma = []
X_p = hier_p[:,np.newaxis,np.newaxis] * all_tr_dat
#X_p[[0,1]] /= 2
#     +/- 3 sigma   * (a+b) / 2
trapezoidal.append(2*3*(1000/6) * (X_p[0,:,:] + X_p[1,:,:])/2)
trapezoidal_sigma.append(2*3*(1000/6) * np.sum((all_tr_dat[:2,:,:] - trapezoidal[-1]) ** 2 * hier_p[:2,np.newaxis,np.newaxis], axis=0) / 2)
start = 2
num = 1
level = 1
while True:
    trapezoidal.append(2*3*(1000/6) / 2**level * np.sum(X_p[start:start+num,:,:], axis=0) + trapezoidal[-1] / 2)
    trapezoidal_sigma.append(2*3*(1000/6) / 2**level * np.sum((all_tr_dat[start:start+num,:,:] - trapezoidal[-1]) ** 2  * hier_p[start:start+num,np.newaxis,np.newaxis], axis=0) + trapezoidal_sigma[-1] / 2)
    trapezoidal_sigma[-1] = trapezoidal_sigma[-1]
    level += 1
    start += num
    num *= 2
    if start+num > X_p.shape[0]: break
print("used first {} data points for trapezoidal rule".format(start))


# t(h)   = I + ch²   + ch⁴
# t(h/2) = I + ch²/4 + ch⁴/16
# s(h)   = I +         ch⁴(1/4 - 1)/(4-1)
#        = I +         ch⁴(-1/4)
# s(h)   = I +         c'h⁴
# s(h/2) = I +         c'h⁴ / 16                                # convergence for last time step:
trapezoidal = np.array(trapezoidal)                             # good h^2 convergence after the first two steps
simpson  = (4*trapezoidal[1:,:] - 1*trapezoidal[:-1,:]) / (4-1) # good h^4 convergence after the first two steps, levels out at 3.7e-9
simpson2 = (16*simpson[1:,:] - 1*simpson[:-1,:]) / (16-1)       # good h^6 convergence after the first two steps, quickly levels out at 3.7e-9
simpson3 = (64*simpson[1:,:] - 1*simpson[:-1,:]) / (64-1)       # only h^6 convergence after the first two steps, quickly levels out at 3.7e-9

trapezoidal_sigma = np.sqrt(np.array(trapezoidal_sigma)) # variance = sigma^2

for ix_show in range(2,11):
    fig = plt.figure("Monte Carlo vs Trapezoidal over Time for {}".format(ix_show), figsize=(20,12))
    fig.clear()
    ax = fig.add_axes([0.1,0.1,0.65,0.65])
    ax.plot(mean[:,1], mean[:,ix_show], 'r', label=r"mean mc $\mu(t)$, $\mu(t)\pm\sigma(t)$")
    ax.fill_between(mean[:,1], mean[:,ix_show] + sigma[:,ix_show], mean[:,ix_show] - sigma[:,ix_show], color='r', alpha=.25)
    ax.plot(mean[:,1], trapezoidal[-1,:,ix_show], 'b', label=r"mean tr $\mu(t)$, $\mu(t)\pm\sigma(t)$") #use most refined trapezoidal
    ax.fill_between(mean[:,1], trapezoidal[-1,:,ix_show] + trapezoidal_sigma[-1,:,ix_show], trapezoidal[-1,:,ix_show] - trapezoidal_sigma[-1,:,ix_show], color='b', alpha=.25)
    ax.legend(loc='best')
    ax.set_xlabel("t")
    axh = fig.add_axes([0.75+0.02,0.1,0.2,0.65])
    axh.yaxis.set_major_formatter(NullFormatter())
    axh.xaxis.set_major_formatter(NullFormatter())
    axh.hist(all_dat[:,-1,ix_show], color='r', bins=np.sqrt(N), alpha=0.25, normed=True, orientation="horizontal", label="histogram for t={:1.1f}s".format(all_dat[0,-1,1]))
    axv = fig.add_axes([0.1,0.75+0.02,0.65,0.17])
    axv.xaxis.set_major_formatter(NullFormatter())
    axv.plot(mean[:,1], sigma[:,ix_show]**2, color='r', label=r"variance mc $\sigma^2(t)$")
    axv.plot(mean[:,1], trapezoidal_sigma[-1,:,ix_show]**2, color='b', label=r"variance tr $\sigma^2(t)$")
    axv.legend(loc='best')

    ylimits = (min(ax.get_ylim()[0], axh.get_ylim()[0]), max(ax.get_ylim()[1], axh.get_ylim()[1]))
    ax.set_ylim(ylimits)
    axh.set_ylim(ylimits)

    x = np.linspace(ylimits[0], ylimits[1], 1000)
    axh.plot(norm.pdf(x, mean[-1,ix_show], sigma[-1,ix_show]), x, 'r', linewidth=2, label=r"fit $\mu={:>.1E}, \sigma={:>.1E}$".format(mean[-1,ix_show], sigma[-1,ix_show]))
    axh.plot(norm.pdf(x, trapezoidal[-1,-1,ix_show], trapezoidal_sigma[-1,-1,ix_show]), x, 'b', linewidth=2, label=r"fit $\mu={:>.1E}, \sigma={:>.1E}$".format(trapezoidal[-1,-1,ix_show], trapezoidal_sigma[-1,-1,ix_show]))
    axh.legend(loc='best')

    ax.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
    axv.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
    show(fig)


fig = plt.figure("Trapezoidal")
fig.clear()
ax = fig.add_subplot(111)
ax.set_xlabel("levels")
ax.set_ylabel("error")
ax.set_xticks(np.arange(trapezoidal.shape[0]))
ax.grid()
ax.semilogy(np.arange(trapezoidal.shape[0]), np.abs(trapezoidal[:,-1,2:] - simpson3[-1,-1,2:]), 'r-', basey=2)
ax.semilogy(np.arange(simpson.shape[0]), np.abs(simpson[:,-1,2:] - simpson3[-1,-1,2:]), 'b-', basey=2)
ax.semilogy(np.arange(simpson2.shape[0]), np.abs(simpson2[:,-1,2:] - simpson3[-1,-1,2:]), 'g-', basey=2)
ax.semilogy(np.arange(simpson3.shape[0]), np.abs(simpson3[:,-1,2:] - simpson3[-1,-1,2:]), 'orange', basey=2)
ax.legend(loc='best')
show(fig)

fig = plt.figure("Trapezoidal Error Difference")
fig.clear()
ax = fig.add_subplot(111)
ax.set_xlabel("levels")
ax.set_ylabel("error")
ax.set_xticks(np.arange(trapezoidal.shape[0]))
ax.grid()
ax.semilogy(np.arange(trapezoidal.shape[0]-1), np.abs(trapezoidal[1:,-1,ix_u1] - trapezoidal[:-1,-1,ix_u1]), 'r-', basey=2)
ax.semilogy(np.arange(simpson.shape[0]-1),     np.abs(simpson[1:,-1,ix_u1] - simpson[:-1,-1,ix_u1]), 'b-', basey=2)
ax.semilogy(np.arange(simpson2.shape[0]-1),    np.abs(simpson2[1:,-1,ix_u1] - simpson2[:-1,-1,ix_u1]), 'g-', basey=2)
ax.semilogy(np.arange(simpson3.shape[0]-1),    np.abs(simpson3[1:,-1,ix_u1] - simpson3[:-1,-1,ix_u1]), 'orange', basey=2)
ax.legend(loc='best')
show(fig)

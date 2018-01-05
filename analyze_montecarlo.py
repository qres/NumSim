import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

# 0 |   .   | 0.5   1    0
# 1 | .   . | 0.25  0.5  1 2
# 2 |. . . .| 0.125 0.25 3 4 5 6
# 3                      7

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

print(all_dat.shape)

mean = np.mean(all_dat, axis=0)
sigma  = np.std(all_dat, axis=0)

last_time_step = all_dat[:,-1,:] # [measurement, ix]
running_mean = np.cumsum(last_time_step,axis=0) / np.arange(1,1+last_time_step.shape[0])[:,np.newaxis] # [measurements, ix]
running_sigma  = np.zeros_like(running_mean) # [measurements, ix]
for n in range(all_dat.shape[0]):
    diffs = last_time_step[:n,:] - running_mean[n,:]
    running_sigma[n,:] = np.sqrt(np.sum(diffs*diffs, axis=0) / (n - 1))

ix_show = ix_vel0

fig = plt.figure("Monte Carlo over Time")
fig.clear()
ax = fig.add_axes([0.1,0.1,0.65,0.65])
ax.plot(mean[:,1], mean[:,ix_show], 'r', label=r"mean $\mu$, $\mu\pm\sigma$")
ax.fill_between(mean[:,1], mean[:,ix_show] + sigma[:,ix_show], mean[:,ix_show] - sigma[:,ix_show], color='r', alpha=.25)
ax.legend(loc='best')
axh = fig.add_axes([0.75+0.02,0.1,0.2,0.65])
axh.yaxis.set_major_formatter(NullFormatter())
axh.xaxis.set_major_formatter(NullFormatter())
axh.hist(all_dat[:,-1,ix_show], bins=np.sqrt(N), orientation="horizontal")
axv = fig.add_axes([0.1,0.75+0.02,0.65,0.2])
axv.xaxis.set_major_formatter(NullFormatter())
axv.plot(mean[:,1], sigma[:,ix_show]**2, color='r', label=r"variance $\sigma^2$")
axv.legend(loc='best')

print("last mean: {}, last sigma: {}".format(mean[-1,ix_show], sigma[-1,ix_show]))

ylimits = (min(ax.get_ylim()[0], axh.get_ylim()[0]), max(ax.get_ylim()[1], axh.get_ylim()[1]))
ax.set_ylim(ylimits)
axh.set_ylim(ylimits)

ax.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
axv.set_xlim(0,all_dat[0,-1,1]) # tight to last time stamp
fig.show()

# Convergence of MC

fig = plt.figure("Monte Carlo Convergence -- Mean")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.loglog(np.abs(running_mean[:,2:] - mean[-1,2:]), basex=2, basey=2)
ax.loglog(np.arange(running_mean.shape[0])**-np.log2(np.sqrt(2)) / 2**6, basex=2, basey=2)
fig.show()

fig = plt.figure("Monte Carlo Convergence -- Std")
fig.clear()
ax = fig.add_subplot(111)
ax.grid()
ax.loglog(np.abs(running_sigma[:,2:] - sigma[-1,2:]), basex=2, basey=2)
ax.loglog(np.arange(running_mean.shape[0])**-np.log2(np.sqrt(2)) / 2**6, basex=2, basey=2)
fig.show()

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

all_tr_dat = np.array(all_tr_dat)
last_tr_time_step = all_tr_dat[:,-1,:] # [measurement, ix]

trapezoidal = []
X_p = hier_p[:,np.newaxis] * last_tr_time_step
#X_p[[0,1]] /= 2
#     +/- 3 sigma   * (a+b) / 2
trapezoidal.append(2*3*(1000/6) * (X_p[0,:] + X_p[1,:])/2)
start = 2
num = 1
level = 1
while True:
    trapezoidal.append(2*3*(1000/6) / 2**level * np.sum(X_p[start:start+num,:], axis=0) + trapezoidal[-1] / 2)
    print(level, start, start+num, num)
    level += 1
    start += num
    num *= 2
    if start+num > X_p.shape[0]: break


# t(h)   = I + ch²   + ch⁴
# t(h/2) = I + ch²/4 + ch⁴/16
# s(h)   = I +         ch⁴(1/4 - 1)/(4-1)
#        = I +         ch⁴(-1/4)
# s(h)   = I +         c'h⁴
# s(h/2) = I +         c'h⁴ / 16
trapezoidal = np.array(trapezoidal)                             # good h^2 convergence after the first two steps
simpson  = (4*trapezoidal[1:,:] - 1*trapezoidal[:-1,:]) / (4-1) # good h^4 convergence after the first two steps, levels out at 3.7e-9
simpson2 = (16*simpson[1:,:] - 1*simpson[:-1,:]) / (16-1)       # good h^6 convergence after the first two steps, quickly levels out at 3.7e-9
simpson3 = (64*simpson[1:,:] - 1*simpson[:-1,:]) / (64-1)       # only h^6 convergence after the first two steps, quickly levels out at 3.7e-9

fig = plt.figure("Trapezoidal")
fig.clear()
ax = fig.add_subplot(111)
ax.set_xlabel("levels")
ax.set_ylabel("error")
ax.set_xticks(np.arange(trapezoidal.shape[0]))
ax.grid()
ax.semilogy(np.arange(trapezoidal.shape[0]), np.abs(trapezoidal[:,2:] - simpson3[-1,2:]), 'r-', basey=2)
ax.semilogy(np.arange(simpson.shape[0]), np.abs(simpson[:,2:] - simpson3[-1,2:]), 'b-', basey=2)
ax.semilogy(np.arange(simpson2.shape[0]), np.abs(simpson2[:,2:] - simpson3[-1,2:]), 'g-', basey=2)
ax.semilogy(np.arange(simpson3.shape[0]), np.abs(simpson3[:,2:] - simpson3[-1,2:]), 'orange', basey=2)
ax.legend(loc='best')
fig.show()

fig = plt.figure("Trapezoidal Error Difference")
fig.clear()
ax = fig.add_subplot(111)
ax.set_xlabel("levels")
ax.set_ylabel("error")
ax.set_xticks(np.arange(trapezoidal.shape[0]))
ax.grid()
ax.semilogy(np.arange(trapezoidal.shape[0]-1), np.abs(trapezoidal[1:,ix_u1] - trapezoidal[:-1,ix_u1]), 'r-', basey=2)
ax.semilogy(np.arange(simpson.shape[0]-1),     np.abs(simpson[1:,ix_u1] - simpson[:-1,ix_u1]), 'b-', basey=2)
ax.semilogy(np.arange(simpson2.shape[0]-1),    np.abs(simpson2[1:,ix_u1] - simpson2[:-1,ix_u1]), 'g-', basey=2)
ax.semilogy(np.arange(simpson3.shape[0]-1),    np.abs(simpson3[1:,ix_u1] - simpson3[:-1,ix_u1]), 'orange', basey=2)
ax.legend(loc='best')
fig.show()

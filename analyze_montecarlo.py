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

N = 1500
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

ix_show = ix_vel0

fig = plt.figure("U and V")
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

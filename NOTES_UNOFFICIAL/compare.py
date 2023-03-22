import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
import logging
logging.getLogger('matplotlib.font_manager').disabled = True
from flystar import match

t1 = np.genfromtxt('01.XYM/MATCHUP_XYMEEE_final')
t2 = np.genfromtxt('02.ARTSTAR/MATCHUP.XYMEEE')
a = np.genfromtxt('02.ARTSTAR/FAKE.UVW')

print('Original detected {0} stars'.format(len(t1)))
print('Added an artificial star, then detected {0} stars'.format(len(t2)))

plt.figure(1)
plt.clf()
plt.plot(a[0], a[1], 'o', label='FAKE', color='red')
plt.plot(t1[:,0], t1[:,1], '+', label='01.XYM', color='purple')
plt.plot(t2[:,0], t2[:,1], 'x', label='02.ARTSTAR', color='green')
plt.xlabel('u (pix)')
plt.ylabel('v (pix)')
plt.legend()
plt.show()

idx1, idx2, dr, dm = match.match(t1[:,0], t1[:,1], t1[:,2],
                                 t2[:,0], t2[:,1], t2[:,2], 0.1)
print('{0} stars matched'.format(len(dr)))

adx = np.where((t2[:,0] < a[0] + 0.1) &
               (t2[:,0] > a[0] - 0.1) &
               (t2[:,1] < a[1] + 0.1) &
               (t2[:,1] > a[1] - 0.1))[0]
if len(adx) != 1:
    import pdb
    pdb.set_trace()

plt.figure(2)
plt.clf()
Q = plt.quiver(t1[idx1,0], t1[idx1,1],
               t1[idx1,0] - t2[idx2,0],
               t1[idx1,1] - t2[idx2,1],
               angles='xy', scale_units='xy', scale=0.00001)
qk = plt.quiverkey(Q, 0.15, 0.93, 0.001, r'$0.001$ pix',
                   coordinates='figure')

Qa = plt.quiver(t2[adx,0], t2[adx,1],
                t2[adx,0] - a[0],
                t2[adx,1] - a[1],
                angles='xy', scale_units='xy', scale=0.00001,
                color='red')
plt.plot(a[0], a[1], 'o', color='red', label='out $-$ in')
plt.xlabel('u (pix)')
plt.ylabel('v (pix)')
plt.title('01.XYM $-$ 02.ARTSTAR')
plt.legend()
plt.show()

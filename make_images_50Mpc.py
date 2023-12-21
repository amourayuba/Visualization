import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import LambdaCDM
import sys
import matplotlib as mpl
from MAS_library import MASL
import base_for_images as bfi

if __name__ == "__main__":
    params = {'legend.fontsize': 7,
              'legend.handlelength': 2}
    mpl.rcParams['figure.dpi'] = 550
    mpl.rcParams['font.family'] = 'serif'
    plt.rcParams.update(params)

    basefolder = '/home/ayuba/scratch/'

    BoxSize=50
    mlim = 1e13
    zmin, zmax = 100, 150
    sim = sys.argv[1]
    om0 = bfi.sims[sim][0]
    s8 = bfi.sims[sim][1]
    cmp = 'cubehelix'
    smin = int(sys.argv[2])
    cosmol = LambdaCDM(H0=100 * 0.7, Om0=om0, Ode0=1 - om0)
    xmin, xmax = 50, 100
    ymin, ymax = 50, 100
    grid = int(sys.argv[3])
    vmax = float(sys.argv[4])

    for s in range(smin, 118):
        z, pos = bfi.read_pos(basefolder, sim, s)

        inds = np.where((pos[:,2]>=zmin) & (pos[:,2]<=zmax) & (pos[:,1]>=ymin) & (pos[:,1]<=ymax)
                                           & (pos[:,0]>=xmin) & (pos[:,0]<=xmax))[0]
        npos = pos[inds]
        npos[:,2]-= zmin

        age = cosmol.age(z).value

        density_field = np.zeros((grid,grid), dtype=np.float32)

        MASL.MA(npos[:,0:2].astype(np.float32), density_field, BoxSize, MAS='CIC')

        density_field /= np.mean(density_field)

        zr = np.arcsinh(density_field.T -1)

        fig1 = plt.figure()
        plt.imshow(zr, interpolation='bicubic', vmax=6.7, cmap=cmp, origin='lower')

        plt.text(0.08*grid, 0.9*grid, 'z={:1.2f}'.format(z), color='white', size=7)
        plt.text(0.8*grid, 0.9*grid, 't={:1.2f} Gyr'.format(age), color='white', size=6)

        plt.text(0.05*grid, 0.15*grid, r'$\Omega_m$={:1.2f}'.format(om0), color='white', size=5)
        plt.text(0.05*grid, 0.1*grid, r'$\sigma_8$={:1.1f}'.format(s8), color='white', size=5)

        plt.annotate("", xy=(0.4*grid, 0.9*grid), xytext=(0.6*grid, 0.9*grid), arrowprops=dict(arrowstyle='<->', linewidth=1, color='white'))
        plt.text(0.5*grid, 0.915*grid, "10Mpc/h", ha='center', color='white', size=5)

        plt.axis('off')
        fig1.savefig('./zoom_new/zoom{}_{}_bins{}_cmap{}_snap{}_vmax{}.png'.format(BoxSize, sim, grid, cmp, s, vmax), bbox_inches='tight', dpi=1000, transparent=True)
        plt.close(fig1)

import numpy as np
import matplotlib.pyplot as plt
import sys
from astropy.cosmology import LambdaCDM
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

    BoxSize = 500

    zmin, zmax = 100, 105
    print(sys.argv[1])
    sim = sys.argv[1]
    om0 = bfi.sims[sim][0]
    s8 = bfi.sims[sim][1]
    cmp = 'cubehelix'
    smin = int(sys.argv[2])
    print(sys.argv[2])
    grid=int(sys.argv[3])
    print(sys.argv[3])
    vm = sys.argv[4]
    print(sys.argv[4])

    cosmol = LambdaCDM(H0=100 * 0.7, Om0=om0, Ode0=1 - om0)

    for s in range(smin, 119):
        fig = plt.figure()
        z, pos = bfi.read_pos(basefolder, sim, s)
        age = cosmol.age(z).value
        inds = np.where((pos[:,2]>=zmin) & (pos[:,2]<=zmax))[0]

        npos = pos[inds]
        npos[:,2]-= zmin
        ngroups=0

        density_field = np.zeros((grid,grid), dtype=np.float32)
        MASL.MA(npos[:,0:2].astype(np.float32), density_field, BoxSize, MAS='CIC')
        density_field /= np.mean(density_field)
        #zr = np.arcsinh(density_field.T - 1)
        zr = np.log10(density_field.T + 1)
        if vm == 'variable':
            vmax = 4 + s*(2.1/59)
        else:
            vmax = float(vm)

        print(np.max(zr), vmax)
        plt.imshow(zr, interpolation='bicubic', cmap=cmp, vmin=0.05, vmax=vmax, origin='lower')

        plt.text(0.08*grid, 0.9*grid, 'z={:1.2f}'.format(z), color='white', size=7)
        plt.text(0.8*grid, 0.9*grid, 't={:1.2f} Gyr'.format(age), color='white', size=6)

        plt.annotate("", xy=(0.4*grid, 0.9*grid), xytext=(0.6*grid, 0.9*grid), arrowprops=dict(arrowstyle='<->', linewidth=1, color='white'))
        plt.text(0.5*grid, 0.915*grid, "100Mpc/h", ha='center', color='white', size=5)

        plt.axis('off')

        plt.savefig('./full_{:1.0f}/full_{}_bins{}_cmap{}_snap{}_vmax{}.png'.format(grid, sim, grid, cmp, s, vm), bbox_inches='tight', dpi=1000,  transparent=True)
        plt.close(fig)

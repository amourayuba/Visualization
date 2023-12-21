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
    BoxSize=500
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
    mlim = float(sys.argv[5])

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
        if s>17:
            hal_x, hal_y, hal_z = bfi.read_hal_pos(basefolder, sim , s, mlim)
            nhal_x = hal_x[(hal_z>zmin)*(hal_z<zmax)]
            nhal_y = hal_y[(hal_z>zmin)*(hal_z<zmax)]
            ngroups = len(nhal_x)
            plt.scatter((nhal_x/BoxSize)*grid, (nhal_y/BoxSize)*grid, s=5*mlim/1e15, c='r', marker='o')

        plt.text(0.08*grid, 0.9*grid, 'z={:1.2f}'.format(z), color='white', size=7)
        plt.text(0.8*grid, 0.9*grid, 't={:1.2f} Gyr'.format(age), color='white', size=6)
        plt.text(0.08*grid, 0.85*grid, 'N={:1.0f}'.format(ngroups), color='white', size=7)
        plt.text(0.05*grid, 0.15*grid, r'$\Omega_m$={:1.2f}'.format(om0), color='white', size=5)
        plt.text(0.05*grid, 0.1*grid, r'$\sigma_8$={:1.1f}'.format(s8), color='white', size=5)

        plt.annotate("", xy=(0.4*grid, 0.9*grid), xytext=(0.6*grid, 0.9*grid), arrowprops=dict(arrowstyle='<->', linewidth=1, color='white'))
        plt.text(0.5*grid, 0.915*grid, "100Mpc/h", ha='center', color='white', size=5)

        plt.axis('off')
        #ax = plt.gca()
        #ax.set_xticks([])
        #ax.set_yticks([])

        plt.savefig('./full_{:1.0f}/full_{}_bins{}_cmap{}_snap{}_vmax{}_wgroups{:1.2e}.png'.format(grid, sim, grid, cmp, s, vm, mlim), bbox_inches='tight', dpi=1000,  transparent=True)
        plt.close(fig)

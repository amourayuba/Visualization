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

    BoxSize=25
    zmin, zmax = 100, 150
    sim = sys.argv[1]
    print(sim)
    om0 = bfi.sims[sim][0]
    s8 = bfi.sims[sim][1]
    cmp = 'cubehelix'
    smin = int(sys.argv[2])
    print(smin)
    cosmol = LambdaCDM(H0=100 * 0.7, Om0=om0, Ode0=1 - om0)
    xmin, xmax = 85, 95
    ymin, ymax = 58, 68
    xl, yl = 75, 50
    vmax = 1.5
    #grid = 512
    grid = 1500

    for s in range(smin, 119):
        z, pos = bfi.read_pos(basefolder, sim, s)

        inds = np.where((pos[:,2]>=zmin) & (pos[:,2]<=zmax) & (pos[:,1]>=ymin) & (pos[:,1]<=ymax)
                                           & (pos[:,0]>=xmin) & (pos[:,0]<=xmax))[0]
        npos = pos[inds]
        npos[:,2] -= zmin

        age = cosmol.age(z).value

        density_field = np.zeros((grid,grid), dtype=np.float32)

        MASL.MA(npos[:,0:2].astype(np.float32), density_field, BoxSize, MAS='CIC')

        #density_field /= np.mean(density_field)
        sdensity_field = density_field[int(grid*(xmin-xl)/BoxSize):int(grid*(xmax-xl)/BoxSize), int(grid*(ymin-yl)/BoxSize):int(grid*(ymax-yl)/BoxSize)]

        zr = np.log10(sdensity_field.T + 1)

        fig1 = plt.figure()
        plt.imshow(zr, interpolation='bicubic', vmin=0.05, vmax=vmax, cmap=cmp, origin='lower')
        print(np.max(zr))
        ngrid = int(grid*(xmax-xmin)/BoxSize)
        plt.text(0.08*ngrid, 0.9*ngrid, 'z={:1.2f}'.format(z), color='white', size=7)
        plt.text(0.8*ngrid, 0.9*ngrid, 't={:1.2f} Gyr'.format(age), color='white', size=6)

        plt.text(0.05*ngrid, 0.15*ngrid, r'$\Omega_m$={:1.2f}'.format(om0), color='white', size=5)
        plt.text(0.05*ngrid, 0.1*ngrid, r'$\sigma_8$={:1.1f}'.format(s8), color='white', size=5)

        plt.annotate("", xy=(0.4*ngrid, 0.9*ngrid), xytext=(0.6*ngrid, 0.9*ngrid), arrowprops=dict(arrowstyle='<->', linewidth=1, color='white'))
        plt.text(0.5*ngrid, 0.915*ngrid, "2Mpc/h", ha='center', color='white', size=5)

        plt.axis('off')
        fig1.savefig('./zoom_new/zoom{}_{}_bins{}_cmap{}_snap{}_vmax{}.png'.format(BoxSize, sim, grid, cmp, s, vmax), bbox_inches='tight', dpi=700, transparent=True)
        plt.close(fig1)

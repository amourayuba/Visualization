import numpy as np
import pandas as pd
import readgadget
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from scipy.ndimage import gaussian_filter
import scipy.ndimage as ndimage
from MAS_library import MASL

sim_names = ['M25S07', 'M25S08', 'M25S09', 'M03S07', 'M03S08', 'M03S09', 'M35S07', 'M35S08', 'M35S09',
             'Illustris', 'bolshoiP', 'bolshoiW', 'M03S08b', 'm25s85', 'm2s8', 'm4s7', 'm4s8', 'm2s9',
             'm3s8_50', 'm3s8', 'm35s75', 'm4s9', 'm3s9', 'm25s75', 'm2s1', 'm3s7', 'm3s85', 'm2s7', 'm25s8', 'm35s8',
             'm25s9', 'm35s85', 'm3s75', 'm35s9', 'm35s7', 'm25s7']
omegas = [0.25, 0.25, 0.25, 0.3, 0.3, 0.3, 0.35, 0.35, 0.35, 0.309, 0.307, 0.27, 0.3, 0.25, 0.2, 0.4, 0.4, 0.2, 0.3
    , 0.3, 0.35, 0.4, 0.3, 0.25, 0.2, 0.3, 0.3, 0.2, 0.25, 0.35, 0.25, 0.35, 0.3, 0.35, 0.35, 0.25]
sigmas = [0.7, 0.8, 0.9, 0.7, 0.8, 0.9, 0.7, 0.8, 0.9, 0.816, 0.82, 0.82, 0.8, 0.85, 0.8, 0.7, 0.8, 0.9, 0.8
    , 0.8, 0.75, 0.9, 0.9, 0.75, 1.0, 0.7, 0.85, 0.7, 0.8, 0.8, 0.9, 0.85, 0.75, 0.9, 0.7, 0.7]
sims = dict(zip(sim_names, list(zip(omegas, sigmas))))


def read_hal_pos(basefolder, sim, snp, mlim=0):
    folder = basefolder + '{}/AHF/'.format(sim)
    with open(folder + sim + '_prefixes.txt') as file:
        prefs = file.read().splitlines()
    hal = pd.read_table(folder + '/halos/' + prefs[118 - snp] + '.AHF_halos', delim_whitespace=True, header=0)
    shal = hal[hal['Mhalo(4)'] > mlim]
    return np.array(shal['Xc(6)']) / 1000, np.array(shal['Yc(7)']) / 1000, np.array(shal['Zc(8)']) / 1000


def read_pos(basefolder, sim, snp):
    folder = basefolder + '{}/'.format(sim)
    # print(folder)

    snap = 'snapdir_{:03}/snapshot_{:03}'.format(snp, snp)
    head = readgadget.header(folder + snap)
    z = head.redshift

    # print(mean_background)
    ptype = [1]
    return z, readgadget.read_block(folder + snap, "POS ", ptype).astype(np.float64)  # Mpc/h


def smooth_function(f, sigma):
    smoothed = gaussian_filter(f, sigma=sigma)
    return smoothed


def smooth_function2(f, x, y, num_points):
    # interp_func = RegularGridInterpolator((x, y), f)  # or 'linear' or 'quintic'
    interp_func = interp2d(x, y, f, kind='cubic')
    x_smooth = np.linspace(min(x), max(x), num_points)
    y_smooth = np.linspace(min(y), max(y), num_points)
    smoothed = interp_func(x_smooth, y_smooth)
    return smoothed


def adaptive_smoothing(density_field, smoothing_scale_min, smoothing_scale_max):
    # Compute the local density gradient using a Sobel filter
    # gradient = ndimage.sobel(density_field, mode='reflect')

    # Normalize the gradient to [0, 1]
    gradient = 1 - (density_field - np.min(density_field)) / np.max(density_field)

    # Compute the smoothing scale based on the gradient
    smoothing_scale = smoothing_scale_min + (smoothing_scale_max - smoothing_scale_min) * gradient

    # Apply the adaptive smoothing using a Gaussian filter
    smoothed_field = np.zeros_like(density_field)
    for i in range(smoothed_field.shape[0]):
        for j in range(smoothed_field.shape[1]):
            smoothed_field[i, j] = ndimage.gaussian_filter(density_field[i, j], smoothing_scale[i, j])

    return smoothed_field


def anisotropic_diffusion(density_field, iterations, delta_t, kappa):
    smoothed_field = density_field.copy()
    for _ in range(iterations):
        # Compute the gradient magnitude
        gradients = np.gradient(smoothed_field)
        gradient_magnitude = np.sqrt(np.sum(np.square(gradients), axis=0))

        # Compute the diffusion coefficient based on the gradient magnitude
        diffusion_coefficient = np.exp(-kappa * np.square(gradient_magnitude))

        # Update the smoothed field using anisotropic diffusion
        laplacian = ndimage.laplace(smoothed_field)
        smoothed_field += delta_t * diffusion_coefficient * laplacian

    return smoothed_field


def calculate_density_field(data, density_threshold, dense_region_size, low_density_region_size):
    density_field = np.zeros_like(data, dtype=float)

    # Iterate over each point in the data
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            density = calculate_density(data, i, j, density_threshold, dense_region_size, low_density_region_size)
            density_field[i, j] = density

    return density_field


def calculate_density(data, x, y, density_threshold, dense_region_size, low_density_region_size):
    point_value = data[x, y]
    neighborhood = data[max(0, x - low_density_region_size // 2):x + low_density_region_size // 2 + 1,
                   max(0, y - low_density_region_size // 2):y + low_density_region_size // 2 + 1]

    if point_value >= density_threshold:
        neighborhood = data[max(0, x - dense_region_size // 2):x + dense_region_size // 2 + 1,
                       max(0, y - dense_region_size // 2):y + dense_region_size // 2 + 1]

    density = np.mean(neighborhood)
    return density


def plot_xyslice(pos, z, sim, grid=2048, zmin=100, zmax=120, vmin=0.2, vmax=None, BoxSize=500, cmap='hot', snap=118,
                 show=False):
    inds = np.where((pos[:, 2] >= zmin) & (pos[:, 2] <= zmax))[0]
    pos = pos[inds]
    pos[:, 2] -= zmin
    # x,y = pos[:,0], pos[:, 1]
    density_field = np.zeros((grid, grid), dtype=np.float32)
    MASL.MA(pos[:, 0:2].astype(np.float32), density_field, BoxSize, MAS='CIC')
    density_field = density_field / np.mean(density_field)
    dmin, dmax = np.min(np.log10(1 + density_field)), np.max(np.log10(1 + density_field))
    print(dmin, dmax)
    plt.figure()
    if type(vmax) == str:
        plt.imshow(np.log10(1 + density_field), interpolation='gaussian', vmin=vmin, vmax=dmax, cmap=cmap)
    else:
        plt.imshow(np.log10(1 + density_field), interpolation='gaussian', vmin=vmin, vmax=vmax, cmap=cmap)
    plt.text(880, 180, 'z={:1.2f}'.format(z), color='white')

    # plt.title(sim, size=25)
    ax = plt.gca()
    ax.set_xticks([])
    ax.set_yticks([])

    plt.savefig('./snapfigs/{}_bins{}_cmap{}_snap{}_vmax{}.png'.format(sim, grid, cmap, snap, vmax),
                bbox_inches='tight', dpi=1000, facecolor='white', transparent=False)
    # plt.savefig('{}_bins{}_cmap{}_snap{}.pdf'.format(sim, grid, cmap, snap), bbox_inches='tight', dpi=1000, facecolor='white', transparent=False)
    if show:
        plt.show()



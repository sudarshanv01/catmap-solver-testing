"""Plot the coverage and rate maps from the CatMAP pickle file."""
import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import pickle
EPSILON = 1e-150 # Anything lower than this is considered zero
from plot_params import get_plot_params
get_plot_params()

def plot_raw_maps(maps, plot_log=True, plot_params={}):
    """Plot the raw data maps."""

    # Read in data
    maps = np.array(maps, dtype=object)
    descriptors = [] 
    plot_maps = []
    for (desc1, desc2), quantity in maps:
        descriptors.append([desc1, desc2])
        plot_maps.append(quantity)
    descriptors = np.array(descriptors).T
    descriptor_1, descriptor_2 = descriptors
    plot_maps = np.array(plot_maps).T

    # Generate plots
    fig, ax = plt.subplots(1, len(plot_maps), 
                           constrained_layout=True,
                           figsize=(4*len(plot_maps), 3), 
                           sharex=True, sharey=True)

    for species_i, quantity in enumerate(plot_maps):
        quantity = np.array(quantity, dtype=float)
        if plot_log:
            quantity[quantity<EPSILON] = EPSILON
            quantity = np.log10(quantity)
        cax = ax[species_i].scatter(descriptor_1, descriptor_2, 
                   c=quantity, s=10, cmap='jet',**plot_params)
        ax[species_i].set_ylabel(r'$\Delta G_2$ (eV)')
        ax[species_i].set_xlabel(r'$\Delta G_1$ (eV)')
        # Plot the colorbar
        fig.colorbar(cax, ax=ax[species_i])

    return fig, ax


if __name__ == '__main__':
    """Plot the coverage and rate maps where are stored
    in the CatMAP pickle file. The idea is to see the 
    differences in coverage and rates without the 
    interpolation that CatMAP does."""

    # Load the pickle file
    with open('data.pkl', 'rb') as f:
        data = pickle.load(f)

    # Uncomment if you would like to know which other maps exist 
    # for item_type, item in data.items():
    #     print(item_type)

    # Plot the raw points
    plot_params = dict(vmin=0, vmax=1)
    figc, axc = plot_raw_maps(data['coverage_map'], plot_log=False, plot_params=plot_params)
    figr, axr = plot_raw_maps(data['rate_map'], plot_log=True)
    figp, axp = plot_raw_maps(data['production_rate_map'], plot_log=True)

    figc.savefig('output/coverage_map.png', dpi=300)
    figr.savefig('output/rate_map.png', dpi=300)
    figp.savefig('output/production_rate_map.png', dpi=300)


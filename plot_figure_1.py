"""Plot Figure 1 of the manuscript."""
import os
from glob import glob
import readline
import numpy as np
import json
import csv
import matplotlib.pyplot as plt
from plot_params import get_plot_params
import pickle
import string
get_plot_params()

if __name__ == '__main__':
    """Plot the initial guess Boltzmann coverage
    of a particular solver and convergence level
    and compare the values against the conventional
    coverage based solver."""

    # Provide the filenames of the different coverages
    files = {'boltzmann': 'initial_guess.csv', 
             'solution': 'solution.csv',
             }

    # read in calculation details    
    with open('details.json', 'r') as f:
        details = json.load(f)

    reaction = 'co_oxidation'
    solver = 'coverages' # 'numbers_fix_xstar'
    conv_level = 'midtight'

    # Get the data from the solution file
    foldername = f'{conv_level}/{reaction}/{solver}'
    print(foldername)

    mkm_file = glob(os.path.join(foldername, '*.mkm'))[0]
    # Get all the quantities from the input file
    exec(open(os.path.join(foldername, f"input.mkm")).read())
    # Read the adsorbate names
    # adsorbate_names_files = os.path.join(foldername, f"adsorbate_names.log")
    # adsorbate_names = np.genfromtxt(adsorbate_names_files, dtype=str) 
    # print(adsorbate_names)

    # titles = [ r'$\theta$ %s*'%adsorbate_names[i] for i in range(len(adsorbate_names)) ] 
    titles = [r'$\theta_{\rm CO}$', r'$\theta_{\rm O}$', r'$\theta_{\rm *}$']

    # Plot separate plots for the initial coverage guess
    # from Boltzmann and the solution coverages
    for k, (handle, filename) in enumerate(files.items()):

        # If not yet converged, skip
        if not os.path.isfile(os.path.join(foldername, filename)):
            raise FileNotFoundError(f'{filename} not found in {foldername}')

        with open(os.path.join(foldername, filename), 'r') as f:
            reader = csv.reader(f)
            data = list(reader)
        
        # Get all the data from the solution file
        data = np.array(data, dtype=float).T

        # Store the descriptors
        desc1 = data[0]
        desc2 = data[1]
        thetas = []

        for i in range(len(desc1)):
            dG1 = desc1[i]
            dG2 = desc2[i]
            x_list = []
            for j in range(len(data[2:])):
                x_list.append(data[2:][j][i])
            thetas.append(x_list)
        thetas = np.array(thetas) 
        sum_thetas = np.sum(thetas, axis=1)
        thetas = thetas.T

        fig, ax = plt.subplots(1, len(thetas),
                               constrained_layout=True,
                               figsize=(6.75, 2.5),
                               sharex=True, sharey=True)

        for j in range(len(thetas)):
            cax = ax[j].scatter(desc1, desc2, c=thetas[j], s=10, cmap='jet', vmin=0, vmax=1)
            ax[j].set_xlabel(r'$\Delta G$ %s* (eV)'%descriptor_names[0].replace('_s', ''))
            if j==0:
                ax[j].set_ylabel(r'$\Delta G$ %s* (eV)'%descriptor_names[1].replace('_s', ''))
            if j == len(thetas)-1:
                # Plot the colorbar
                fig.colorbar(cax, ax=ax[j])

        # Label in alphabetical order the axes on the top left above the plot
        for j, a in enumerate(ax):
            a.text(0.05, 1.15, string.ascii_lowercase[j]+')',
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=a.transAxes,
                    fontsize=10) 

        for i, a in enumerate(ax):
            a.set_title(titles[i], fontsize=10)
        fig.savefig(f'figures/{handle}_{reaction}_{solver}.png', dpi=300) 
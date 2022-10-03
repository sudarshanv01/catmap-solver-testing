"""Plot the solution of CO coverages with different solver parameters."""
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
import json
import os
import string
from glob import glob
from plot_params import get_plot_params
get_plot_params()

def plot_coverage_maps(foldername, ax, ti=0):
    """Plot the initial guesses for the coverage
    and the final converage after the solution for
    the different cases supplied through folder_names."""

    # Find the name of the mkm file in this folder
    print(f'Plotting {foldername}')
    mkm_file = glob(os.path.join(foldername, '*.mkm'))[0]
    case_name = mkm_file.split('/')[-1].split('.')[0]

    # Get information about the descriptors
    try:
        exec(open(os.path.join(foldername, f"input.log")).read())
        titles = [ r'$\theta$ %s*'%adsorbate_names[i].replace('_s', '') for i in range(len(adsorbate_names)) ] 
        site_names_withoutg = [site_name for site_name in site_names if 'g' not in site_name] 
        titles += [r'$\theta %s*$'%site_name for site_name in site_names_withoutg]
        use_titles = True
    except FileNotFoundError:
        print('No log file found')
        adsorbate_names = ['1', '2']
        descriptor_names = ['1', '2']
        use_titles = False

    for k, (handle, filename) in enumerate(files.items()):

        # If not yet converged, skip
        if not os.path.isfile(os.path.join(foldername, filename)):
            continue

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

        cax = ax.scatter(desc1, desc2, c=thetas[ti], s=10, cmap='jet', vmin=0, vmax=1)

        return cax


if __name__ == '__main__':
    """Plot the initial coverage and the solution."""

    # Provide the filenames of the different coverages
    files = {'solution': 'solution.csv'}

    # read in calculation details    
    with open('details.json', 'r') as f:
        details = json.load(f)

    reaction = 'co_oxidation'
    solvers = details['solvers']
    convergence_levels = details['convergence_levels']
    # Remove tightplus from convergence_levels
    convergence_levels.remove('tightplus')

    # Loop over the different solvers
    for i, solver in enumerate(solvers): 
        fig, ax = plt.subplots(1, len(convergence_levels), 
                              figsize=(6.75, 1.5), constrained_layout=True)
        for j, conv_level in enumerate(convergence_levels):
            cax = plot_coverage_maps(f'{conv_level}/{reaction}/{solver}', ax[j])
            ax[j].set_xlabel(r'$\Delta G$ O* (eV)', fontsize=8)
            if j == 0:
                ax[j].set_ylabel(r'$\Delta G$ CO* (eV)', fontsize=8)
            if j == len(convergence_levels) -1 :
                fig.colorbar(cax, ax=ax[j])
            ax[j].set_title(conv_level, fontsize=8)

        # Label in alphabetical order the axes on the top left above the plot
        for j, a in enumerate(ax):
            a.text(0.05, 1.15, string.ascii_lowercase[j]+')',
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=a.transAxes,
                    fontsize=6) 

        fig.savefig(f'figures/solver_parameters_{solver}.png', dpi=300)

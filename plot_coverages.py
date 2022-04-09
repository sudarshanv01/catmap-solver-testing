"""Plot the initial guess for debugging the new solver."""
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
import json
import os
from glob import glob
from plot_params import get_plot_params
get_plot_params()

def plot_coverage_maps(foldername):
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

        fig, ax = plt.subplots(1, len(thetas),
                               constrained_layout=True,
                               figsize=(3.5*len(thetas), 3.),
                               sharex=True, sharey=True)

        for j in range(len(thetas)):
            cax = ax[j].scatter(desc1, desc2, c=thetas[j], s=10, cmap='jet', vmin=0, vmax=1)
            ax[j].set_xlabel(r'$\Delta G$ %s* (eV)'%descriptor_names[0].replace('_s', ''))
            if j==0:
                ax[j].set_ylabel(r'$\Delta G$ %s* (eV)'%descriptor_names[1].replace('_s', ''))
            if j == len(thetas)-1:
                # Plot the colorbar
                fig.colorbar(cax, ax=ax[j])
            if use_titles:
                ax[j].set_title(titles[j])

        fig.savefig('output/'+handle+'_'+foldername.replace('/', '_')+'.png', dpi=300) 

if __name__ == '__main__':
    """Plot the initial coverage and the solution."""

    # Provide the filenames of the different coverages
    files = {'boltzmann': 'initial_guess.csv', 
             'solution': 'solution.csv',
             }

    # read in calculation details    
    with open('details.json', 'r') as f:
        details = json.load(f)

    reactions = details['reactions'] 
    solvers = details['solvers']
    convergence_levels = details['convergence_levels']

    # Loop over the different solvers
    for conv_level in convergence_levels:
        for r, reaction in enumerate(reactions):
            for i, solver in enumerate(solvers): 
                plot_coverage_maps(f'{conv_level}/{reaction}/{solver}')

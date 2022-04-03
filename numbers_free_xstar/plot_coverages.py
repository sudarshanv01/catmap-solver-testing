"""Plot the initial guess for debugging the new solver."""
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
from glob import glob
from plot_params import get_plot_params
get_plot_params()

if __name__ == '__main__':
    """Read in the csv file initial_guess.csv the first two 
    columns of the csv file are the descriptors, the others
    are all the initial x0 values. Convert the x0 values 
    to theta and see if the coverages are what is expected."""

    # Generate all the variables 
    # Find the name of the mkm file in this folder
    mkm_file = glob('*.mkm')[0]
    case_name = mkm_file.split('.')[0]
    exec(open(f"{case_name}.log").read())

    files = {'boltzmann': 'initial_guess.csv', 
             'solution': 'solution.csv',
             }
    titles = [ r'$\theta$ %s*'%adsorbate_names[i].replace('_s', '') for i in range(len(adsorbate_names)) ] 
    site_names_withoutg = [site_name for site_name in site_names if 'g' not in site_name] 
    titles += [r'$\theta %s*$'%site_name for site_name in site_names_withoutg]

    # Generate a contour plot of the coverages
    # each coverage has a separate plot

    for k, (handle, filename) in enumerate(files.items()):

        with open(filename, 'r') as f:
            reader = csv.reader(f)
            data = list(reader)

        data = np.array(data, dtype=float).T

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

        fig, ax = plt.subplots(1, len(thetas), constrained_layout=True, figsize=(3.5*len(thetas), 5), sharex=True, sharey=True)

        for j in range(len(thetas)):
            cax = ax[j].scatter(desc1, desc2, c=thetas[j], s=10, cmap='jet', vmin=0, vmax=1)
            ax[j].set_xlabel(r'$\Delta G$ %s* (eV)'%descriptor_names[0].replace('_s', ''))
            if j==0:
                ax[j].set_ylabel(r'$\Delta G$ %s* (eV)'%descriptor_names[1].replace('_s', ''))
            if j == len(thetas)-1:
                # Plot the colorbar
                fig.colorbar(cax, ax=ax[j])
            ax[j].set_title(titles[j])
        
        fig.savefig(f'output/{handle}_coverage.png', dpi=300)
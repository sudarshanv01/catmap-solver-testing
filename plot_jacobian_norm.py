"""Compare the analytical and numerical norm of the Jacobian.""" 
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
from pprint import pprint
from plot_params import get_plot_params
import os
get_plot_params()
import matplotlib as mpl
mpl.rcParams['lines.markersize'] = 4
import string

def prepare_data(data):
    """Prepare the data for plotting."""
    results = []
    for i, dat in enumerate(data):
        results.append([ [data[i][0], data[i][1]], data[i][2], data[i][3], data[i][4] ])
    return results

def plot_data(results, ax, color='tab:red'):
    """Plot the numerical and analytical norm of the Jacobian."""

    prev_desc = results[0][0]
    iter_data = []

    for desc, iteration, norm_analytical, norm_numerical in results:
        if prev_desc == desc:
            # Store the descriptor and the error
            iter_data.append([iteration, norm_analytical, norm_numerical])
        else:
            # Create a new cycle and plot what was already there
            iter_data = np.array(iter_data)
            # Sort iter_data by the first column
            # iter_data = iter_data[iter_data[:,0].argsort()]
            # Split the array into chunks such that the first column
            # is always monotonically increasing
            chunks = np.split(iter_data, np.where(np.diff(iter_data[:,0]) < 0)[0]+1)

            for iter_data in chunks:
                ax.plot(iter_data[:,1], iter_data[:,2], 'o', alpha=0.5, color=color)

            iter_data = []
            iter_data.append([iteration, norm_analytical, norm_numerical])

        prev_desc = desc

if __name__ == '__main__':
    """Compare the norm of the Jacobian computed analytically and
    numerically for different types of solver."""

    comparison_solvers_all = [
        [ 'numbers_fix_xstar', 'numbers_free_xstar'],
        [ 'numbers_ads_ads_fix_xstar', 'numbers_ads_ads_free_xstar'],
    ]
    labels = ['fix-x', 'free-x']
    SOLVER_COLOR = ['tab:blue', 'tab:green']
    titles = ['CO oxidation', 'CO oxidation with ads-ads']

    fig, ax = plt.subplots(1, len(comparison_solvers_all), figsize=(5.5, 2.5), constrained_layout=True)

    for i, (solver_1, solver_2 ) in enumerate(comparison_solvers_all):
        
        # Import the Jacobian error for the first solver
        with open(os.path.join(solver_1, 'jacobian_norm.csv'), 'r') as f:
            reader = csv.reader(f)
            data_1 = list(reader)
        data_1 = np.array(data_1, dtype=float)
        # Import the Jacobian error for the second solver
        with open(os.path.join(solver_2, 'jacobian_norm.csv'), 'r') as f:
            reader = csv.reader(f)
            data_2 = list(reader)
        data_2 = np.array(data_2, dtype=float)

        # Prepare the data for plotting
        data_1 = prepare_data(data_1)
        data_2 = prepare_data(data_2)

        # Plot the data
        plot_data(data_1, ax[i], color=SOLVER_COLOR[0])
        plot_data(data_2, ax[i], color=SOLVER_COLOR[1])

    # Plot axis on log10 scale
    for a in ax.flatten():
        a.set_yscale('log') 
        a.set_xscale('log') 
        a.set_xlabel(r'$\left| \mathrm{J}(x)\right |$ (Analytical)')
        a.set_ylabel(r'$\left| \mathrm{J}(x)\right |$ (Numerical)')
        a.set_aspect('equal')
        # Make a parity line
        xlims = a.get_xlim()
        ylims = a.get_ylim()
        a.plot([xlims[0], xlims[1]], [xlims[0], xlims[1]], 'k--')
    for label, color in zip(labels, SOLVER_COLOR):
        ax[0].plot([], [], 'o', color=color, label=label)
    for title, a in zip(titles, ax):
        a.set_title(title)
    
    # Annotate an alphabet to indicate the subplot number
    # at the top left of the figure
    for i, a in enumerate(ax):
        a.text(0.05, 0.95, string.ascii_lowercase[i] + ')', transform=a.transAxes, 
                va='top')

    ax[0].legend(loc='lower right')

    fig.savefig(f'output/norm_jacobian.png', dpi=300)
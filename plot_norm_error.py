"""Make a parity plot of the errors form two solvers."""
import numpy as np
import matplotlib.pyplot as plt
import os
import csv
import mpmath as mp
from pprint import pprint
from plot_params import get_plot_params
from collections import defaultdict
get_plot_params()
ACCURACY = 1e-75
import matplotlib as mpl
mpl.rcParams['lines.markersize'] = 4
import string

def prepare_data(data):
    """Prepare the data for plotting."""
    results = []
    for i, dat in enumerate(data):
        results.append([ [data[i][0], data[i][1]], data[i][2:] ])
    return results

def get_data_to_plot(results):
    """Get the exact points to plot from the results."""
    # Initialise the variables
    prev_desc = results[0][0]
    iter_data = []
    # Store data per descriptor
    data_to_plot = defaultdict(lambda: defaultdict(list))
    # Iterate over the results
    for desc, (iteration, error) in results:
        if prev_desc == desc:
            # Store the descriptor and the error
            iter_data.append([iteration, error])
        else:
            # Create a new cycle and plot what was already there
            iter_data = np.array(iter_data)
            # Sort iter_data by the first column
            # iter_data = iter_data[iter_data[:,0].argsort()]
            # Split the array into chunks such that the first column
            # is always monotonically increasing
            chunks = np.split(iter_data, np.where(np.diff(iter_data[:,0]) < 0)[0]+1)

            for iter_data in chunks:
                if np.min(iter_data[:,1]) <= ACCURACY:
                    memo_desc = str(desc[0]) + '&' + str(desc[1])
                    data_to_plot[memo_desc]['success'].extend(iter_data.T)
                else:
                    memo_desc = str(desc[0]) + '&' + str(desc[1])
                    data_to_plot[memo_desc]['failure'].extend(iter_data.T)

            iter_data = []
            iter_data.append([iteration, error])

        prev_desc = desc

    return data_to_plot

def plot_data(data_to_plot_1, data_to_plot_2, data_to_plot_3, ax):
    """Plot the data to compare two solvers."""
    success_runs = [0, 0, 0]
    failure_runs = [0, 0, 0]
    for desc, error_1 in data_to_plot_1.items():
        if len(error_1['success']) > 0:
            ax[0,0].plot(error_1['success'][0], error_1['success'][1], 'o',
                       markerfacecolor='none',  alpha=0.2 )
            success_runs[0] += 1
        if len(error_1['failure']) > 0:
            failure_runs[0] += 1
            ax[1,0].plot(error_1['failure'][0], error_1['failure'][1], 'o',
                markerfacecolor='none',  alpha=0.2)

    for desc, error_2 in data_to_plot_2.items():
        if len(error_2['success']) > 0:
            ax[0,1].plot(error_2['success'][0], error_2['success'][1], 'o',
                markerfacecolor='none',  alpha=0.2 )
            success_runs[1] += 1
        if len(error_2['failure']) > 0:
            failure_runs[1] += 1
            ax[1,1].plot(error_2['failure'][0], error_2['failure'][1], 'o',
                markerfacecolor='none', alpha=0.2)
    
    for desc, error_3 in data_to_plot_3.items():
        if len(error_3['success']) > 0:
            ax[0,2].plot(error_3['success'][0], error_3['success'][1], 'o',
                markerfacecolor='none',  alpha=0.2)
            success_runs[2] += 1
        if len(error_3['failure']) > 0:
            failure_runs[2] += 1
            ax[1,2].plot(error_3['failure'][0], error_3['failure'][1], 'o',
                markerfacecolor='none', alpha=0.2)

    print('Success runs:', success_runs)
    print('Failure runs:', failure_runs)
    # Annotate the failed run number on the plot with a bbox
    for i, ax_ in enumerate(ax[-1,:].flatten()):
        if failure_runs[i] > 0:
            ax_.annotate('Failed points: ' + str(failure_runs[i]), xy=(0.4, 0.8),
                bbox=dict(boxstyle="round", fc="0.8"),
                xycoords='axes fraction', color='k')
        ax_.set_title('Failed: ' + labels[i], fontsize=9)
    for i, ax_ in enumerate(ax[0,:].flatten()):
        ax_.set_title('Successful: ' + labels[i], fontsize=9)
    
if __name__ == '__main__':
    """Compare two solvers on the basis of their
    errors. Each point is an iteration, if the 
    error from the first solver is larger than the
    error from the second solver, the point will be 
    on the left side of the parity plot."""

    comparison_solvers_all = [
        ['coverages', 'numbers_fix_xstar', 'numbers_free_xstar'],
        ['coverages_ads_ads', 'numbers_ads_ads_fix_xstar', 'numbers_ads_ads_free_xstar'],
        ['coverages_steam_reforming', 'numbers_fix_xstar_steam_reforming', 'numbers_free_xstar_steam_reforming'],
        ['coverages_etoh', 'numbers_fix_xstar_etoh', 'numbers_free_xstar_etoh'],
    ]
    labels = ['Coverages', 'fix-x', 'free-x']

    for i, (solver_1, solver_2, solver_3) in enumerate(comparison_solvers_all):
        # Plot separate figures for the comparison of the solvers
        fig, ax = plt.subplots(2, 3, figsize=(6.75,4.25), constrained_layout=True)

        # Import data from csv files
        with open(os.path.join(solver_1, 'error_log.csv'), 'r') as f:
            reader = csv.reader(f)
            data_1 = list(reader)
        with open(os.path.join(solver_2, 'error_log.csv'), 'r') as f:
            reader = csv.reader(f)
            data_2 = list(reader)
        with open(os.path.join(solver_3, 'error_log.csv'), 'r') as f:
            reader = csv.reader(f)
            data_3 = list(reader)

        # Prepare the data for plotting
        data_1 = np.array(data_1, dtype=float)
        data_2 = np.array(data_2, dtype=float)
        data_3 = np.array(data_3, dtype=float)

        # Prepare the data for plotting
        data_1 = prepare_data(data_1)
        data_2 = prepare_data(data_2)
        data_3 = prepare_data(data_3)

        # Get the data to plot
        data_to_plot_1 = get_data_to_plot(data_1)
        data_to_plot_2 = get_data_to_plot(data_2)
        data_to_plot_3 = get_data_to_plot(data_3)

        # Plot the result
        plot_data(data_to_plot_1, data_to_plot_2, data_to_plot_3, ax)

        # Label the axes
        for a in ax.flatten():
            a.set_xlabel('Iteration')
            a.set_ylabel('Error')

        # Plot axis on log10 scale
        for a in ax.flatten():
            a.set_yscale('log') 

        for j, a in enumerate(ax.flatten()):
            a.text(-0.1, 1.1, string.ascii_lowercase[j] + ')', transform=a.transAxes, 
                    va='top')

        fig.savefig(f'output/norm_error_case_{i}.png', dpi=300)
         
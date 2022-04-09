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
import matplotlib as mpl
import json
mpl.rcParams['lines.markersize'] = 4
import string

def prepare_data(data):
    """Prepare the data for plotting."""
    results = []
    for i, dat in enumerate(data):
        results.append([ [data[i][0], data[i][1]], data[i][2:] ])
    return results

def get_data_to_plot(results, accuracy):
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
                if np.min(iter_data[:,1]) <= accuracy:
                    memo_desc = str(desc[0]) + '&' + str(desc[1])
                    data_to_plot[memo_desc]['success'].extend(iter_data.T)
                else:
                    memo_desc = str(desc[0]) + '&' + str(desc[1])
                    data_to_plot[memo_desc]['failure'].extend(iter_data.T)

            iter_data = []
            iter_data.append([iteration, error])

        prev_desc = desc

    return data_to_plot

def plot_data(data_to_plot_1, ax):
    """Plot the data to compare two solvers."""
    success_runs = 0 
    failure_runs = 0 
    for desc, error_1 in data_to_plot_1.items():
        if len(error_1['success']) > 0:
            ax[0].plot(error_1['success'][0], error_1['success'][1], 'o',
                       markerfacecolor='none',  alpha=0.2 )
            success_runs += 1
        if len(error_1['failure']) > 0:
            failure_runs += 1
            ax[1].plot(error_1['failure'][0], error_1['failure'][1], 'o',
                markerfacecolor='none',  alpha=0.2)

    print('Success runs:', success_runs)
    print('Failure runs:', failure_runs)
    return success_runs, failure_runs

if __name__ == '__main__':
    """Compare two solvers on the basis of their
    errors. Each point is an iteration, if the 
    error from the first solver is larger than the
    error from the second solver, the point will be 
    on the left side of the parity plot."""

    # read in calculation details    
    with open('details.json', 'r') as f:
        details = json.load(f)

    reactions = details['reactions'] 
    solvers = details['solvers']
    convergence_levels = details['convergence_levels']

    labels = ['Coverages', 'fix-x', 'free-x']

    convergence_data = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(int))))

    for conv_level in convergence_levels:
        for r, reaction in enumerate(reactions):

            # Plot separate figures for the comparison of the solvers
            fig, ax = plt.subplots(2, 3, figsize=(6.75,4.25), constrained_layout=True)

            for i, solver in enumerate(solvers): 

                # Import data from csv files
                if not os.path.isfile(os.path.join(conv_level, reaction, solver, 'error_log.csv')):
                    continue

                with open(os.path.join(conv_level, reaction, solver, 'error_log.csv'), 'r') as f:
                    reader = csv.reader(f)
                    data = list(reader)
                
                print('Plotting:', reaction, solver, conv_level)
                
                # Get the accuracy of the solver
                with open(os.path.join(conv_level, reaction, solver, 'solver_specifics.json')) as f:
                    solver_specifics = json.load(f)
                ACCURACY = solver_specifics['tolerance']

                # Prepare the data for plotting
                data = np.array(data, dtype=float)

                # Prepare the data for plotting
                data = prepare_data(data)

                # Get the data to plot
                data_to_plot = get_data_to_plot(data, accuracy=ACCURACY)

                # Plot the result
                success, fail = plot_data(data_to_plot, ax[:,i])

                convergence_data[conv_level][reaction][solver]['success'] = success
                convergence_data[conv_level][reaction][solver]['failure'] = fail
                

            # Label the axes
            for a in ax.flatten():
                a.set_xlabel('Iteration')
                a.set_ylabel('Error')
                a.set_yscale('log') 

            for j, a in enumerate(ax.flatten()):
                a.text(-0.1, 1.1, string.ascii_lowercase[j] + ')', transform=a.transAxes, 
                        va='top')

            fig.savefig(f'output/norm_error_{conv_level}_{reaction}.png', dpi=300)
            plt.close(fig)

    # Save the convergence data
    with open('output/convergence_data.json', 'w') as f:
        json.dump(convergence_data, f, indent=4) 
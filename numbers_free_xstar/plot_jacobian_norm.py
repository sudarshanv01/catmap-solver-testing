"""Plot the norm of the Jacobian.""" 
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
from pprint import pprint
from plot_params import get_plot_params
get_plot_params()

if __name__ == '__main__':
    """Read the data from the jacobian_norm.csv file and plot
    the error as a function of the number of steps."""
    with open('jacobian_norm.csv', 'r') as f:
        reader = csv.reader(f)
        data = list(reader)

    data = np.array(data, dtype=float)
    results = []
    for i, dat in enumerate(data):
        results.append([ [data[i][0], data[i][1]], data[i][2:] ])

    # Plot two separate figures with the first graph showing
    # the datapoints where the error is monotically decreasing
    # and the second shows all the points where this condition does 
    # not hold, or the solver breaks within the first iteration
    fig, ax = plt.subplots(1, 1, figsize=(6,5.5), constrained_layout=True)

    prev_desc = results[0][0]
    iter_data = []
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
                ax.plot(iter_data[:,0], iter_data[:,1], '.-', alpha=0.5)

            iter_data = []
            iter_data.append([iteration, error])

        prev_desc = desc

    # Plot axis on log10 scale
    ax.set_yscale('log') 
    ax.set_xlabel('Number of iterations')
    ax.set_ylabel(r'$\left| \mathrm{J}(x)\right |$')
    fig.savefig('output/norm_jacobian.png')
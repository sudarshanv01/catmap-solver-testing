"""Plot the error through the norm of the f=0 function."""
import numpy as np
import matplotlib.pyplot as plt
import csv
import mpmath as mp
from pprint import pprint
from plot_params import get_plot_params
get_plot_params()

if __name__ == '__main__':
    """Read the data from the error_log.csv file and plot
    the error as a function of the number of steps."""
    with open('error_log.csv', 'r') as f:
        reader = csv.reader(f)
        data = list(reader)

    ACCURACY = 1e-75

    data = np.array(data, dtype=float)
    results = []
    for i, dat in enumerate(data):
        results.append([ [data[i][0], data[i][1]], data[i][2:] ])

    # Plot two separate figures with the first graph showing
    # the datapoints where the error is monotically decreasing
    # and the second shows all the points where this condition does 
    # not hold, or the solver breaks within the first iteration
    fig, ax = plt.subplots(1, 2, figsize=(12,5.5), constrained_layout=True)

    prev_desc = results[0][0]
    iter_data = []
    number_of_error_points = 0
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
                    ax[0].plot(iter_data[:,0], iter_data[:,1], '.-', alpha=0.5)
                else:
                    number_of_error_points += 1
                    ax[1].plot(iter_data[:,0], iter_data[:,1], '.-', alpha=0.5)

            iter_data = []
            iter_data.append([iteration, error])

        prev_desc = desc
        # Count the number of points
    
    print('Number of points that dont converge:', number_of_error_points)
    # Plot axis on log10 scale
    for a in ax:
        a.set_yscale('log') 
        a.set_xlabel('Number of iterations')
        a.set_ylabel('Error')
        a.axhline(y=ACCURACY, color='k', linestyle='--')
    fig.savefig('output/norm_error.png')
         
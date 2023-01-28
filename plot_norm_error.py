import string
import os
import csv
import json
from collections import defaultdict
import itertools

from typing import List, Tuple, Dict, Any, Union, Optional

import numpy as np

import pandas as pd

import mpmath as mp

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from matplotlib.lines import Line2D

from plot_params import get_plot_params

get_plot_params()

import logging

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)


def evaulate_data(
    data: pd.DataFrame,
    accuracy: float,
    ax: plt.Axes,
    success_settings: dict = {},
    failure_settings: dict = {},
) -> Tuple[float]:
    """Plot the data to compare two solvers."""

    # Remove rows where the error is equal to 0
    data = data[data["error"] != 0]

    # Store success and failure data
    success = 0
    failure = 0

    for (desc1, desc2), dat in data.groupby(["descriptor1", "descriptor2"]):

        # Split the pandas dataframe into a series of dataframes where the
        # "iteration" column is not strictly monotonically increasing
        chunks = np.split(dat, np.where(np.diff(dat["iteration"]) < 0)[0] + 1)

        # Iterate over the chunks
        for chunk in chunks:
            # Make sure that the chunks are ordered with respect to iteration
            chunk = chunk.sort_values(by=["iteration"])
            # Check if the error of the last row is lower than the accuracy
            if chunk["error"].iloc[-1] < accuracy: 
                success += 1
                ax.plot(chunk["iteration"], chunk["error"], **success_settings)
            else:
                # Make sure that none of the errors are lower than the accuracy
                assert not any(chunk["error"] < accuracy)
                failure += 1
                ax.plot(chunk["iteration"], chunk["error"], **failure_settings)

    return success, failure


if __name__ == "__main__":
    """Plot the norm error for the solvers."""

    # Read in the details of the comparison
    with open(os.path.join("config", "details.json"), "r") as f:
        details = json.load(f)
    logger.info("Details: %s", details)

    # Store the relevant data to find the directories to look in
    reactions = details["reactions"]
    solvers = details["solvers"]
    convergence_levels = details["convergence_levels"]

    # Labels for the solvers that are being compared
    labels = [r"Coverages solver, $\theta$", r"Numbers solver, $x$"]

    # Create a dataframe to store the data
    convergence_data = pd.DataFrame(
        columns=["solver", "reaction", "convergence_level", "success", "failure"]
    )

    # Plot settings for success and failure, remove fill for points and no lines
    success_settings = {
        "color": "tab:green",
        "alpha": 0.5,
        "marker": "o",
        "fillstyle": "none",
        "linestyle": "--",
    }
    failure_settings = {
        "color": "tab:red",
        "alpha": 0.5,
        "marker": "o",
        "fillstyle": "none",
        "linestyle": "--",
    }

    # Reaction name
    reaction_name = {
        'methanol_synthesis': 'Methanol synthesis',
        'co_oxidation': 'CO oxidation',
        'co_hydrogenation': 'CO hydrogenation',
        'co_oxidation_ads_ads': 'CO oxidation (ads-ads)',
    }

    # Using itertools to iterate over covergence_levels and reactions
    for conv_level, reaction in itertools.product(convergence_levels, reactions):

        # Plot separate figures for the comparison of the solvers
        fig, ax = plt.subplots(1, 2, figsize=(4.5, 2.5), constrained_layout=True)

        # Iterate over the solvers
        for i, solver in enumerate(solvers):

            # Import data from csv files
            if not os.path.isfile(
                os.path.join(conv_level, reaction, solver, "error_log.csv")
            ):
                logger.warning(
                    f"File {os.path.join(conv_level, reaction, solver, 'error_log.csv')} does not exist."
                )
                continue

            # Use pandas to read the csv file
            data = pd.read_csv(
                os.path.join(conv_level, reaction, solver, "error_log.csv")
            )
            # Set labels for data
            data.columns = ["descriptor1", "descriptor2", "iteration", "error"]
            logger.info("Plotting: %s, %s, %s", reaction, solver, conv_level)

            # Get the accuracy of the solver
            with open(
                os.path.join(conv_level, reaction, solver, "solver_specifics.json")
            ) as f:
                solver_specifics = json.load(f)
            ACCURACY = solver_specifics["tolerance"]
            logger.info("Accuracy: %s", ACCURACY)

            # Get the data to plot
            success, failure = evaulate_data(
                data,
                accuracy=ACCURACY,
                ax=ax[i],
                success_settings=success_settings,
                failure_settings=failure_settings,
            )

            # Store the data in the dataframe
            data_to_df = {
                "solver": solver,
                "reaction": reaction,
                "convergence_level": conv_level,
                "success": success,
                "failure": failure,
            }

            # Concatenate the data to the dataframe
            convergence_data = pd.concat(
                [convergence_data, pd.DataFrame(data_to_df, index=[0])]
            )

            for j, a in enumerate(ax.flatten()):
                a.text(
                    -0.1,
                    1.1,
                    string.ascii_lowercase[j] + ")",
                    transform=a.transAxes,
                    va="top",
                    fontsize=10,
                )

            # Label the figure
            ax[i].set_title(labels[i])

            # Make a dashed line to show the accuracy
            ax[i].axhline(y=ACCURACY, color="black", linestyle="--")

        # Set the x and y labels
        ax[0].set_xlabel("Iteration")
        ax[1].set_xlabel("Iteration")
        ax[0].set_ylabel("Maximum error")
        # Plot y on a log scale
        ax[0].set_yscale("log")
        ax[1].set_yscale("log")
        # Set the log axis ticks on the minor ticks
        ax[0].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))
        ax[1].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))

        # Make super title for the figure with the reaction name
        fig.suptitle(reaction_name[reaction], fontsize=11)

        # Set legend for the figure on axis 1, where red is the failure and green is the success
        # Make it on the right hand side of the figure
        ax[1].plot([], [], color='tab:green', label="Success")
        ax[1].plot([], [], color='tab:red', label="Failure")
        ax[1].legend(
            loc="center left",
            frameon=False,
            bbox_to_anchor=(1.04, 1),
            borderaxespad=0
        )


        fig.savefig(f"output/norm_error_{conv_level}_{reaction}.png", dpi=300)
        plt.close(fig)

    # Sum up the success and failure for each solver and store it as the sum 
    # of points column
    convergence_data["sum_of_points"] = convergence_data["success"] + convergence_data["failure"]
    # Take a difference between the success and failure for each solver and store it as the difference
    # between points column
    convergence_data["difference_between_points"] = convergence_data["success"] - convergence_data["failure"]
    # Save the dataframe
    convergence_data.to_csv("output/convergence_data.csv", index=False)
    logger.info(convergence_data)

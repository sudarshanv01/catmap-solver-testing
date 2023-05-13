import os

import shutil

import json

import glob

from string import Template

import argparse

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

import yaml

from typing import Dict, List, Tuple, Any

import numpy as np

import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

from template_reactions import templates

from plot_params import get_plot_params

get_plot_params()


def get_cli_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--resolution", type=int, default=10, help="Number of points to plot."
    )
    return parser.parse_args()


def evaulate_data(
    data: pd.DataFrame,
    accuracy: float,
    ax: plt.Axes,
    success_settings: dict = {
        "color": "tab:green",
        "marker": ".",
        "linestyle": "-",
        "fillstyle": "none",
        "alpha": 0.5,
    },
    failure_settings: dict = {
        "color": "tab:red",
        "marker": ".",
        "linestyle": "-",
        "fillstyle": "none",
        "alpha": 0.5,
    },
) -> Tuple[float]:
    """Plot the data to compare two solvers."""

    # Remove rows where the error is equal to 0
    data = data[data["error"] != 0]

    # Store success and failure data
    success = 0
    failure = 0
    total = 0

    for (desc1, desc2), dat in data.groupby(["descriptor1", "descriptor2"]):

        total += 1

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
                # If it is a succes then there is no need to look at the
                # other chunks (for things like rate calculations)
                break
            else:
                # Make sure that none of the errors are lower than the accuracy
                assert not any(chunk["error"] < accuracy)
                failure += 1
                ax.plot(chunk["iteration"], chunk["error"], **failure_settings)
                break

    return success, failure, total


def parse_run(datadir: str, solver: str, reaction: str, ax, i, convergence_data):
    """Parse the calculation"""

    if not os.path.isfile(os.path.join(datadir, "error_log.csv")):
        logger.warning(f"File {os.path.join(datadir, 'error_log.csv')} does not exist.")
        print(f"File {os.path.join(datadir, 'error_log.csv')} does not exist.")

        return

    # Use pandas to read the csv file
    # Read the data as strings
    data = pd.read_csv(os.path.join(datadir, "error_log.csv"), dtype=str, delimiter=",", header=None)
    # Convert the strings to floats
    data = data.astype(float)
    # Set labels for data
    data.columns = ["descriptor1", "descriptor2", "iteration", "error"]

    # Get the accuracy of the solver
    with open(os.path.join(datadir, "solver_specifics.json")) as f:
        solver_specifics = json.load(f)
    ACCURACY = solver_specifics["tolerance"]
    logger.info("Accuracy: %s", ACCURACY)

    # Get the data to plot
    success, failure, total = evaulate_data(
        data,
        accuracy=ACCURACY,
        ax=ax[i],
    )

    # Store the data in the dataframe
    data_to_df = {
        "solver": solver,
        "reaction": reaction,
        "success": success,
        "failure": failure,
        "total": total,
    }

    # Concatenate the data to the dataframe
    convergence_data = pd.concat(
        [convergence_data, pd.DataFrame(data_to_df, index=[0])]
    )

    return convergence_data


def run_and_parse_calculation(
    templates,
    reaction_name: str,
    desc1: float,
    desc2: float,
    basedir: str,
    solver_specifics: Dict[str, Any],
    solvers: List[str],
    ax: plt.Axes,
    convergence_data: pd.DataFrame,
):
    """Run the single point calculation if the folder is not already present."""

    # Read the common template files for the reactions
    template = Template(templates["mkm_" + reaction_name])
    mkm_job_template = Template(templates["mkm_job"])

    # Create the input file for the reaction
    input_file = template.substitute(
        desc1=desc1,
        desc2=desc2,
    )

    mkm_job = mkm_job_template.substitute()

    # Create the directory for the reaction
    for idx, solver in enumerate(solvers):

        datadir = os.path.join(
            basedir, reaction_name, solver, f"desc1_{desc1}_desc2_{desc2}"
        )

        if not os.path.exists(datadir):
            os.makedirs(datadir)
        else:
            logger.info(f"Skipping {reaction_name} {desc1} {desc2} {solver}")
            convergence_data = parse_run(
                datadir, solver, reaction_name, ax, idx, convergence_data
            )
            continue

        # Write out a json file with the desc1 and desc2 values
        with open(os.path.join(datadir, "desc1_desc2.json"), "w") as f:
            json.dump({"desc1": desc1, "desc2": desc2}, f)

        # Write out the input file from the template
        with open(os.path.join(datadir, "input.mkm"), "w") as f:
            f.write(input_file)

        # Write out the solver specifics file
        with open(os.path.join(datadir, "solver_specifics.json"), "w") as f:
            json.dump(solver_specifics, f)

        # Write out the mkm_job.py file from the template
        with open(os.path.join(datadir, "mkm_job.py"), "w") as f:
            f.write(mkm_job)

        # Copy the energies file for the reaction to the required directory
        # Make all letters in reaction_name lowercase
        shutil.copy(
            os.path.join("input_files", f"energies_{reaction_name.lower()}.txt"),
            os.path.join(datadir, "energies.txt"),
        )

        # Execute the mkm_job.py file
        current_dir = os.getcwd()
        os.chdir(datadir)
        os.system("python mkm_job.py")
        os.chdir(current_dir)

        convergence_data = parse_run(
            datadir, solver, reaction_name, ax, idx, convergence_data
        )

    return convergence_data


if __name__ == "__main__":
    """Plot the maximum error for each individually run point."""

    descriptor_ranges = yaml.safe_load(open("config/descriptor_ranges.yaml"))

    args = get_cli_args()

    _basedir = os.path.join(os.getcwd(), "individual_runs")

    precision_folders = glob.glob(os.path.join(_basedir, "precision_*"))    

    solver_specifics = json.load(
        open(os.path.join("input_files", "solver_specifics_individual.json"))
    )

    solvers = ["coverages", "numbers_free_xstar"]

    solver_statistics = pd.DataFrame(
        columns=["precision", "solver", "reaction", "sum of success", "sum of failure"]
    )

    for basedir in precision_folders:

        precision = basedir.split("_")[-1].replace("precision_", "")
        precision = int(precision)

        for reaction_name, desc12_range in descriptor_ranges.items():

            logger.info(f"Plotting and running {reaction_name}.")

            desc1_range = np.linspace(
                desc12_range["desc1_min"], desc12_range["desc1_max"], args.resolution
            )
            desc2_range = np.linspace(
                desc12_range["desc2_min"], desc12_range["desc2_max"], args.resolution
            )

            mesh_desc1, mesh_desc2 = np.meshgrid(desc1_range, desc2_range)

            fig, ax = plt.subplots(
                1, 2, figsize=(4.5, 2.5), constrained_layout=True, sharey=True
            )

            # Create a dataframe to store the data
            convergence_data = pd.DataFrame(
                columns=["solver", "reaction", "success", "failure", "total"]
            )

            for desc1, desc2 in zip(mesh_desc1.flatten(), mesh_desc2.flatten()):

                convergence_data = run_and_parse_calculation(
                    templates,
                    reaction_name,
                    desc1,
                    desc2,
                    basedir,
                    solver_specifics,
                    solvers,
                    ax,
                    convergence_data,
                )

            # Sum up the successes of the numbers solver
            numbers_success = convergence_data.loc[
                convergence_data["solver"] == "numbers_free_xstar", "success"
            ].sum()
            # Sum up the failures of the numbers solver
            numbers_failure = convergence_data.loc[
                convergence_data["solver"] == "numbers_free_xstar", "failure"
            ].sum()
            # Sum up the successes of the coverages solver
            coverages_success = convergence_data.loc[
                convergence_data["solver"] == "coverages", "success"
            ].sum()
            # Sum up the failures of the coverages solver
            coverages_failure = convergence_data.loc[
                convergence_data["solver"] == "coverages", "failure"
            ].sum()
            # Get the total
            numbers_total = convergence_data.loc[
                convergence_data["solver"] == "numbers_free_xstar", "total"
            ].sum()
            coverages_total = convergence_data.loc[
                convergence_data["solver"] == "coverages", "total"
            ].sum()

            # Store the solver statistics
            solver_statistics = pd.concat(
                [
                    solver_statistics,
                    pd.DataFrame(
                        {
                            "precision": [precision],
                            "solver": ["numbers_free_xstar"],
                            "reaction": [reaction_name],
                            "sum of success": [numbers_success],
                            "sum of failure": [numbers_failure],
                            "total": [numbers_total],
                        }
                    ),
                    pd.DataFrame(
                        {
                            "precision": [precision],
                            "solver": ["coverages"],
                            "reaction": [reaction_name],
                            "sum of success": [coverages_success],
                            "sum of failure": [coverages_failure],
                            "total": [coverages_total],
                        }
                    ),
                ]
            )

            # Set the x and y labels
            ax[0].set_xlabel("Iteration")
            ax[1].set_xlabel("Iteration")
            ax[0].set_ylabel("Least squares error")
            ax[0].set_title("Coverages solver")
            ax[1].set_title("Numbers solver")
            # Plot y on a log scale
            ax[0].set_yscale("log")
            ax[1].set_yscale("log")
            # Set the log axis ticks on the minor ticks
            ax[0].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))
            ax[1].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))

            fig.suptitle(reaction_name.replace('_', ' '), fontsize=11)

            ax[1].plot([], [], color="tab:green", label="Success")
            ax[1].plot([], [], color="tab:red", label="Failure")
            ax[1].legend(
                loc="center left", frameon=False, bbox_to_anchor=(1.04, 1), borderaxespad=0
            )

            fig.savefig(os.path.join(basedir, reaction_name, "convergence.png"), dpi=300)
            convergence_data.to_csv(
                os.path.join(basedir, reaction_name, "convergence_data.csv")
            )
            plt.close(fig)
        
    # Sort the dataframe with precision
    solver_statistics = solver_statistics.sort_values(by="precision")

    solver_statistics.to_latex(
        os.path.join("solver_statistics.tex"),
        index=False,
        float_format="{:0.2f}".format,
    )
    solver_statistics.to_markdown(
        os.path.join("solver_statistics.md"), index=False
    )

    # Plot the precision vs success rate for each reaction
    unique_reactions = solver_statistics["reaction"].unique()
    fig, ax = plt.subplots(1, len(unique_reactions), figsize=(8, 2), constrained_layout=True)
    for idx, reaction_name in enumerate(descriptor_ranges.keys()):
        numbers_success = solver_statistics.loc[
            (solver_statistics["solver"] == "numbers_free_xstar")
            & (solver_statistics["reaction"] == reaction_name),
            "sum of success",
        ]
        coverage_success = solver_statistics.loc[
            (solver_statistics["solver"] == "coverages")
            & (solver_statistics["reaction"] == reaction_name),
            "sum of success",
        ]
        precision = solver_statistics.loc[
            (solver_statistics["solver"] == "numbers_free_xstar")
            & (solver_statistics["reaction"] == reaction_name),
            "precision",
        ]
        ax[idx].plot(precision, numbers_success, '-o', color="tab:blue", label="Numbers")
        ax[idx].plot(precision, coverage_success, '-o', color="tab:orange", label="Coverages")
        ax[idx].set_title(reaction_name.replace('_', ' '), fontsize=11)

    for _ax in ax:
        _ax.set_xlabel("Precision")
        _ax.set_ylim(-5, 105)
    ax[0].set_ylabel("Successful runs")
    ax[-1].legend(loc="center left", frameon=False, bbox_to_anchor=(1.04, 1), borderaxespad=0)

    fig.savefig(os.path.join("outputs/precision_vs_success_rate.png"), dpi=300)
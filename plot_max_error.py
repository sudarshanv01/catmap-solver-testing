import os

import shutil

import json

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


def parse_run(datadir: str, solver: str, reaction: str, ax, i, convergence_data):
    """Parse the calculation"""

    if not os.path.isfile(os.path.join(datadir, "error_log.csv")):
        logger.warning(f"File {os.path.join(datadir, 'error_log.csv')} does not exist.")

        return

    # Use pandas to read the csv file
    data = pd.read_csv(os.path.join(datadir, "error_log.csv"))
    # Set labels for data
    data.columns = ["descriptor1", "descriptor2", "iteration", "error"]

    # Get the accuracy of the solver
    with open(os.path.join(datadir, "solver_specifics.json")) as f:
        solver_specifics = json.load(f)
    ACCURACY = solver_specifics["tolerance"]
    logger.info("Accuracy: %s", ACCURACY)

    # Get the data to plot
    success, failure = evaulate_data(
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

    basedir = os.path.join(os.getcwd(), "individual_runs")

    solver_specifics = json.load(
        open(os.path.join("input_files", "solver_specifics_individual.json"))
    )

    solvers = ["coverages", "numbers_free_xstar"]

    for reaction_name, desc12_range in descriptor_ranges.items():
        logger.info(f"Plotting and running {reaction_name}.")

        desc1_range = np.linspace(
            desc12_range["desc1_min"], desc12_range["desc1_max"], args.resolution
        )
        desc2_range = np.linspace(
            desc12_range["desc2_min"], desc12_range["desc2_max"], args.resolution
        )

        mesh_desc1, mesh_desc2 = np.meshgrid(desc1_range, desc2_range)

        fig, ax = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

        # Create a dataframe to store the data
        convergence_data = pd.DataFrame(
            columns=["solver", "reaction", "success", "failure"]
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
        print(convergence_data)

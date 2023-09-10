import argparse
import copy
import glob
import json
import logging
import os
import shutil
from string import Template
from typing import Any, Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml
from matplotlib.ticker import LogLocator

from plot_params import get_plot_params
from template_reactions import templates

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
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
    data = pd.read_csv(
        os.path.join(datadir, "error_log.csv"), dtype=str, delimiter=",", header=None
    )
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


def compute_success_failure_total(
    solver_name: str, df: pd.DataFrame
) -> Tuple[int, int, int]:
    """Compute the number of successful and failed runs for a given solver."""

    # Sum up the successes of the numbers solver
    success = df.loc[df["solver"] == solver_name, "success"].sum()
    # Sum up the failures of the numbers solver
    failure = df.loc[df["solver"] == solver_name, "failure"].sum()
    # Get the total
    total = df.loc[df["solver"] == solver_name, "total"].sum()

    return success, failure, total


def create_concat_df(
    precision: float, solver: str, reaction: str, success: int, failure: int, total: int
) -> pd.DataFrame:
    """Create a dataframe with the solver statistics."""
    return pd.DataFrame(
        {
            "precision": [precision],
            "solver": [solver],
            "reaction": [reaction],
            "sum of success": [success],
            "sum of failure": [failure],
            "total": [total],
        }
    )


if __name__ == "__main__":
    """Plot the maximum error for each individually run point."""
    descriptor_ranges = yaml.safe_load(open("config/descriptor_ranges.yaml"))
    args = get_cli_args()
    _basedir = os.path.join(os.getcwd(), "individual_runs")
    precision_folders = glob.glob(os.path.join(_basedir, "precision_*"))
    _solver_specifics = json.load(
        open(os.path.join("input_files", "solver_specifics_individual.json"))
    )
    solvers = ["coverages", "numbers_free_xstar", "numbers_fix_xstar"]
    solver_statistics = pd.DataFrame(
        columns=["precision", "solver", "reaction", "sum of success", "sum of failure"]
    )

    for basedir in precision_folders:
        solver_specifics = copy.deepcopy(_solver_specifics)
        precision = basedir.split("_")[-1].replace("precision_", "")
        precision = int(precision)
        solver_specifics["decimal_precision"] = precision
        for reaction_name, desc12_range in descriptor_ranges.items():
            logger.info(f"Plotting and running {reaction_name}.")
            desc1_range = np.linspace(
                desc12_range["desc1_min"],
                desc12_range["desc1_max"],
                num=args.resolution,
            )
            desc2_range = np.linspace(
                desc12_range["desc2_min"],
                desc12_range["desc2_max"],
                num=args.resolution,
            )
            mesh_desc1, mesh_desc2 = np.meshgrid(desc1_range, desc2_range)
            fig, ax = plt.subplots(
                1, 3, figsize=(4.5, 2.5), constrained_layout=True, sharey=True
            )
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
            numbers_free_x_stats = compute_success_failure_total(
                "numbers_free_xstar", convergence_data
            )
            numbers_free_x_star_success = numbers_free_x_stats[0]
            numbers_free_x_star_failure = numbers_free_x_stats[1]
            numbers_free_x_star_total = numbers_free_x_stats[2]
            numbers_free_x_star_df = create_concat_df(
                precision,
                "numbers_free_xstar",
                reaction_name,
                numbers_free_x_star_success,
                numbers_free_x_star_failure,
                numbers_free_x_star_total,
            )
            numbers_fix_x_stats = compute_success_failure_total(
                "numbers_fix_xstar", convergence_data
            )
            numbers_fix_x_star_success = numbers_fix_x_stats[0]
            numbers_fix_x_star_failure = numbers_fix_x_stats[1]
            numbers_fix_x_star_total = numbers_fix_x_stats[2]
            numbers_fix_x_star_df = create_concat_df(
                precision,
                "numbers_fix_xstar",
                reaction_name,
                numbers_fix_x_star_success,
                numbers_fix_x_star_failure,
                numbers_fix_x_star_total,
            )
            coverages_stats = compute_success_failure_total(
                "coverages", convergence_data
            )
            coverages_success = coverages_stats[0]
            coverages_failure = coverages_stats[1]
            coverages_total = coverages_stats[2]
            coverages_df = create_concat_df(
                precision,
                "coverages",
                reaction_name,
                coverages_success,
                coverages_failure,
                coverages_total,
            )
            solver_statistics = pd.concat(
                [
                    solver_statistics,
                    numbers_free_x_star_df,
                    numbers_fix_x_star_df,
                    coverages_df,
                ]
            )

            ax[0].set_xlabel("Iteration")
            ax[1].set_xlabel("Iteration")
            ax[0].set_ylabel("Least squares error")
            for _idx_ax, _solver in enumerate(solvers):
                ax[_idx_ax].set_title(_solver.replace("_", " "))
            ax[0].set_yscale("log")
            ax[1].set_yscale("log")
            ax[0].set_xlim(0, 20)
            ax[1].set_xlim(0, 20)
            ax[0].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))
            ax[1].yaxis.set_minor_locator(LogLocator(base=10.0, subs="all"))
            fig.suptitle(reaction_name.replace("_", " "), fontsize=11)
            ax[-1].plot([], [], color="tab:green", label="Success")
            ax[-1].plot([], [], color="tab:red", label="Failure")
            ax[-1].legend(
                loc="center left",
                frameon=False,
                bbox_to_anchor=(1.04, 1),
                borderaxespad=0,
            )
            fig.savefig(
                os.path.join(basedir, reaction_name, "convergence.png"), dpi=300
            )
            convergence_data.to_csv(
                os.path.join(basedir, reaction_name, "convergence_data.csv")
            )
            plt.close(fig)

    solver_statistics = solver_statistics.sort_values(by="precision")
    solver_statistics = solver_statistics.reset_index(drop=True)
    solver_statistics.to_latex(
        os.path.join("solver_statistics.tex"),
        index=False,
        float_format="{:0.2f}".format,
    )
    solver_statistics.to_markdown(os.path.join("solver_statistics.md"), index=False)
    fig, ax = plt.subplots(
        1, len(descriptor_ranges), figsize=(8, 2), constrained_layout=True
    )
    for idx, reaction_name in enumerate(descriptor_ranges.keys()):
        df_reaction = solver_statistics.loc[
            solver_statistics["reaction"] == reaction_name
        ]
        df_reaction = df_reaction.sort_values(by="solver")
        sns.lineplot(
            x="precision",
            y="sum of success",
            hue="solver",
            data=df_reaction,
            ax=ax[idx],
            markers=True,
            style="solver",
            dashes=False,
            alpha=0.5,
        )
        ax[idx].set_title(reaction_name.replace("_", " "))
    for _ax in ax:
        _ax.set_xlabel("Precision")
        _ax.set_ylim(-5, 105)
    ax[0].set_ylabel("Successful runs")
    ax[-1].legend(
        loc="center left", frameon=False, bbox_to_anchor=(1.04, 1), borderaxespad=0
    )
    for _ax in ax[:-1]:
        _ax.get_legend().remove()
    for _ax in ax[1:]:
        _ax.set_ylabel("")
    fig.savefig(os.path.join("outputs/precision_vs_success_rate.png"), dpi=300)

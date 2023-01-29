import os

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


def get_cli_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--resolution", type=int, default=10, help="Number of points to plot."
    )
    return parser.parse_args()


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

        for desc1, desc2 in zip(mesh_desc1.flatten(), mesh_desc2.flatten()):

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
            for solver in solvers:
                datadir = os.path.join(
                    basedir, reaction_name, solver, f"desc1_{desc1}_desc2_{desc2}"
                )
                if not os.path.exists(datadir):
                    os.makedirs(datadir)

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

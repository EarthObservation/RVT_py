"""
This module contains functions to create requirements.txt (pip) and environment.yaml (conda) from conda environment.

Anaconda system variable should be configured.

To create conda env use:
 - "conda env create -f environment.yml" or
 - "conda create --name <env_name> --file requirements.txt"

To install dependencies with pip:
 - "pip install -r requirements.txt"

In development process the dependencies are added/updated/removed and on every change project's requirements.txt and
environment.yaml should also be updated. Updating (overwriting) those is purpose of this module.
"""
import os
from pathlib import Path

CONDA_ENV_NAME = "rvt39"

def create_requirements_txt() -> None:
    # remove existing file
    requirements_path = Path(__file__).parent.parent / "requirements.txt"
    if requirements_path.exists():
        os.remove(requirements_path)

    # activate conda env and create file
    os.system(f"activate {CONDA_ENV_NAME} && pip freeze > {requirements_path}")

    if not requirements_path.exists():
        raise FileNotFoundError(f"File {requirements_path} doesn't exist, something went wrong!")


def create_environment_yaml() -> None:
    # remove existing file
    environment_path = Path(__file__).parent.parent / "environment.yml"
    if environment_path.exists():
        os.remove(environment_path)

    # activate conda env and create file
    os.system(f"activate {CONDA_ENV_NAME} && conda env export > {environment_path}")

    if not environment_path.exists():
        raise FileNotFoundError(f"File {environment_path} doesn't exist, something went wrong!")


if __name__ == '__main__':
    create_requirements_txt()
    create_environment_yaml()

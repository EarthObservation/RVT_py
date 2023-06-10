"""
In development process the dependencies are added/updated/removed and on every dependency change project's
requirements.txt and environment.yaml should also be updated.
If dependencies were added/removed/updated run this module to update environment.yaml (conda) and requirements.txt (pip)
from conda environment.
These two files define dependencies and versions of python packages that are needed and used for running rvt and
development of rvt.

Anaconda system variable should be configured.

To create conda env use:
 - "conda env create -f environment.yml" or
 - "conda create --name <env_name> --file requirements.txt"

To install dependencies with pip:
 - "pip install -r requirements.txt"
"""
import os
from dev_tools import CONDA_ENV_NAME, CONDA_ENV_YAML_PATH, REQUIREMENTS_TXT_PATH


def create_requirements_txt() -> None:
    # remove existing file
    if REQUIREMENTS_TXT_PATH.exists():
        os.remove(REQUIREMENTS_TXT_PATH)

    # activate conda env and create file
    os.system(f"activate {CONDA_ENV_NAME} && pip freeze > {REQUIREMENTS_TXT_PATH}")

    if not REQUIREMENTS_TXT_PATH.exists():
        raise FileNotFoundError(f"File {REQUIREMENTS_TXT_PATH} doesn't exist, something went wrong!")


def create_conda_environment_yaml() -> None:
    # remove existing file
    if CONDA_ENV_YAML_PATH.exists():
        os.remove(CONDA_ENV_YAML_PATH)

    # activate conda env and create file
    os.system(f"activate {CONDA_ENV_NAME} && conda env export > {CONDA_ENV_YAML_PATH}")

    if not CONDA_ENV_YAML_PATH.exists():
        raise FileNotFoundError(f"File {CONDA_ENV_YAML_PATH} doesn't exist, something went wrong!")


if __name__ == '__main__':
    create_requirements_txt()
    create_conda_environment_yaml()

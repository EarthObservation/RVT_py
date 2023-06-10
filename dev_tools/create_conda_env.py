"""
To start working on the project you need to set up the anaconda environment from environemnt.yaml file
(which defines dependencies).

Anaconda system variable should be configured.
"""
import os
from dev_tools import CONDA_ENV_YAML_PATH, CONDA_ENV_NAME


def create_conda_environment_from_yaml() -> None:
    if not CONDA_ENV_YAML_PATH.exists():
        raise FileNotFoundError(f"File {CONDA_ENV_YAML_PATH} doesn't exist, something went wrong!")

    os.system(
        f"conda env create -f {CONDA_ENV_YAML_PATH}"
    )


if __name__ == '__main__':
    create_conda_environment_from_yaml()

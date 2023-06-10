"""
Project at some state (branch, commit hash) can depend on different python libraries. Therefore, before start working
on the project it is necessary to update the conda environment (defined in environment.yaml) corresponding to that
state.

Anaconda system variable should be configured.
"""
import os
from dev_tools import CONDA_ENV_YAML_PATH, CONDA_ENV_NAME

def update_conda_environment_from_yaml() -> None:
    if not CONDA_ENV_YAML_PATH.exists():
        raise FileNotFoundError(f"File {CONDA_ENV_YAML_PATH} doesn't exist, something went wrong!")

    os.system(
        f"conda env update --name {CONDA_ENV_NAME} --file {CONDA_ENV_YAML_PATH} --prune"
    )


if __name__ == '__main__':
    update_conda_environment_from_yaml()

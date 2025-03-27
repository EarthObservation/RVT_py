import subprocess
import shutil
import os

"""
Versions are built from PyPi rvt-py.
First from PyPi we build win-64 for each python version.
Then from win-64 version we create versions for other platforms (for each python version).
"""

BUILDS_DIR = "./rvt-py_builds"  # Builds for each python version
PLATFORM_BUILDS_DIR = "./rvt-py_platforms"  # Builds for each python version and platform
PYTHON_VERSIONS_TO_BUILD = [3.8, 3.9, 3.10, 3.11]
CONDA_CLOUD_USERNAME = "rvtpy"
CONDA_CLOUD_PASSWORD = "MY_PASSWORD"  # do not commit real password!

# subprocess.run("activate base", shell=True, check=True)

shutil.rmtree(BUILDS_DIR, ignore_errors=True)
shutil.rmtree(PLATFORM_BUILDS_DIR, ignore_errors=True)
os.mkdir(BUILDS_DIR)
os.mkdir(PLATFORM_BUILDS_DIR)

for python_version in PYTHON_VERSIONS_TO_BUILD:
    subprocess.run(f"conda build rvt-py --python {python_version} --output-folder {BUILDS_DIR}", shell=True, check=True)


builds = [f"{BUILDS_DIR}/win-64/{file}" for file in os.listdir(BUILDS_DIR + "/win-64") if file.endswith(".conda")]

for build in builds:
    subprocess.run(
        "conda convert --platform all {} -o {}".format(build, PLATFORM_BUILDS_DIR), shell=True, check=True
    )
    if not os.path.exists(f"{PLATFORM_BUILDS_DIR}/win-64"):
        os.mkdir(f"{PLATFORM_BUILDS_DIR}/win-64")
    shutil.copy(build, f"{PLATFORM_BUILDS_DIR}/win-64/{os.path.basename(build)}")

platform_version_builds = os.listdir(PLATFORM_BUILDS_DIR)

subprocess.run(
    "anaconda login --username {} --password {}".format(CONDA_CLOUD_USERNAME, CONDA_CLOUD_PASSWORD),
    shell=True,
    check=True
)

for platform_version_build in platform_version_builds:
    platform_version_build_path = os.path.join(PLATFORM_BUILDS_DIR, platform_version_build)
    file_name = os.listdir(platform_version_build_path)[0]
    subprocess.run(
        "anaconda upload {}".format(os.path.join(platform_version_build_path, file_name)), shell=True, check=True
    )

shutil.rmtree(BUILDS_DIR, ignore_errors=True)
shutil.rmtree(PLATFORM_BUILDS_DIR, ignore_errors=True)

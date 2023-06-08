"""
Run this module to do mypy check.
"""

import os
from pathlib import Path


def mypy_check() -> int:
    rvt_library_dir_path = Path(__file__).parent.parent / "rvt"
    exit_code = os.system(f"mypy {rvt_library_dir_path}")
    if exit_code != 0:
        raise Exception


if __name__ == '__main__':
    mypy_check()

"""
Run this module to do black formatting.
"""

import os
from pathlib import Path


def black_formatting() -> int:
    rvt_library_dir_path = Path(__file__).parent.parent / "rvt"
    exit_code = os.system(f"black {rvt_library_dir_path}")
    if exit_code != 0:
        raise Exception


if __name__ == '__main__':
    black_formatting()

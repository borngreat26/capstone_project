"""Shared helpers: repository paths and the PLINK2 executable."""
import os
import shutil
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def plink2_path():
    """Find PLINK2: the PLINK2_BIN variable, then plink2 on PATH, then a local ./plink2_binary."""
    configured = os.environ.get('PLINK2_BIN')
    if configured:
        if not Path(configured).is_file():
            raise FileNotFoundError(f'PLINK2_BIN points to a missing file: {configured}')
        return configured
    on_path = shutil.which('plink2')
    if on_path:
        return on_path
    local = REPO_ROOT / 'plink2_binary'
    if local.is_file():
        return str(local)
    raise FileNotFoundError('PLINK2 was not found. Install it (see README.md) or set PLINK2_BIN to its path.')


def run_plink(arguments):
    """Run PLINK2 with the project's thread and memory limits; show its output if it fails."""
    command = [plink2_path(), *map(str, arguments), '--threads', '2', '--memory', '1024']
    finished = subprocess.run(command, capture_output=True, text=True)
    if finished.returncode != 0:
        raise RuntimeError(f"PLINK2 failed: {' '.join(command)}\n{finished.stdout}\n{finished.stderr}")

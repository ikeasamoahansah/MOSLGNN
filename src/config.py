import os
from pathlib import Path

ROOT_DIR = os.environ.get('', None)

if ROOT_DIR is not None:
    ROOT_DIR = Path(ROOT_DIR)

if ROOT_DIR is None or not ROOT_DIR.exists() or not ROOT_DIR.is_dir():
    ROOT_DIR = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT_DIR / 'data'
RESULTS_DIR = ROOT_DIR / 'results'
SEED = 42

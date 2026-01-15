"""
Test script to verify all imports and basic functionality.
Run this before using the main scripts to check if all dependencies are installed.
"""

import sys
from pathlib import Path

print("="*70)
print("Testing AG_BenchMark_Builder Dependencies")
print("="*70)

# Test basic imports
print("\n1. Testing basic imports...")
try:
    import numpy as np
    print("   numpy: OK")
except ImportError as e:
    print(f"   numpy: FAILED - {e}")

try:
    import matplotlib.pyplot as plt
    print("   matplotlib: OK")
except ImportError as e:
    print(f"   matplotlib: FAILED - {e}")

try:
    import netCDF4
    print("   netCDF4: OK")
except ImportError as e:
    print(f"   netCDF4: FAILED - {e}")

try:
    import pandas as pd
    print("   pandas: OK")
except ImportError as e:
    print(f"   pandas: FAILED - {e}")

try:
    from sklearn.mixture import BayesianGaussianMixture
    print("   scikit-learn: OK")
except ImportError as e:
    print(f"   scikit-learn: FAILED - {e}")

try:
    from scipy.optimize import minimize
    print("   scipy: OK")
except ImportError as e:
    print(f"   scipy: FAILED - {e}")

try:
    from tqdm import tqdm
    print("   tqdm: OK")
except ImportError as e:
    print(f"   tqdm: FAILED - {e}")

# Test library imports
print("\n2. Testing library imports...")
sys.path.insert(0, str(Path(__file__).parent / 'library'))

try:
    from goci2_reader import read_goci2, get_flag_by_name, GOCI2Data
    print("   goci2_reader: OK")
except ImportError as e:
    print(f"   goci2_reader: FAILED - {e}")

# Test MargModel import
print("\n3. Testing MargModel...")
sys.path.insert(0, str(Path(__file__).parent / 'MargModel'))

try:
    import MargModel
    print("   MargModel: OK")
except ImportError as e:
    print(f"   MargModel: FAILED - {e}")

# Check folder structure
print("\n4. Checking folder structure...")
base_dir = Path(__file__).parent

folders = ['library', 'MargModel', 'data', 'nfis', 'results']
for folder in folders:
    folder_path = base_dir / folder
    if folder_path.exists():
        print(f"   {folder}/: OK")
    else:
        print(f"   {folder}/: MISSING")

# Check for data files
print("\n5. Checking for data files...")
goci2_dir = base_dir / 'data' / 'goci2'
nfis_dir = base_dir / 'data' / 'nfis'

nc_files = list(goci2_dir.glob('*.nc')) if goci2_dir.exists() else []
if nc_files:
    print(f"   Found {len(nc_files)} GOCI-II NetCDF files in data/goci2/")
else:
    print(f"   WARNING: No GOCI-II NetCDF files found in data/goci2/")
    print(f"   Please place GOCI-II L2 data files in: {goci2_dir}")

nfis_files = list(nfis_dir.glob('*.csv')) + list(nfis_dir.glob('*.xlsx')) if nfis_dir.exists() else []
if nfis_files:
    print(f"   Found {len(nfis_files)} NFIS files in data/nfis/")
else:
    print(f"   Note: No NFIS files found in data/nfis/ (optional)")

print("\n" + "="*70)
print("Test Complete!")
print("="*70)


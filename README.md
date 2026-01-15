# 🌊 AG Benchmark Builder

**Benchmark Builder for Algal Bloom Detection using GOCI-II Satellite Data**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## 📋 Overview

A comprehensive toolkit for building algal bloom benchmark datasets from GOCI-II (Geostationary Ocean Color Imager-II) satellite data. This package provides interactive GUI tools for region selection, clustering analysis, and bio-optical model-based fluorescence correction.

### Main Features

1. **P01 - Region Selector**: Interactive GUI for selecting and cropping regions from GOCI-II imagery
2. **P02 - Cluster Selector**: Online clustering with global cluster database management
3. **P03 - cFLH Calculator**: Corrected Fluorescence Line Height (cFLH) using bio-optical modeling

## Structure
```
AG_BenchMark_Builder/
├── library/                           # GOCI-II reader library
│   ├── __init__.py
│   └── goci2_reader.py               # NetCDF reader + flag handling
├── MargModel/                         # Bio-optical model (CCRR-RT)
│   ├── __init__.py
│   └── MargModel.py
├── data/
│   ├── goci2/                         # GOCI-II NetCDF files (.nc)
│   └── nfis/                          # NFIS observation data (CSV/Excel)
├── results/                           # Output folder
│   └── GUI_archive/
│       ├── regions.json               # Saved regions
│       ├── cropped_data/              # Cropped NPZ files
│       ├── clustering_results/        # Local clustering results
│       ├── cluster_archive/           # Archived clusters
│       ├── global_clusters.pkl        # Global cluster database
│       └── scene_results.json         # Scene processing log
├── P01_GUI_region_selector.py        # Region selection GUI
├── P02_GUI_cluster_selector.py       # Clustering GUI  
├── P03_add_cFLH_to_npz.py            # Add cFLH to archived clusters
├── test_imports.py                    # Dependency test script
├── requirements.txt                   # Python dependencies
└── README.md                          # This file
```

## Installation

### 1. Install Python dependencies
```bash
pip install -r requirements.txt
```

Or manually:
```bash
pip install numpy matplotlib netCDF4 pandas scikit-learn scipy openpyxl tqdm
```

### 2. Test installation
```bash
python test_imports.py
```

### 3. Prepare data
- Place GOCI-II L2 NetCDF files (`.nc`) in `data/goci2/` folder
- (Optional) Place NFIS observation files (`.csv` or `.xlsx`) in `data/nfis/` folder

## Usage

### P01: Region Selector
Interactive GUI for selecting and cropping rectangular regions from GOCI-II data.

**Features:**
- Display GOCI-II RhoC RGB composite
- Draw rectangular regions with mouse
- Save region definitions to JSON
- Load and display NFIS observation points
- Crop regions to NPZ files (includes RhoC, Rrs, FLH, lat/lon, flags)

**Usage:**
```bash
python P01_GUI_region_selector.py
```

**Workflow:**
1. Select background GOCI-II file from dropdown
2. Draw rectangle with mouse drag
3. Click "Save" to name and save region
4. Click "Crop" to export selected region
5. (Optional) Click "Add NFIS" to load observation points

**Output:**
- `results/GUI_archive/regions.json` - Region definitions
- `results/GUI_archive/cropped_data/*.npz` - Cropped data files

---

### P02: Cluster Selector
Interactive GUI for performing online clustering and managing global cluster database.

**Features:**
- Load cropped GOCI-II data
- Perform online clustering (Bayesian GMM + spatial splitting)
- View clusters colored by mean FLH
- Interactive Rrs spectrum display
- Archive selected clusters
- Manage global cluster database across multiple scenes

**Usage:**
```bash
python P02_GUI_cluster_selector.py
```

**Workflow:**
1. Select date from dropdown (extracted from cropped files)
2. Select file from list
3. Click "Clustering" to perform online clustering
4. Click on clusters to view Rrs spectra
5. Use "Hold" checkbox to compare multiple clusters
6. Click "Archive" to save selected clusters
7. Use "Display Archived" to view previously archived clusters
8. Use "Delete Cluster" to remove clusters from archive

**Clustering Algorithm:**
- **Feature extraction**: MNF (4 components) + NDVI + FLH
- **Spectral clustering**: Bayesian GMM (K=50, Dirichlet Process prior)
- **Spatial splitting**: Connected component analysis
- **Global matching**: Spectral Angle Mapper + Symmetric KL divergence
- **Online update**: Exponentially weighted moving average (ρ=0.1)

**Output:**
- `results/GUI_archive/clustering_results/` - Local clustering results
- `results/GUI_archive/global_clusters.pkl` - Global cluster database
- `results/GUI_archive/scene_results.json` - Processing log
- `results/GUI_archive/cluster_archive/clusters_YYYYMMDD.npz` - Archived clusters
- `results/GUI_archive/cluster_archive/clusters_YYYYMMDD.txt` - Summary text

---

### P03: Add cFLH to NPZ
Batch processing script to calculate corrected FLH (cFLH) for all archived cluster pixels.

**Features:**
- Bio-optical model inversion (MargModel/CCRR-RT)
- Fluorescence modeling
- Mineral-corrected FLH calculation
- Peak wavelength detection (660-745nm range)
- Progress tracking with tqdm

**Usage:**
```bash
python P03_add_cFLH_to_npz.py
```

**cFLH Calculation Method:**
1. Optimize water quality parameters (Chla, Mineral, CDOM) to match observed Rrs
2. Generate hyperspectral Rrs (400-750nm, 0.5nm resolution) with fluorescence
3. Recalculate Rrs with mineral=0 (corrected)
4. Find peak between 660-745nm from corrected Rrs
5. Calculate cFLH as peak height above baseline (660-745nm)

**Optimization Modes:**
- **relative** (default): Match spectral slopes + absolute at 709/745nm
- **strict**: Match absolute Rrs values directly

**Note:** This process is computationally intensive (~20-30 seconds per pixel). For 2665 pixels, expect ~15-20 hours.

**Output:**
- `results/GUI_archive/cluster_archive/clusters_YYYYMMDD_cFLH.npz` - Original data + cFLH values

---

## Data Format

### Cropped NPZ Files
```python
data = np.load('cropped_file.npz')
# Contains:
# - RhoC_380, RhoC_412, ..., RhoC_865  (12 bands)
# - Rrs_380, Rrs_412, ..., Rrs_865    (12 bands)
# - FLH_RhoC, FLH_Rrs                  (FLH from both)
# - latitude, longitude                 (Geolocation)
# - flag                                (Quality flags)
```

### Archived Cluster NPZ Files
```python
data = np.load('clusters_20220909.npz', allow_pickle=True)
clusters = data['archived_clusters'].item()  # Dictionary

for cluster_id, cluster_info in clusters.items():
    # cluster_info contains:
    # - global_cluster_id: int
    # - n_pixels: int
    # - latitude, longitude: arrays
    # - rows, cols: arrays (pixel coordinates)
    # - Rrs_380, ..., Rrs_865: arrays (12 bands)
    # - RhoC_380, ..., RhoC_865: arrays (12 bands)
    # - FLH: array
    # - cFLH: array (if P03 has been run)
    # - source_file, date, timestamp: metadata
```

## Dependencies

- **numpy**: Array operations
- **matplotlib**: Visualization
- **netCDF4**: GOCI-II file reading
- **pandas**: NFIS data handling
- **scikit-learn**: Bayesian GMM clustering
- **scipy**: Optimization, linear algebra
- **openpyxl**: Excel file support
- **tqdm**: Progress bars

## Notes

1. **FLH Baseline**: All FLH calculations use 660-745nm baseline (not 680-745nm)
2. **Coordinate System**: All pixel coordinates are absolute (relative to original image)
3. **Clustering**: Global clusters persist across scenes and dates
4. **Performance**: P03 is slow (~20s/pixel). Consider running overnight for large datasets.
5. **Data Quality**: P03 filters out pixels with optimization failures (returns NaN)

## 📊 Example Dataset

This repository includes sample data:
- **GOCI-II L2 data**: September 9, 2022 scene (example file)
- **NFIS red tide observations**: Matched observation points

## 📚 Citation

If you use this tool, please cite:
- **GOCI-II data**: Korea Ocean Satellite Center (KOSC)
- **NFIS data**: National Institute of Fisheries Science (NIFS), Korea

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details

## 👨‍💻 Author

**wonk19**

## 🙏 Acknowledgments

- Korea Ocean Satellite Center (KOSC) for GOCI-II data
- National Institute of Fisheries Science (NIFS) for red tide observation data

---

**Note**: This tool is designed for research purposes. For large-scale processing, consider optimizing P03 (cFLH calculation) which is computationally intensive.


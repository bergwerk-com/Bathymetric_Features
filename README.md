# Bathymetric Feature Delineation and Analysis

A complete workflow for detecting and analyzing bathymetric features using topographic prominence analysis and multi-level contour extraction.

## Overview

This pipeline processes bathymetric data to:
1. Identify bathymetric peaks and calculate their topographic prominence
2. Extract multi-level contours around each feature
3. Generate GeoPackage outputs for visualization in QGIS

## Repository Structure

```
knoll/
├── code/                           # Python processing scripts
│   ├── mountains/                  # Git submodule (not included in repo)
│   ├── organize_prominence_results.py
│   ├── find_ctrs.py
│   └── ctrs_cluster.py
├── data/                           # Data directory (input/output)
│   ├── input/                      # Input bathymetry files (user-provided)
│   └── output/                     # Processing results
├── environment.yml                 # Conda environment specification
├── run_GEBCO_2014_30s.sh          # GEBCO 2014 workflow
└── run_GEBCO_2025_15s.sh          # GEBCO 2025 workflow
```

**Note:** The `code/mountains/` directory is managed as a Git submodule and points to the [Mountains](https://github.com/akirmse/mountains) repository by Adam Kirmse.

## Quick Start

For users cloning this repository for the first time:

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/bergwerk-com/Bathymetric_Features.git
cd Bathymetric_Features

# 2. Set up conda environment
conda env create -f environment.yml
conda activate BathyFeatures

# 3. Compile Mountains binary (here for OSX and gcc, see [Mountains](https://github.com/akirmse/mountains) repository for Windows)
cd code/mountains/code && make
cd ../../../../  # Return to repository root

# 4. Run a workflow
./run_GEBCO_2014_30s.sh
```

## Installation

### 1. Clone the Repository with Submodules

To clone this repository along with the Mountains prominence calculation tool, use:

```bash
# Clone the repository and initialize the mountains submodule in one command
git clone --recurse-submodules https://github.com/bergwerk-com/Bathymetric_Features.git

# Or if you've already cloned without --recurse-submodules:
cd Bathymetric_Features
git submodule update --init --recursive
```

This will automatically clone the [Mountains](https://github.com/akirmse/mountains) repository into `code/mountains/`.

### 2. Install Software Dependencies

The easiest way to install all dependencies is using conda:

```bash
# Create the environment from the provided file
conda env create -f environment.yml

# Activate the environment
conda activate BathyFeaturesCtrs
```

The `environment.yml` file includes:
- Python 3.11 with geospatial packages:
  - `geopandas` - Geospatial data manipulation
  - `rasterio` - Raster I/O with GDAL
  - `gdal` - GDAL command-line tools and Python bindings
  - `contourpy` - Fast contour generation
  - `shapely` - Geometric operations
  - `pyproj` - Coordinate transformations
  - `numpy` - Numerical operations
  - `pandas` - Data manipulation
  - `xarray` - NetCDF file handling
  - `netcdf4` - NetCDF backend for xarray
  - `tqdm` - Progress bars

**Additional requirement:**
- Mountains prominence calculation binary (see [Compile Mountains Binary](#3-compile-mountains-binary) below)

### System Requirements
- 64GB+ RAM recommended for large datasets
- Multi-core CPU (8+ cores recommended)
- Sufficient disk space for intermediate files (~3x input data size)
- CMake 3.10+ and C++ compiler (for building Mountains binary)

### 3. Compile Mountains Binary

After cloning with submodules, compile the Mountains prominence calculation tool:

```bash
cd code/mountains/code
mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j$(nproc)
```

The compiled binary will be available at `code/mountains/code/release/`, which is where the bash scripts expect it.

**Verify the installation:**

```bash
# From the release directory
./divide_tree --help

# Or from the repository root
code/mountains/code/release/divide_tree --help
```

You should see the help message for the Mountains prominence calculation tool.

**About the Mountains Tool:**

This pipeline uses the [Mountains](https://github.com/akirmse/mountains) prominence calculation tool by Adam Kirmse for detecting topographic peaks and calculating their prominence.

- **Repository**: https://github.com/akirmse/mountains
- **Documentation**: See the Mountains repository README
- **Citation**: Kirmse, A. (2017). "Calculating the prominence and isolation of every mountain in the world"

**Alternative: Use Existing Installation**

If you already have Mountains compiled elsewhere, you can either:
1. Create a symbolic link: `ln -s /path/to/your/mountains code/mountains`
2. Modify the `--binary_dir` parameter in the bash scripts to point to your installation

## Pipeline Workflow

### Step 1: Bathymetry Preprocessing

Convert and filter bathymetric data to appropriate format:

```bash
# For NetCDF input (GEBCO datasets)
gdal_translate -of EHdr -ot Float32 NETCDF:"input.nc":elevation output_full.flt

# Filter to bathymetry only (values <= 0)
gdal_calc.py -A output_full.flt \
    --outfile=output_bathy.flt \
    --calc="A*(A<=0)" \
    --NoDataValue=-32768 \
    --format=EHdr \
    --overwrite

# Build pyramids for visualization
gdaladdo output_bathy.flt 2 4 8 16 32
```

For projected data (e.g., EPSG:25833), reproject to WGS84 (EPSG:4326):

```bash
gdalwarp -s_srs EPSG:25833 -t_srs EPSG:4326 \
    -tr 0.001411784399965 0.000470594799988 \
    -r bilinear input.tif output_wgs84.tif
```

### Step 2: Peak Detection & Prominence Calculation

Calculate topographic prominence using the Mountains algorithm:

```bash
python code/mountains/scripts/run_prominence.py \
    --binary_dir code/mountains/code/release \
    --threads 18 \
    --degrees_per_tile 1 \
    --samples_per_tile 120 \
    --skip_boundary \
    --min_prominence 300 \
    --bathymetry \
    input_bathy.flt
```

**Key Parameters:**
- `--samples_per_tile`: Grid points per degree (120 for GEBCO 30s, 240 for GEBCO 15s)
- `--min_prominence`: Minimum prominence threshold in meters
- `--bathymetry`: Flag for underwater features

### Step 3: Organize Results, Filter and Relocate Peaks

Sort, filter and process prominence calculation outputs, then relocate peaks to exact local maxima (corrects for tiling artifacts):

```bash
python code/organize_prominence_results.py \
    --folder_prominence_results prominence/ \
    --folder_organised_results data/output/OUTPUT_FOLDER/ \
    --bathymetry_file input_bathy.flt \
    --window_size 5
```

This step:
- Converts results.txt to geopackage format
- Filters invalid peaks (islands at elevation = 0, continent margin artifacts)
- Relocates peaks to local bathymetry maxima within a search window
- Saves final peak locations with relocation metadata

### Step 4: Multi-Level Contour Extraction

Extract contours at multiple prominence levels (100%, 90%, 75%, 50%, 25%):

```bash
python code/find_ctrs.py \
    --dem_file input.nc \
    --gpkg_file data/output/OUTPUT_FOLDER/bathymetry_features_peaks.gpkg \
    --window_size_min 1 \
    --window_size_max 25 \
    --enable_dateline_stitching \
    --batch_size 5000 \
    --num_workers 18
```

**Key Parameters:**
- `--window_size_min/max`: Adaptive window sizing (degrees) for contour search
- `--enable_dateline_stitching`: Enable for features crossing ±180° longitude (can lead to visualization issues on UTM maps)
- `--batch_size`: Number of features per batch (memory optimization)
- `--num_workers`: Parallel processing groups/threads

**Supported Formats:**
- NetCDF (`.nc`)
- GeoTIFF (`.tif`, `.tiff`)
- VRT (`.vrt`)

**⚠️ CRITICAL: Window Size Parameter**

The `--window_size_max` parameter is **critical for runtime performance**. The value of 25° was used for global datasets to ensure large plateau-like features are captured, but this results in significantly longer processing times. For regional/local analysis, decrease this value to 5-10° adequatly. The algorithm adaptively searches from `window_size_min` to `window_size_max` until all contours are found. Larger max values mean more iterations and exponentially more data loading per feature.

### Step 5: Nesting Analysis & Hierarchical Clustering

Detect feature nesting relationships and prepare for visualization:

```bash
python code/ctrs_cluster.py \
    --input_gpkg data/output/OUTPUT_FOLDER/contours.gpkg \
    --output_gpkg data/output/OUTPUT_FOLDER/bathymetry_features_contours.gpkg \
    --tile_buffer_degrees 45
```

This step:
- Identifies which features are nested within larger features
- Sorts contours for proper rendering in QGIS (largest → smallest)
- Adds `nested_on_feature_id` attribute for hierarchical relationships

## Example Workflows

Two complete workflow scripts are provided in the repository root. These scripts demonstrate the entire pipeline from bathymetry preprocessing to final outputs using GEBCO global datasets.

### GEBCO 2014 (30 arc-second resolution)

```bash
./run_GEBCO_2014_30s.sh
```

**Configuration:**
- Resolution: 30 arc-seconds (~900m at equator)
- Samples per tile: 120 points/degree
- Output: `data/output/GEBCO2014_30s_1deg120pts/`

### GEBCO 2025 (15 arc-second resolution)

```bash
./run_GEBCO_2025_15s.sh
```

**Configuration:**
- Resolution: 15 arc-seconds (~450m at equator)
- Samples per tile: 240 points/degree
- Output: `data/output/GEBCO2025_15s_1deg240pts/`

**Note:** Before running these scripts, ensure you have:
1. Activated the conda environment: `conda activate BathyFeaturesCtrs`
2. Compiled the Mountains binary (see [Installation](#installation))
3. Placed your input bathymetry data in the appropriate location (modify paths in scripts as needed)

## Output Files

Each workflow produces the following outputs in the specified folder:

### Final Outputs
- `bathymetry_features_peaks.gpkg` - Peak locations with prominence values
- `bathymetry_features_contours.gpkg` - Multi-level contours with nesting relationships


### Output Attributes

**Peaks GeoPackage:**
- `peak_id` - Unique feature identifier
- `longitude`, `latitude` - Peak location (WGS84, after relocation, 5 decimal precision)
- `elevation` - Peak elevation (m, negative for bathymetry, 2 decimal precision)
- `prominence` - Topographic prominence (m, 2 decimal precision)
- `key_saddle_latitude`, `key_saddle_longitude` - Key saddle location (5 decimal precision)
- `original_lon`, `original_lat` - Position before relocation (5 decimal precision)
- `original_elevation` - Elevation before relocation (2 decimal precision)
- `error_with_contours` - Boolean flag for contour extraction errors
- `error_no_contours` - Boolean flag for missing contours

**Contours GeoPackage:**
- `feature_id` - Links to parent peak
- `peak_longitude`, `peak_latitude` - Associated peak location (5 decimal precision)
- `peak_elevation` - Peak elevation (m, 2 decimal precision)
- `peak_prominence` - Peak prominence (m, 2 decimal precision)
- `contour_name` - Display label (e.g., "feature id 123, ctr 75%")
- `prominence_percentage` - Contour level (100%, 90%, 75%, 50%, 25%)
- `contour_elevation` - Contour elevation (m, 2 decimal precision)
- `nested_on_feature_id` - Comma-separated IDs of parent features (if nested)
- `area_sq_km` - Contour area (km², 2 decimal precision)
- `circularity_percent` - Shape circularity (0-100%, 2 decimal precision)
- `bbox_orientation_deg` - Minimum rotated rectangle orientation (0-180°, 2 decimal precision)
- `bbox_length_m` - Length of minimum rotated rectangle (m, 2 decimal precision)
- `bbox_width_m` - Width of minimum rotated rectangle (m, 2 decimal precision)
- `mean_slope_deg`, `min_slope_deg`, `max_slope_deg` - Slope statistics (2 decimal precision)
- `centroid_lon`, `centroid_lat` - Contour centroid (5 decimal precision)
- `window_size_used` - Search window size where contour was found (degrees, 2 decimal precision)
- `geometry` - Polygon geometry with vertices at 5 decimal precision

## Visualization in QGIS

1. Open the final `bathymetry_features_contours.gpkg` in QGIS
2. Contours are pre-sorted for proper rendering (largest at bottom, smallest on top)
3. Style by `prominence_percentage` to visualize feature hierarchy
4. Use `nested_on_feature_id` to identify feature relationships
5. Filter by `area_sq_km` or `circularity_percent` for feature classification

## Performance Optimization

### Memory Management
- Adjust `--batch_size` based on available RAM (larger = faster but more memory)
- The pipeline processes features in batches to limit memory usage
- Intermediate files are deleted after each step

### Parallel Processing
- Set `--num_workers` to match available CPU cores
- GDAL cache and thread settings are pre-configured in scripts
- For very large datasets (>100M points), consider spatial subsetting

### Adaptive Window Sizing
- Contour extraction uses adaptive windows (min → max)
- Only increases window size when contours are not found
- Tracks window size per contour level for diagnostics

## Troubleshooting

### Peak Extraction Errors
- Some peaks may fail contour extraction due to invalid nesting (e.g., 25% contour larger than 100%)
- These are flagged in `bathymetry_features_peaks.gpkg` with error attributes
- Common causes: Data pixel artifacts

### Memory Issues
- Reduce `--batch_size` parameter
- Reduce `--num_workers` to lower parallel memory usage
- Process data in geographic subsets

### Missing Contours
- Increase `--window_size_max` for large features
- Check that prominence values are sufficient (>300m recommended)
- Verify input data has sufficient resolution

## File Cleanup

Temporary files are automatically removed:
- `prominence/` - Raw prominence calculation outputs
- `tiles/` - Intermediate prominence tiling data


## Recent Changes


## Citation

If using this pipeline for research, please cite XX

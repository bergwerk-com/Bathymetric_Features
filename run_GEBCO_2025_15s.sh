#!/bin/bash

# Increase the argument limit for prompt. Can be an issue with length limitation running the data merge step in prominence.
ulimit -s 600000

# Input file paths
input_nc="../../Bathymetry_data/gebco_2025/GEBCO_2025.nc"
input_flt_full="../../Bathymetry_data/gebco_2025/GEBCO_2025_full.flt"
input_flt_bathy="../../Bathymetry_data/gebco_2025/GEBCO_2025_bathy_only.flt"

# Output folder results
folder_results="data/output/GEBCO2025_15s_1deg240pts"

# Bathymetry preprocessing
# Convert NetCDF directly to .flt format first
gdal_translate \
    -of EHdr \
    -ot Float32 \
    NETCDF:"$input_nc":elevation "$input_flt_full"
# Filter to bathymetry only (values <= 0) using gdal_calc on the .flt
gdal_calc.py \
    -A "$input_flt_full" \
    --outfile="$input_flt_bathy" \
    --calc="A*(A<=0)" \
    --NoDataValue=-32768 \
    --format=EHdr \
    --overwrite
# Build overviews for visualization
gdaladdo "$input_flt_bathy" 2 4 8 16 32


# Execute mountains script to find peaks and corresponding prominences
# sample_per_tiles argument set to 240 points for 1 degree, corresponding to the resolution of GEBCO 15s dataset. 
python code/mountains/scripts/run_prominence.py  \
      --binary_dir code/mountains/code/release \
      --threads 18  \
      --degrees_per_tile 1 \
      --samples_per_tile 240 \
      --skip_boundary \
      --min_prominence 300 \
      --bathymetry \
    "$input_flt_bathy"

# Sort, organize, filter and relocate peaks from mountains calculations
echo "Organize and relocate peaks"
output_prominence="prominence/"
python code/organize_prominence_results.py \
    --folder_prominence_results "$output_prominence" \
    --folder_organised_results "$folder_results" \
    --bathymetry_file "$input_flt_bathy" \
    --window_size 5

# Clean up before extracting the contours 
echo "Cleaning up mountains calculation files and folders"
[ -d "prominence" ] && rm -rf prominence/ && echo "Removed prominence/"
[ -d "tiles" ] && rm -rf tiles/ && echo "Removed tiles/"
[ -f "$folder_results/all_results.gpkg" ] && rm -f "$folder_results/all_results.gpkg" && echo "Removed temporary all_results.gpkg file"
[ -f "$folder_results/all_results_processed.gpkg" ] && rm -f "$folder_results/all_results_processed.gpkg" && echo "Removed temporary all_results_processed.gpkg file"

# Execute the base contours extraction for each topographic feature and save to geopackage
echo "Contours base exctraction"
python code/find_ctrs.py --dem_file="$input_nc" \
                      --gpkg_file="$folder_results/bathymetry_peaks.gpkg" \
                      --window_size_min=1 \
                      --window_size_max=25 \
                      --enable_dateline_stitching \
                      --batch_size=5000 \
                      --num_workers=18 

# Peak hierarchy and cluster visualization 
python code/cluster_ctrs.py --input_gpkg="$folder_results/contours.gpkg" \
                       --output_gpkg="$folder_results/bathymetry_features.gpkg" \
                       --tile_buffer_degrees=45

# Final clean up of temporary files
[ -f "$folder_results/contours.gpkg" ] && rm -f "$folder_results/contours.gpkg" && echo "Removed temporary contours.gpkg file"

echo "Done. All results are in the folder $folder_results."  
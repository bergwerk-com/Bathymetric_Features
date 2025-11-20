#!/usr/bin/env python3
"""
Script to post-process prominence analysis results by:
- Converting results.txt to geopackages
- Filtering invalid peaks adjacent to NoData cells (islands and continents)
- Relocating peaks to local bathymetry maxima
"""

import argparse
from pathlib import Path
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
from tqdm import tqdm


def relocate_peak_to_local_maximum(peak_lon, peak_lat, bathymetry_src, window_size=3):
    """
    Relocate a peak to the local maximum within a window around the original position.
    Also filters out peaks if NoData/NaN values are present in the relocation window.

    Parameters:
    -----------
    peak_lon : float
        Peak longitude
    peak_lat : float
        Peak latitude
    bathymetry_src : rasterio dataset
        Open bathymetry raster
    window_size : int
        Size of search window (e.g., 3 for 3x3 window)

    Returns:
    --------
    tuple: (new_lon, new_lat, elevation, is_valid)
        - new_lon, new_lat: relocated coordinates (or original if failed)
        - elevation: elevation at relocated position (or None if failed)
        - is_valid: False if NoData/NaN found in window, True otherwise
    """
    transform = bathymetry_src.transform
    nodata = bathymetry_src.nodata

    # Convert lon/lat to pixel coordinates
    col, row = ~transform * (peak_lon, peak_lat)
    col, row = int(round(col)), int(round(row))

    # Check if point is within raster bounds
    if not (0 <= row < bathymetry_src.height and 0 <= col < bathymetry_src.width):
        return peak_lon, peak_lat, None, False

    # Define window around the peak
    half_window = window_size // 2
    window_row_start = max(0, row - half_window)
    window_row_end = min(bathymetry_src.height, row + half_window + 1)
    window_col_start = max(0, col - half_window)
    window_col_end = min(bathymetry_src.width, col + half_window + 1)

    # Read the window
    window = rasterio.windows.Window(
        window_col_start, window_row_start,
        window_col_end - window_col_start,
        window_row_end - window_row_start
    )

    try:
        data = bathymetry_src.read(1, window=window)

        # Check for NoData/NaN values in the window
        has_nodata = False
        if nodata is not None:
            has_nodata = np.any(data == nodata)
        has_nan = np.any(np.isnan(data))

        # If NoData/NaN found in window, mark peak as invalid
        if has_nodata or has_nan:
            return peak_lon, peak_lat, None, False

        # Mask NoData values for finding maximum
        if nodata is not None:
            data_masked = np.ma.masked_equal(data, nodata)
        else:
            data_masked = np.ma.masked_invalid(data)

        # Find maximum value (for bathymetry, maximum means least negative = shallowest/highest)
        if data_masked.count() == 0:
            # No valid data in window
            return peak_lon, peak_lat, None, False

        # Find the position of maximum value within the window
        max_idx = np.unravel_index(data_masked.argmax(), data_masked.shape)
        max_elevation = float(data_masked.max())

        # Convert window-relative position to absolute pixel coordinates
        max_row = window_row_start + max_idx[0]
        max_col = window_col_start + max_idx[1]

        # Convert pixel coordinates to lon/lat (use pixel center)
        new_lon, new_lat = transform * (max_col + 0.5, max_row + 0.5)

        return new_lon, new_lat, max_elevation, True

    except Exception as e:
        print(f"Warning: Could not relocate peak at ({peak_lon}, {peak_lat}): {e}")
        return peak_lon, peak_lat, None, False


def round_gdf_precision(gdf):
    """
    Round coordinate and float values to reduce file size.
    - Longitude/latitude: 5 decimal places (~1.1m precision)
    - Other floats: 2 decimal places
    """
    # Longitude/latitude columns: 5 decimal places
    lon_lat_cols = ['latitude', 'longitude', 'key_saddle_latitude', 'key_saddle_longitude',
                    'original_lat', 'original_lon']

    for col in lon_lat_cols:
        if col in gdf.columns:
            gdf[col] = gdf[col].round(5)

    # Other float columns: 2 decimal places
    float_cols_2decimals = ['elevation', 'prominence', 'original_elevation']

    for col in float_cols_2decimals:
        if col in gdf.columns:
            gdf[col] = gdf[col].round(2)

    # Round geometry coordinates to 5 decimal places
    gdf.geometry = gdf.geometry.apply(lambda geom: gpd.points_from_xy([round(geom.x, 5)], [round(geom.y, 5)])[0])

    return gdf


def relocate_peaks(gdf, bathymetry_path, window_size=3):
    """
    Relocate all peaks in a GeoDataFrame to local maxima in the bathymetry.
    Also filters out peaks with NoData/NaN values in their relocation window.

    Parameters:
    -----------
    gdf : GeoDataFrame
        Input GeoDataFrame with peaks
    bathymetry_path : str
        Path to bathymetry raster (FLT, TIFF, VRT, or NetCDF)
    window_size : int
        Size of search window (default: 3 for 3x3)

    Returns:
    --------
    GeoDataFrame with relocated peaks (invalid peaks removed)
    """
    print(f"\nOpening bathymetry data: {bathymetry_path}")
    with rasterio.open(bathymetry_path) as src:
        print(f"Bathymetry grid: {src.width} x {src.height} pixels")
        print(f"Pixel size: {abs(src.transform[0]):.6f}° x {abs(src.transform[4]):.6f}°")
        print(f"NoData value: {src.nodata}")

        # Store original coordinates
        gdf['original_lon'] = gdf.geometry.x
        gdf['original_lat'] = gdf.geometry.y

        # Relocate each peak and track validity
        new_lons = []
        new_lats = []
        relocated_elevations = []
        relocation_distances = []
        valid_mask = []

        print(f"Relocating peaks using {window_size}x{window_size} window (filtering NoData/NaN)...")
        for idx, row in tqdm(gdf.iterrows(), total=len(gdf), desc="Relocating peaks"):
            peak_lon = row.geometry.x
            peak_lat = row.geometry.y

            new_lon, new_lat, elevation, is_valid = relocate_peak_to_local_maximum(
                peak_lon, peak_lat, src, window_size
            )

            new_lons.append(new_lon)
            new_lats.append(new_lat)
            relocated_elevations.append(elevation if elevation is not None else row.get('elevation', np.nan))
            valid_mask.append(is_valid)

            # Calculate relocation distance in degrees
            distance = np.sqrt((new_lon - peak_lon)**2 + (new_lat - peak_lat)**2)
            relocation_distances.append(distance)

        # Update geometry with relocated coordinates
        gdf.geometry = gpd.points_from_xy(new_lons, new_lats)

        # Update elevation if available
        if 'elevation' in gdf.columns:
            gdf['original_elevation'] = gdf['elevation']
            gdf['elevation'] = relocated_elevations

        # Filter out invalid peaks (those with NoData/NaN in window)
        initial_count = len(gdf)
        gdf = gdf[valid_mask].copy()
        removed_count = initial_count - len(gdf)
        print(f"Removed {removed_count} peaks with NoData/NaN in {window_size}x{window_size} relocation window")

        # Calculate statistics (only for valid peaks)
        valid_distances = [relocation_distances[i] for i in range(initial_count) if valid_mask[i]]
        moved_peaks = sum(d > 0 for d in valid_distances)
        avg_distance = np.mean([d for d in valid_distances if d > 0]) if moved_peaks > 0 else 0
        max_distance = max(valid_distances) if valid_distances else 0

        print(f"\nRelocation summary:")
        print(f"  Valid peaks after filtering: {len(gdf)}")
        print(f"  Peaks relocated: {moved_peaks} ({moved_peaks/len(gdf)*100:.1f}%)")
        print(f"  Peaks unchanged: {len(gdf) - moved_peaks}")
        print(f"  Average relocation distance: {avg_distance:.6f}°")
        print(f"  Maximum relocation distance: {max_distance:.6f}°")

    return gdf


def main():
    parser = argparse.ArgumentParser(description='Post-process prominence analysis results')
    parser.add_argument('--folder_prominence_results', required=True,
                       help='Path to folder containing prominence/results.txt')
    parser.add_argument('--folder_organised_results', required=True,
                       help='Path to output folder for organized results')
    parser.add_argument('--bathymetry_file', required=True,
                       help='Path to bathymetry file (FLT, TIFF, VRT) for NoData filtering and peak relocation')
    parser.add_argument('--window_size', type=int, default=3,
                       help='Search window size for peak relocation (default: 3 for 3x3 window)')

    args = parser.parse_args()

    # Set up paths
    prominence_dir = Path(args.folder_prominence_results)
    output_dir = Path(args.folder_organised_results)
    results_file = prominence_dir / 'results.txt'

    # Create output directory
    output_dir.mkdir(exist_ok=True)

    # Read results file
    print("Reading prominence results...")
    column_names = ['latitude', 'longitude', 'elevation', 'key_saddle_latitude', 'key_saddle_longitude', 'prominence']
    df = pd.read_csv(results_file, header=None, names=column_names, sep=',')
    print(f"Loaded {len(df)} peaks from results.txt")

    # Convert to GeoDataFrame
    # Note: The prominence binary works on a resampled grid (samples_per_tile resolution)
    # which is different from the input bathymetry grid. Peak coordinates come from
    # this finer internal grid, not the original bathymetry pixels.
    geometry = gpd.points_from_xy(df.longitude, df.latitude)
    gdf = gpd.GeoDataFrame(df, geometry=geometry, crs="EPSG:4326")

    # Save unprocessed results
    print(f"\nSaving unprocessed results to {output_dir / 'all_results.gpkg'}")
    gdf_rounded = round_gdf_precision(gdf.copy())
    gdf_rounded.to_file(output_dir / 'all_results.gpkg', driver="GPKG")

    # Add peak_id column (start at 1)
    gdf['peak_id'] = gdf.index + 1
    # Place as first column
    cols = [col for col in gdf.columns if col != 'peak_id']
    gdf = gdf[['peak_id'] + cols]

    # Remove islands at elevation=0
    print("\nFiltering invalid peaks...")
    initial_count = len(gdf)
    gdf = gdf[gdf['elevation'] < 0]
    print(f"Removed {initial_count - len(gdf)} peaks at elevation >= 0 (islands)")

    # Remove bad points at continent margin (high prominence values)
    initial_count = len(gdf)
    gdf = gdf[gdf['prominence'] < 9000]
    gdf = gdf[gdf['key_saddle_latitude'] != 0]
    gdf = gdf[gdf['key_saddle_longitude'] != 0]
    print(f"Removed {initial_count - len(gdf)} peaks with invalid prominence or saddle coordinates")

    # Save processed (but not relocated) results
    print(f"\nSaving filtered results to {output_dir / 'all_results_processed.gpkg'}")
    gdf_rounded = round_gdf_precision(gdf.copy())
    gdf_rounded.to_file(output_dir / 'all_results_processed.gpkg', driver="GPKG")

    # Relocate peaks to local bathymetry maxima (also filters NoData/NaN within relocation window)
    gdf = relocate_peaks(gdf, args.bathymetry_file, args.window_size)

    # Rename elevation to depth
    gdf = gdf.rename(columns={"elevation": "depth"})

    # Save final relocated peaks
    output_path = output_dir / 'bathymetry_peaks.gpkg'
    print(f"\nSaving relocated peaks to {output_path}")
    gdf_rounded = round_gdf_precision(gdf.copy())
    gdf_rounded.to_file(output_path, driver="GPKG")
    print(f"Saved {len(gdf)} relocated peaks")

    print("\nDone!")


if __name__ == "__main__":
    main()

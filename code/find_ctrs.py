#!/usr/bin/env python3
"""
Basic contour extraction: Multi-level contours with slope analysis
Extracts contours without nesting analysis (use find_ctrs_cluster.py for nesting)

- Adaptive crop window sizing to find the contours
- Optional dateline stitching for longitude at ±180° (default: disabled)
"""

import argparse
from tqdm import tqdm
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import rasterio
from rasterio.windows import from_bounds
import contourpy
from shapely.geometry import Polygon, Point
from shapely.ops import transform
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import pyproj
import gc
from functools import lru_cache
import threading
import psutil

# Thread-local storage for persistent file handles
_thread_local = threading.local()


def get_rasterio_handle(file_path):
    """
    Get a persistent rasterio file handle for the current process/thread.
    Reuses the same handle across multiple calls to reduce I/O overhead.
    """
    if not hasattr(_thread_local, 'raster_handle'):
        _thread_local.raster_handle = None
        _thread_local.raster_path = None

    # Open new handle if needed
    if _thread_local.raster_path != file_path:
        # Close old handle if exists
        if _thread_local.raster_handle is not None:
            try:
                _thread_local.raster_handle.close()
            except:
                pass

        # Open new handle
        _thread_local.raster_handle = rasterio.Env(
            GDAL_CACHEMAX=1024,  # 1GB cache per worker
            GDAL_NUM_THREADS='ALL_CPUS',
            GDAL_DISABLE_READDIR_ON_OPEN='EMPTY_DIR'
        )
        _thread_local.raster_handle.__enter__()
        _thread_local.raster_src = rasterio.open(file_path)
        _thread_local.raster_path = file_path

    return _thread_local.raster_src


def get_xarray_handle(file_path):
    """
    Get a persistent xarray dataset handle for the current process/thread.
    Reuses the same handle across multiple calls.
    """
    if not hasattr(_thread_local, 'xr_handle'):
        _thread_local.xr_handle = None
        _thread_local.xr_path = None

    if _thread_local.xr_path != file_path:
        # Close old handle if exists
        if _thread_local.xr_handle is not None:
            try:
                _thread_local.xr_handle.close()
            except:
                pass

        # Open new handle with optimized chunking
        _thread_local.xr_handle = xr.open_dataset(
            file_path,
            chunks={'lat': 1000, 'lon': 1000},
            cache=True
        )
        _thread_local.xr_path = file_path

    return _thread_local.xr_handle


def calculate_slope_along_polygon_path(polygon, lon_grid, lat_grid, dem_data):
    """
    Calculate mean, min, and max slope along the polygon perimeter
    Memory-optimized version: compute gradients once and clean up intermediate arrays
    """
    try:
        # Get polygon boundary coordinates
        boundary_coords = np.array(polygon.exterior.coords)

        if len(boundary_coords) < 3:
            return 0.0, 0.0, 0.0

        # Get the bounds and resolution
        lon_coords = lon_grid[0, :]
        lat_coords = lat_grid[:, 0]
        lon_res = np.abs(lon_coords[1] - lon_coords[0])
        lat_res = np.abs(lat_coords[1] - lat_coords[0])

        # Convert to geographic scaling
        center_lat = np.mean(lat_coords)
        meters_per_deg_lat = 111000
        meters_per_deg_lon = 111000 * np.cos(np.radians(center_lat))

        # Calculate gradients and convert to slope in one step to reduce memory
        grad_y, grad_x = np.gradient(dem_data)

        # Compute slope in-place to avoid extra arrays
        slope_x = grad_x / (lon_res * meters_per_deg_lon)
        del grad_x  # Free memory immediately
        slope_y = grad_y / (lat_res * meters_per_deg_lat)
        del grad_y  # Free memory immediately

        # Calculate slope magnitude in-place
        np.square(slope_x, out=slope_x)
        np.square(slope_y, out=slope_y)
        slope_magnitude = slope_x + slope_y  # Reuse slope_x array
        del slope_x, slope_y  # Free memory
        np.sqrt(slope_magnitude, out=slope_magnitude)

        # Convert to degrees
        slope_degrees = np.degrees(np.arctan(slope_magnitude))
        del slope_magnitude  # Free memory

        # Interpolate slope values along polygon path
        slopes_along_path = []

        for coord in boundary_coords:
            lon, lat = coord[0], coord[1]

            # Convert to grid indices
            col = np.interp(lon, lon_coords, np.arange(len(lon_coords)))
            row = np.interp(lat, lat_coords, np.arange(len(lat_coords)))

            # Bilinear interpolation
            row = max(0, min(row, slope_degrees.shape[0] - 1))
            col = max(0, min(col, slope_degrees.shape[1] - 1))

            row_i = int(row)
            col_i = int(col)
            row_f = row - row_i
            col_f = col - col_i

            # Ensure indices are within bounds
            row_i = max(0, min(row_i, slope_degrees.shape[0] - 1))
            col_i = max(0, min(col_i, slope_degrees.shape[1] - 1))
            row_i1 = min(row_i + 1, slope_degrees.shape[0] - 1)
            col_i1 = min(col_i + 1, slope_degrees.shape[1] - 1)

            # Bilinear interpolation
            if row_i == row_i1 and col_i == col_i1:
                slope_value = slope_degrees[row_i, col_i]
            else:
                v00 = slope_degrees[row_i, col_i]
                v01 = slope_degrees[row_i, col_i1]
                v10 = slope_degrees[row_i1, col_i]
                v11 = slope_degrees[row_i1, col_i1]

                slope_value = (v00 * (1 - row_f) * (1 - col_f) +
                             v01 * (1 - row_f) * col_f +
                             v10 * row_f * (1 - col_f) +
                             v11 * row_f * col_f)

            if not np.isnan(slope_value) and slope_value >= 0 and slope_value < 90:
                slopes_along_path.append(slope_value)

        # Calculate mean, min, and max slope along path
        if slopes_along_path:
            mean_slope = np.mean(slopes_along_path)
            min_slope = np.min(slopes_along_path)
            max_slope = np.max(slopes_along_path)
        else:
            mean_slope = min_slope = max_slope = 0.0

        return (
            mean_slope if not np.isnan(mean_slope) else 0.0,
            min_slope if not np.isnan(min_slope) else 0.0,
            max_slope if not np.isnan(max_slope) else 0.0
        )

    except Exception:
        return 0.0, 0.0, 0.0


def calculate_polygon_orientation_bbox(projected_polygon):
    """
    Calculate polygon orientation and bounding box dimensions using minimum bounding rectangle method.
    Uses the same utm projected polygon as used for area/circularity calculations.
    Returns:
        bbox_orientation_deg: angle in degrees using geographic convention
            0° = North, 90° = East, 180° = South, 270° = West
            Range: 0-180° (since orientation is bidirectional)
        bbox_length_m: Length of the longer side in meters
        bbox_width_m: Width of the shorter side in meters
    """
    try:
        if projected_polygon is None:
            return None, None, None

        # Get minimum rotated rectangle (from Shapely)
        min_rect = projected_polygon.minimum_rotated_rectangle

        # Get rectangle coordinates
        rect_coords = list(min_rect.exterior.coords[:-1])  # Remove duplicate last point

        if len(rect_coords) < 4:
            return None, None, None

        # Calculate edge vectors (in projected coordinates - meters)
        edge1 = np.array(rect_coords[1]) - np.array(rect_coords[0])
        edge2 = np.array(rect_coords[2]) - np.array(rect_coords[1])

        # Calculate edge lengths
        edge1_length = np.linalg.norm(edge1)
        edge2_length = np.linalg.norm(edge2)

        # Length is the longer edge, width is the shorter edge
        bbox_length_m = max(edge1_length, edge2_length)
        bbox_width_m = min(edge1_length, edge2_length)

        # Choose the longer edge as the primary orientation
        if edge1_length >= edge2_length:
            primary_edge = edge1
        else:
            primary_edge = edge2

        # Calculate angle in mathematical convention (for projected coordinates)
        # In UTM/projected coordinates: +X is East, +Y is North
        angle_rad = np.arctan2(primary_edge[1], primary_edge[0])
        angle_deg_math = np.degrees(angle_rad)

        # Convert to geographic convention (0° = North)
        # Mathematical: 0°=East, 90°=North
        # Geographic: 0°=North, 90°=East
        angle_deg_geo = (90 - angle_deg_math) % 360

        # Normalize to 0-180° range (since orientation is bidirectional)
        if angle_deg_geo > 180:
            angle_deg_geo -= 180

        return angle_deg_geo, bbox_length_m, bbox_width_m

    except Exception:
        return None, None, None


# Cache transformer objects to avoid repeated creation
@lru_cache(maxsize=128)
def get_transformer(lat, lon):
    """
    Get a cached pyproj Transformer for the given location.
    LRU cache avoids creating duplicate transformers for nearby seamounts.
    """
    # Determine appropriate projection
    if abs(lat) > 80:
        # For polar regions, use Lambert Azimuthal Equal Area
        if lat > 0:
            proj_string = '+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
        else:
            proj_string = '+proj=laea +lat_0=-90 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs'
    else:
        # Use UTM for most locations
        utm_zone = int((lon + 180) / 6) + 1
        if lat >= 0:
            proj_string = f'+proj=utm +zone={utm_zone} +datum=WGS84 +units=m +no_defs'
        else:
            proj_string = f'+proj=utm +zone={utm_zone} +south +datum=WGS84 +units=m +no_defs'

    return pyproj.Transformer.from_crs('EPSG:4326', proj_string, always_xy=True)


def calculate_area_circularity_and_centroid(polygon, lon, lat):
    """
    Calculate area, circularity and centroid using proper projection
    Returns area, circularity, centroid coordinates, and projected polygon
    """
    try:
        # Get cached transformer (rounded to 0.1 degree for better cache hits)
        transformer = get_transformer(round(lat, 1), round(lon, 1))

        # Project polygon to meters
        projected_polygon = transform(transformer.transform, polygon)

        # Calculate area in square meters
        area_sq_m = projected_polygon.area

        # Calculate circularity using projected polygon
        if projected_polygon.length == 0:
            circularity = 0.0
        else:
            circularity = (4 * np.pi * projected_polygon.area) / (projected_polygon.length ** 2)
            circularity = min(circularity * 100, 100.0)  # Convert to percentage

        # Centroid in original coordinates
        centroid = polygon.centroid
        centroid_lon = centroid.x
        centroid_lat = centroid.y

        return area_sq_m, circularity, centroid_lon, centroid_lat, projected_polygon

    except Exception:
        return None, None, None, None, None


def normalize_longitude(lon):
    """
    Normalize longitude to [-180, 180] range
    """
    return ((lon + 180) % 360) - 180


def check_longitude_wraparound(lon_min, lon_max):
    """
    Check if a longitude window crosses the 180°/-180° boundary
    Returns True if wraparound occurs
    """
    # Normalize both values first
    lon_min_norm = normalize_longitude(lon_min)
    lon_max_norm = normalize_longitude(lon_max)

    # Check if we cross the dateline
    return lon_min_norm > lon_max_norm


def load_data_with_wraparound(file_dem, lon_min, lon_max, lat_min, lat_max, enable_stitching=True, save_crop=False, crop_output_path=None):
    """
    Load DEM data handling longitude wraparound at 180°/-180° dateline boundary

    NOTE: Only stitches in LONGITUDE (E-W direction at ±180° dateline).
    Never stitches in LATITUDE (N-S direction) - latitude is always clipped to [-90, 90].

    Args:
        file_dem: Path to DEM file
        lon_min, lon_max: Longitude bounds (can cross ±180°)
        lat_min, lat_max: Latitude bounds (will be clipped to [-90, 90])
        enable_stitching: If True, stitch data across dateline; if False, clip to [-180, 180]
        save_crop: If True, save the cropped window as a TIF file
        crop_output_path: Path to save the cropped TIF
    """
    file_ext = os.path.splitext(file_dem.lower())[1]

    # LATITUDE: Always clip to valid range [-90, 90], never stitch
    lat_min = max(-90, min(90, lat_min))
    lat_max = max(-90, min(90, lat_max))

    # LONGITUDE: Normalize and check for dateline crossing (E-W stitching only)
    lon_min_norm = normalize_longitude(lon_min)
    lon_max_norm = normalize_longitude(lon_max)

    # Check if we cross the dateline (lon_min_norm > lon_max_norm after normalization)
    crosses_dateline = lon_min_norm > lon_max_norm and enable_stitching

    if crosses_dateline:
        # Split into two windows: western and eastern parts
        # Western part: [lon_min_norm, 180]
        # Eastern part: [-180, lon_max_norm]

        print(f"Dateline crossing detected: [{lon_min:.2f}, {lon_max:.2f}] -> stitching data")

        # Load western part
        west_data, west_lon_grid, west_lat_grid = _load_single_window_optimized(
            file_dem, lon_min_norm, 180, lat_min, lat_max, file_ext
        )

        # Load eastern part
        east_data, east_lon_grid, east_lat_grid = _load_single_window_optimized(
            file_dem, -180, lon_max_norm, lat_min, lat_max, file_ext
        )

        # Determine shift direction based on original peak position
        # Peak near +180°: shift eastern data to positive side (+360)
        # Peak near -180°: shift western data to negative side (-360)
        peak_lon_orig = normalize_longitude(lon_min + (lon_max - lon_min) / 2)  # approximate peak position

        if peak_lon_orig > 0:  # Peak on positive side, near +180°
            # Shift eastern data to positive side
            if west_data is not None and east_data is not None:
                east_lon_grid_shifted = east_lon_grid + 360
                dem_data = np.concatenate([west_data, east_data], axis=1)
                lon_grid = np.concatenate([west_lon_grid, east_lon_grid_shifted], axis=1)
                lat_grid = np.concatenate([west_lat_grid, east_lat_grid], axis=1)
            elif west_data is not None:
                dem_data, lon_grid, lat_grid = west_data, west_lon_grid, west_lat_grid
            elif east_data is not None:
                dem_data, lat_grid = east_data, east_lat_grid
                lon_grid = east_lon_grid + 360
        else:  # Peak on negative side, near -180°
            # Shift western data to negative side
            if west_data is not None and east_data is not None:
                west_lon_grid_shifted = west_lon_grid - 360
                dem_data = np.concatenate([west_data, east_data], axis=1)
                lon_grid = np.concatenate([west_lon_grid_shifted, east_lon_grid], axis=1)
                lat_grid = np.concatenate([west_lat_grid, east_lat_grid], axis=1)
            elif west_data is not None:
                dem_data, lat_grid = west_data, west_lat_grid
                lon_grid = west_lon_grid - 360
            elif east_data is not None:
                dem_data, lon_grid, lat_grid = east_data, east_lon_grid, east_lat_grid

        if dem_data is None:
            raise ValueError("No valid data found in dateline crossing region")

    else:
        # Normal case: no dateline crossing (or stitching disabled)
        if not enable_stitching and lon_min_norm > lon_max_norm:
            # Dateline crossing but stitching disabled - clip to [-180, 180]
            print(f"Dateline crossing detected but stitching disabled: clipping [{lon_min:.2f}, {lon_max:.2f}] to [-180, 180]")
            # Just load the western part (arbitrary choice - could also use eastern)
            lon_min_norm = max(-180, lon_min_norm)
            lon_max_norm = 180

        dem_data, lon_grid, lat_grid = _load_single_window_optimized(
            file_dem, lon_min_norm, lon_max_norm, lat_min, lat_max, file_ext,
            save_crop=save_crop, crop_output_path=crop_output_path
        )

    return dem_data, lon_grid, lat_grid


def _load_single_window_optimized(file_dem, lon_min, lon_max, lat_min, lat_max, file_ext, save_crop=False, crop_output_path=None):
    """
    Load a single rectangular window of DEM data
    Open/close per request - OS cache + GDAL cache handle the optimization

    Args:
        save_crop: If True, save the cropped window as a TIF file
        crop_output_path: Path to save the cropped TIF (required if save_crop=True)
    """
    try:
        if file_ext == '.nc':
            # NetCDF format
            ds = xr.open_dataset(file_dem, chunks={'lat': 1000, 'lon': 1000})
            ds['elevation'] = ds['elevation'].where(ds['elevation'] < 0, 0).astype(np.float32)

            cropped = ds.sel(lon=slice(lon_min, lon_max), lat=slice(lat_min, lat_max)).load()
            dem_data = cropped.elevation.values.copy()  # Make explicit copy

            lon_ = cropped.lon.values.copy()
            lat_ = cropped.lat.values.copy()
            lon_grid, lat_grid = np.meshgrid(lon_, lat_)

            # Explicitly clean up
            cropped.close()
            ds.close()
            del cropped, ds

            return dem_data, lon_grid, lat_grid

        elif file_ext in ['.tif', '.tiff', '.vrt']:
            # Raster format
            with rasterio.open(file_dem) as src:
                window = from_bounds(lon_min, lat_min, lon_max, lat_max, src.transform)
                dem_data_raw = src.read(1, window=window, masked=True)

                # Convert masked array to regular numpy array immediately
                if hasattr(dem_data_raw, 'mask'):
                    dem_data = np.where(dem_data_raw.mask, np.nan, dem_data_raw.data).astype(np.float32)
                    del dem_data_raw  # Release masked array
                else:
                    dem_data = dem_data_raw.astype(np.float32)

                dem_data = np.where(dem_data < 0, dem_data, 0)

                # FIX: Calculate the ACTUAL window that was read (accounting for clipping at file boundaries)
                # When the requested window extends beyond file bounds, rasterio clips it automatically,
                # but we need to use the actual clipped window dimensions for coordinate generation
                from rasterio.windows import Window
                height, width = dem_data.shape

                # Calculate actual window bounds (clipped to file extent)
                actual_row_off = max(0, int(window.row_off))
                actual_col_off = max(0, int(window.col_off))
                actual_row_end = min(src.height, int(window.row_off + window.height))
                actual_col_end = min(src.width, int(window.col_off + window.width))
                actual_height = actual_row_end - actual_row_off
                actual_width = actual_col_end - actual_col_off

                # Create window object for the actual data that was read
                actual_window = Window(actual_col_off, actual_row_off, actual_width, actual_height)
                window_transform = src.window_transform(actual_window)

                # Create coordinate arrays for pixel centers
                lon_coords = []
                lat_coords = []

                for col in range(width):
                    lon, _ = window_transform * (col + 0.5, 0.5)  # +0.5 for pixel center
                    lon_coords.append(lon)

                for row in range(height):
                    _, lat = window_transform * (0.5, row + 0.5)  # +0.5 for pixel center
                    lat_coords.append(lat)

                # Create meshgrids
                lon_grid, lat_grid = np.meshgrid(lon_coords, lat_coords)

                # Save cropped window if requested
                if save_crop and crop_output_path:
                    # Use the ACTUAL window transform (not the requested bounds)
                    # This ensures the saved TIF has correct georeferencing
                    with rasterio.open(
                        crop_output_path,
                        'w',
                        driver='GTiff',
                        height=height,
                        width=width,
                        count=1,
                        dtype=dem_data.dtype,
                        crs='EPSG:4326',
                        transform=window_transform,  # Use actual window transform
                        nodata=-9999
                    ) as dst:
                        dst.write(dem_data, 1)

            return dem_data, lon_grid, lat_grid
        else:
            raise ValueError(f"Unsupported format: {file_ext}")

    except Exception as e:
        print(f"Warning: Failed to load window [{lon_min:.2f}, {lon_max:.2f}, {lat_min:.2f}, {lat_max:.2f}]: {e}")
        return None, None, None


def haversine_distance(lon1, lat1, lon2, lat2):
    """
    Calculate distance between two points in meters
    """
    # Convert to radians
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])

    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(np.sqrt(a))

    # Radius of Earth in meters
    R = 6371000
    return R * c


def _process_single_contour(cont_gen, target_depth, percentage, lon_grid, lat_grid, dem_data,
                           target_lon, lat, lon):
    """
    Helper function to process a single contour level.
    Returns a result dict if successful, None otherwise.
    """
    try:
        lines = cont_gen.lines(target_depth)

        if len(lines) == 0:
            return None

        # Find the contour that contains the peak
        selected_line = None
        for line in lines:
            if len(line) < 3:
                continue

            # Check if naturally closed
            is_naturally_closed = np.allclose(line[0], line[-1], atol=1e-10)
            if not is_naturally_closed:
                continue

            # Test if peak is inside (use adjusted longitude for wraparound)
            test_polygon = Polygon(line)
            if test_polygon.contains(Point(target_lon, lat)):
                selected_line = line
                break  # Found the contour containing the peak

        # If no contour found, skip this level
        if selected_line is None:
            return None

        # Process the selected contour
        if selected_line is not None:
            # Create polygon directly without shifting coordinates
            polygon = Polygon(selected_line)

            # Fix topology if needed
            if not polygon.is_valid:
                try:
                    polygon = polygon.buffer(0)
                    if polygon.geom_type == 'MultiPolygon':
                        # Take the largest polygon
                        polygon = max(polygon.geoms, key=lambda x: x.area)
                except:
                    return None

            # Calculate metrics
            area_sq_m, circularity, centroid_lon, centroid_lat, projected_polygon = calculate_area_circularity_and_centroid(polygon, lon, lat)

            if circularity is not None:
                # Calculate slope metrics (use filled data to avoid NaN issues)
                mean_slope, min_slope, max_slope = calculate_slope_along_polygon_path(
                    polygon, lon_grid, lat_grid, dem_data
                )

                # Calculate orientation and bbox dimensions using projected polygon
                bbox_orientation_deg, bbox_length_m, bbox_width_m = calculate_polygon_orientation_bbox(projected_polygon)

                return {
                    'polygon': polygon,
                    'depth': target_depth,
                    'prominence_percentage': percentage,
                    'area_sq_m': area_sq_m,
                    'circularity_percent': circularity,
                    'bbox_orientation_deg': bbox_orientation_deg,
                    'bbox_length_m': bbox_length_m,
                    'bbox_width_m': bbox_width_m,
                    'mean_slope_deg': mean_slope,
                    'min_slope_deg': min_slope,
                    'max_slope_deg': max_slope,
                    'centroid_lon': centroid_lon,
                    'centroid_lat': centroid_lat
                }

    except Exception:
        pass

    return None


def extract_multiple_contours(lon_grid, lat_grid, dem_data, peak_depth, prominence, lon, lat, contour_percentages):
    """
    Extract contours at multiple levels using contourpy
    Process from top (25%) to base (100%) for efficiency
    Generates all contour lines at once using a single contour generator, then processes each level if not found.
    """
    results = []

    # Determine target longitude based on data coordinate system
    # If data spans dateline, we may need to adjust peak longitude
    lon_grid_min = np.nanmin(lon_grid)
    lon_grid_max = np.nanmax(lon_grid)

    # If data spans > 180 degrees, it means we have stitched data across dateline
    if (lon_grid_max - lon_grid_min) > 180:
        # Adjust peak longitude to match stitched coordinate system
        if lon_grid_max > 180:  # Data shifted to positive side (peak near +180°)
            target_lon = lon if lon >= 0 else lon + 360
        elif lon_grid_min < -180:  # Data shifted to negative side (peak near -180°)
            target_lon = lon if lon < 0 else lon - 360
        else:
            target_lon = normalize_longitude(lon)
    else:
        # Normal case: use normalized longitude
        target_lon = normalize_longitude(lon)

    # Sort contour percentages in ascending order (25%, 50%, 75%, ..., 100%)
    # Process from top (25% near peak) to base (100%)
    sorted_percentages = sorted(contour_percentages)

    # Calculate target elevations for valid percentages only
    valid_targets = []
    dem_min = np.nanmin(dem_data)
    dem_max = np.nanmax(dem_data)

    for percentage in sorted_percentages:
        target_depth = peak_depth - (percentage / 100.0 * prominence)
        if dem_min <= target_depth <= dem_max:
            valid_targets.append((target_depth, percentage))

    if len(valid_targets) == 0:
        return results

    # Keep NaN as NaN - do not fill, Contourpy will not generate contours in NaN regions
    dem_data_filled = dem_data

    try:
        # Create contour generator ONCE and reuse it for all levels
        cont_gen = contourpy.contour_generator(lon_grid, lat_grid, dem_data_filled)

        # Process each target elevation in ascending order (top to base)
        for target_depth, percentage in valid_targets:
            result = _process_single_contour(
                cont_gen, target_depth, percentage,
                lon_grid, lat_grid, dem_data,
                target_lon, lat, lon
            )
            if result is not None:
                results.append(result)

    except Exception:
        pass

    # Sort results by prominence_percentage in descending order to ensure proper ordering
    results.sort(key=lambda x: x['prominence_percentage'], reverse=True)

    return results


def check_memory_usage(threshold_percent=85):
    """
    Check current memory usage and return True if below threshold
    threshold_percent: Maximum allowed memory usage percentage (default 85%)
    """
    try:
        memory_info = psutil.virtual_memory()
        return memory_info.percent < threshold_percent
    except:
        return True  # If check fails, assume it's okay to continue


class MemoryLimitExceeded(Exception):
    """Custom exception to signal memory limit was reached"""
    pass


def crop_and_process_peak(peak_data):
    """Worker function to process a single peak with adaptive window sizing"""
    idx, lon, lat, depth, prominence, file_dem, window_size_min, window_size_max, contour_percentages, enable_stitching, output_dir = peak_data

    try:
        # Check memory before processing - HARD STOP if at limit
        memory_info = psutil.virtual_memory()
        if memory_info.percent >= 90:
            # Raise exception to stop the entire run
            raise MemoryLimitExceeded(f"Memory usage at {memory_info.percent:.1f}% (limit: 90%)")
        # Adaptive window sizing with incremental contour search
        # Start with minimum window size and increase only when needed
        # Only search for missing contours in larger windows
        # Track the window size where each contour was found
        current_window_size = window_size_min
        contour_results = []
        contour_window_sizes = {}  # Maps prominence_percentage -> first window size where ctr found
        contour_areas = {}  # Maps prominence_percentage -> area in sq_m
        all_contours_found = False
        error_with_contours = False  # Flag to mark peaks with invalid contour nesting

        while current_window_size <= window_size_max and not all_contours_found and not error_with_contours:
            try:
                # Periodic memory check during window expansion - HARD STOP if at limit
                memory_info = psutil.virtual_memory()
                if memory_info.percent >= 90:
                    # Raise exception to stop the entire run
                    raise MemoryLimitExceeded(f"Memory usage at {memory_info.percent:.1f}% (limit: 90%)")

                # Calculate window for current size
                lat_rad = np.radians(lat)
                lon_window = current_window_size
                lat_window = current_window_size * np.cos(lat_rad)

                lon_min, lon_max = lon - lon_window / 2, lon + lon_window / 2
                lat_min, lat_max = lat - lat_window / 2, lat + lat_window / 2

                # Load data with optional dateline stitching (longitude only, never latitude)
                dem_data, lon_grid, lat_grid = load_data_with_wraparound(
                    file_dem, lon_min, lon_max, lat_min, lat_max, enable_stitching,
                    save_crop=False, crop_output_path=None
                )

                # Skip if no valid data
                if dem_data is None or np.isnan(dem_data).all():
                    # Clean up before continuing
                    del dem_data, lon_grid, lat_grid
                    gc.collect()
                    current_window_size += 1.0
                    continue

                # Determine which contour percentages are still missing
                if len(contour_results) > 0:
                    found_percentages = set(result['prominence_percentage'] for result in contour_results)
                    missing_percentages = [pct for pct in contour_percentages if pct not in found_percentages]
                else:
                    missing_percentages = contour_percentages

                # Only extract missing contours (not all contours)
                if missing_percentages:
                    new_contours = extract_multiple_contours(
                        lon_grid, lat_grid, dem_data, depth, prominence, lon, lat, missing_percentages
                    )

                    # Add newly found contours to results and record the window size and area
                    for contour in new_contours:
                        pct = contour['prominence_percentage']
                        area_sq_m = contour['area_sq_m']

                        # Only record window size if this contour wasn't found before
                        if pct not in contour_window_sizes:
                            contour_window_sizes[pct] = current_window_size
                            contour_areas[pct] = area_sq_m

                            # VALIDATION: Check if a larger contour (higher %) exists with smaller area
                            # This violates nesting rules and indicates a problem
                            for existing_pct, existing_area in contour_areas.items():
                                # If we found a lower contour level (smaller %) but it's LARGER in area
                                # than a higher contour level (larger %), this is invalid
                                if pct < existing_pct and area_sq_m > existing_area:
                                    error_with_contours = True
                                    break
                                # Or if we found a higher contour level (larger %) but it's SMALLER in area
                                # than a lower contour level (smaller %), this is also invalid
                                elif pct > existing_pct and area_sq_m < existing_area:
                                    error_with_contours = True
                                    break

                        if error_with_contours:
                            break

                    # Only extend results if no error detected
                    if not error_with_contours:
                        contour_results.extend(new_contours)

                # Check if all target contour levels were found
                if len(contour_results) > 0:
                    found_percentages = set(result['prominence_percentage'] for result in contour_results)
                    target_percentages = set(contour_percentages)
                    all_contours_found = found_percentages == target_percentages

                    # VALIDATION: If we have some contours but 25% is missing AND we've reached max window
                    # this is an error (25% should always exist if deeper contours exist)
                    # But if we haven't reached max window yet, continue searching
                    if (25 in target_percentages and 25 not in found_percentages and
                        len(found_percentages) > 0 and current_window_size >= window_size_max):
                        error_with_contours = True
                else:
                    all_contours_found = False

                # If not all contours found, try larger window
                if not all_contours_found and not error_with_contours:
                    # Clean up data before next iteration
                    del dem_data, lon_grid, lat_grid
                    gc.collect()
                    current_window_size += 1.0
                else:
                    # Clean up data when loop is done
                    del dem_data, lon_grid, lat_grid
                    gc.collect()

            except Exception:
                # Clean up on error
                try:
                    del dem_data, lon_grid, lat_grid
                    gc.collect()
                except:
                    pass
                # If there's an error with current window size, try larger
                current_window_size += 1.0
                continue

        # If error detected, discard all contours for this peak
        if error_with_contours:
            # Return empty results with error flag
            return idx, [], True

        # Format results with peak_id attribute and per-contour window sizes
        formatted_results = []
        for result in contour_results:
            pct = result['prominence_percentage']
            formatted_results.append({
                'peak_id': idx,
                'polygon': result['polygon'],
                'depth': result['depth'],
                'prominence_percentage': result['prominence_percentage'],
                'area_sq_km': result['area_sq_m']/1e6, # conversion 1e6 m2 == 1 km2
                'circularity_percent': result['circularity_percent'],
                'centroid_lon': result['centroid_lon'],
                'centroid_lat': result['centroid_lat'],
                'mean_slope_deg': result['mean_slope_deg'],
                'min_slope_deg': result['min_slope_deg'],
                'max_slope_deg': result['max_slope_deg'],
                'bbox_orientation_deg': result['bbox_orientation_deg'],
                'bbox_length_m': result['bbox_length_m'],
                'bbox_width_m': result['bbox_width_m'],
                'window_size_used': contour_window_sizes.get(pct, None),  # First window size where this contour was found
            })

        return idx, formatted_results, False

    except Exception:
        return idx, [], False
    finally:
        # Force garbage collection after each peak
        gc.collect()


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Extract seamount contours with slope analysis')

    parser.add_argument('--dem_file', type=str, required=True,
                       help='Path to the DEM/bathymetry file (.nc, .tif, .vrt)')
    parser.add_argument('--gpkg_file', type=str, required=True,
                       help='Path to the seamounts GeoPackage file')
    parser.add_argument('--batch_size', type=int, default=2000,
                       help='Number of seamounts per batch (default: 2000)')
    parser.add_argument('--num_workers', type=int, default=18,
                       help='Number of parallel workers (default: 18)')
    parser.add_argument('--window_size_min', type=float, default=1.0,
                       help='Minimum window size in degrees (default: 1.0)')
    parser.add_argument('--window_size_max', type=float, default=15.0,
                       help='Maximum window size in degrees (default: 15.0)')
    parser.add_argument('--limit_seamounts', type=int, default=None,
                       help='Limit number of seamounts for testing')
    parser.add_argument('--output_name', type=str, default='contours',
                       help='Base name for output files (default: contours)')
    parser.add_argument('--enable_dateline_stitching', action='store_true', default=False,
                       help='Enable stitching of data across the ±180° dateline (longitude only, never latitude)')

    return parser.parse_args()


def main():
    args = parse_arguments()

    # Calculate optimal GDAL cache size dynamically
    # Strategy: Reserve memory for system, leave buffer for peak processing
    num_workers = args.num_workers
    total_memory_gb = psutil.virtual_memory().total / (1024**3)

    # Reserve 30% for system + leave 30% buffer for peak data processing
    # Remaining 40% can be used for GDAL cache across all workers
    available_for_gdal = total_memory_gb * 0.4
    cache_per_worker_mb = int((available_for_gdal * 1024) / num_workers)

    # Clamp between reasonable bounds (min 128MB, max 2GB per worker)
    cache_per_worker_mb = max(128, min(cache_per_worker_mb, 2048))

    print(f"Memory configuration:")
    print(f"  Total system memory: {total_memory_gb:.1f} GB")
    print(f"  Number of workers: {num_workers}")
    print(f"  GDAL cache per worker: {cache_per_worker_mb} MB")
    print(f"  Total GDAL cache: {(cache_per_worker_mb * num_workers) / 1024:.1f} GB")

    # Optimized GDAL environment variables
    os.environ['GDAL_CACHEMAX'] = str(cache_per_worker_mb)
    os.environ['GDAL_MAX_DATASET_POOL_SIZE'] = '100'  # Reduced to prevent too many open files
    os.environ['GDAL_NUM_THREADS'] = '1'  # Limit threads per worker to avoid oversubscription
    os.environ['GDAL_DISABLE_READDIR_ON_OPEN'] = 'EMPTY_DIR'
    os.environ['CPL_VSIL_CURL_ALLOWED_EXTENSIONS'] = '.tif,.tiff,.vrt'

    # Configuration
    name_file = args.output_name
    file_dem = args.dem_file
    gpkg_file = args.gpkg_file
    batch_size = args.batch_size
    window_size_min = args.window_size_min
    window_size_max = args.window_size_max
    limit_seamounts = args.limit_seamounts
    enable_dateline_stitching = args.enable_dateline_stitching

    # Contour levels: 100% (base), then 90%, 75%, 50%, 25%
    contour_percentages = [100, 90, 75, 50, 25]

    # Validate window size parameters
    if window_size_min >= window_size_max:
        raise ValueError("window_size_min must be less than window_size_max")
    if window_size_min <= 0:
        raise ValueError("window_size_min must be positive")

    print(f"Configuration:")
    print(f"  DEM file: {file_dem}")
    print(f"  GPKG file: {gpkg_file}")
    print(f"  Batch size: {batch_size}")
    print(f"  Workers: {num_workers}")
    print(f"  Window size: {window_size_min}° to {window_size_max}° (adaptive)")
    print(f"  Contour levels: {contour_percentages}% of prominence")
    print(f"  Dateline stitching: {'ENABLED' if enable_dateline_stitching else 'DISABLED'} (longitude only)")
    print(f"  GDAL cache max: {os.environ['GDAL_CACHEMAX']} MB per worker")
    print(f"  GDAL threads per worker: {os.environ['GDAL_NUM_THREADS']}")
    print(f"  Transformer caching: ENABLED (LRU cache)")
    if limit_seamounts:
        print(f"  Limit: {limit_seamounts} seamounts")

    # Load data
    print("\nLoading seamounts...")
    df = gpd.read_file(gpkg_file)

    if limit_seamounts is not None:
        df = df.head(limit_seamounts)
        print(f"Processing {len(df)} seamounts")
    else:
        print(f"Processing {len(df)} seamounts")

    # Validate dataset
    print("Validating dataset...")
    file_ext = os.path.splitext(file_dem.lower())[1]

    if file_ext == '.nc':
        ds_test = xr.open_dataset(file_dem)
        ds_test.close()
    elif file_ext in ['.tif', '.tiff', '.vrt']:
        with rasterio.open(file_dem) as src:
            print(f"  Format: {src.driver}")
            print(f"  Shape: {src.shape}")
    else:
        raise ValueError(f"Unsupported format: {file_ext}")

    print("Dataset validated")

    # Determine output file path
    output_file = os.path.join(os.path.dirname(gpkg_file), f"{name_file}.gpkg")

    # Process in batches with incremental writes
    total_batches = (len(df) + batch_size - 1) // batch_size
    first_batch = True
    total_contours_processed = 0
    all_error_peak_ids = []  # Track all peaks with contour nesting errors across all batches
    all_nocontour_peak_ids = []  # Track all peaks that produced no contours (other failures)

    for batch_num in range(total_batches):
        start_idx = batch_num * batch_size
        end_idx = min(start_idx + batch_size, len(df))
        batch_df = df.iloc[start_idx:end_idx]

        print(f"\n--- Batch {batch_num + 1}/{total_batches} ---")
        print(f"Processing peaks {start_idx} to {end_idx-1}")

        # BATCH-SPECIFIC results list (cleared each batch)
        batch_contour_results = []
        batch_error_peak_ids = []  # Track peaks with contour errors

        # Prepare data
        peak_data_list = []
        output_dir = os.path.dirname(gpkg_file)  # Get output directory from gpkg path
        for i, (idx, row) in enumerate(batch_df.iterrows()):
            # Use the peak_id from the input geopackage if it exists, otherwise use sequential numbering
            if 'peak_id' in row and not pd.isna(row['peak_id']):
                peak_id = int(row['peak_id'])
            else:
                # Fallback to sequential position if peak_id doesn't exist
                peak_id = start_idx + i + 1
            peak_data_list.append((
                peak_id, row['longitude'], row['latitude'],
                row['depth'], row['prominence'],
                file_dem, window_size_min, window_size_max, contour_percentages,
                enable_dateline_stitching, output_dir
            ))

        # Process in parallel
        with ProcessPoolExecutor(max_workers=num_workers) as executor:
            future_to_idx = {executor.submit(crop_and_process_peak, peak_data): peak_data[0]
                           for peak_data in peak_data_list}

            for future in tqdm(as_completed(future_to_idx), total=len(peak_data_list),
                              desc=f"Batch {batch_num + 1}"):
                try:
                    idx, contour_results, has_error = future.result()
                    if has_error:
                        # Track peaks with contour nesting errors
                        batch_error_peak_ids.append(idx)
                    elif contour_results:
                        batch_contour_results.extend(contour_results)
                except MemoryLimitExceeded as e:
                    # CRITICAL: Memory limit reached - stop everything
                    print(f"\n{'='*60}")
                    print(f"MEMORY LIMIT REACHED: {e}")
                    print(f"STOP AND RE-RUN WITH LOWER --num_workers")
                    print(f"{'='*60}")
                    # Shutdown executor immediately
                    executor.shutdown(wait=False, cancel_futures=True)
                    import sys
                    sys.exit(1)
                except Exception as e:
                    print(f"Error processing peak: {e}")

        print(f"Batch {batch_num + 1} completed - found {len(batch_contour_results)} contours")

        # Track peaks with errors (nesting errors)
        if batch_error_peak_ids:
            print(f"WARNING: {len(batch_error_peak_ids)} peaks had contour nesting errors (contours discarded)")
            all_error_peak_ids.extend(batch_error_peak_ids)

        # Track peaks that produced NO contours (not nesting errors, just failed to extract)
        batch_peak_ids_with_contours = set(result['peak_id'] for result in batch_contour_results)
        batch_peak_ids_processed = set(peak_data[0] for peak_data in peak_data_list)
        batch_peak_ids_without_contours = batch_peak_ids_processed - batch_peak_ids_with_contours - set(batch_error_peak_ids)

        if batch_peak_ids_without_contours:
            print(f"WARNING: {len(batch_peak_ids_without_contours)} peaks produced no contours (extraction failed)")
            all_nocontour_peak_ids.extend(batch_peak_ids_without_contours)

        # WRITE BATCH RESULTS TO DISK IMMEDIATELY
        if batch_contour_results:
            print(f"Writing batch {batch_num + 1} results to disk...")

            # Convert to DataFrame
            contour_data = []
            for result in batch_contour_results:
                # Look up peak by peak_id if it exists in the dataframe
                if 'peak_id' in df.columns:
                    peak_row = df[df['peak_id'] == result['peak_id']].iloc[0]
                else:
                    # Fallback to iloc for position-based indexing (peak_id is 1-based, iloc is 0-based)
                    peak_row = df.iloc[result['peak_id'] - 1]

                # Create display label for QGIS
                display_label = f"feature id {result['peak_id']}, ctr {result['prominence_percentage']:.0f}%"

                contour_data.append({
                    'contour name': display_label,
                    'feature_id': result['peak_id'],
                    'peak_id': result['peak_id'],  # Keep peak_id for compatibility
                    'nested_on_feature_id': '',  # Empty - will be populated by find_ctrs_cluster.py
                    'peak_longitude': peak_row['longitude'],
                    'peak_latitude': peak_row['latitude'],
                    'peak_depth': peak_row['depth'],
                    'peak_prominence': peak_row['prominence'],
                    'contour_depth': result['depth'],
                    'prominence_percentage': result['prominence_percentage'],
                    'area_sq_km': result['area_sq_km'],
                    'circularity_percent': result['circularity_percent'],
                    'bbox_orientation_deg': result['bbox_orientation_deg'],
                    'bbox_length_m': result['bbox_length_m'],
                    'bbox_width_m': result['bbox_width_m'],
                    'mean_slope_deg': result['mean_slope_deg'],
                    'min_slope_deg': result['min_slope_deg'],
                    'max_slope_deg': result['max_slope_deg'],
                    'centroid_lon': result['centroid_lon'],
                    'centroid_lat': result['centroid_lat'],
                    'window_size_used': result['window_size_used'],
                    'geometry': result['polygon']
                })

            # Create GeoDataFrame with proper type conversion
            gdf_batch = gpd.GeoDataFrame(contour_data, crs="EPSG:4326")

            # Sort with base contours (100% level) at bottom, building up to 25% at top
            # Within each seamount (feature_id), sort from base (100%) to top (25%)
            gdf_batch = gdf_batch.sort_values(['feature_id', 'prominence_percentage'], ascending=[True, False])

            # Reset index
            gdf_batch = gdf_batch.reset_index(drop=True)

            # Round coordinates and float values to reduce file size
            # Longitude/latitude coordinates: 5 decimal places (~1.1m precision)
            # Other float values: 2 decimal places
            lon_lat_cols = ['peak_longitude', 'peak_latitude', 'centroid_lon', 'centroid_lat']
            float_cols_2decimals = ['peak_depth', 'peak_prominence', 'contour_depth',
                                   'area_sq_km', 'circularity_percent', 'bbox_orientation_deg',
                                   'bbox_length_m', 'bbox_width_m', 'mean_slope_deg',
                                   'min_slope_deg', 'max_slope_deg', 'window_size_used']

            for col in lon_lat_cols:
                if col in gdf_batch.columns:
                    gdf_batch[col] = gdf_batch[col].round(5)

            for col in float_cols_2decimals:
                if col in gdf_batch.columns:
                    gdf_batch[col] = gdf_batch[col].round(2)

            # Round geometry coordinates (polygon vertices) to 5 decimal places
            def round_coords(geom):
                if geom is None:
                    return geom
                return Polygon([(round(x, 5), round(y, 5)) for x, y in geom.exterior.coords])

            gdf_batch['geometry'] = gdf_batch['geometry'].apply(round_coords)

            # Convert data types for GeoPackage compatibility
            for col in gdf_batch.columns:
                if col != gdf_batch.geometry.name:
                    dtype = gdf_batch[col].dtype
                    if dtype.kind in ['i', 'u']:
                        gdf_batch[col] = gdf_batch[col].astype('int64')
                    elif dtype.kind == 'f':
                        gdf_batch[col] = gdf_batch[col].astype('float64')
                        gdf_batch[col] = gdf_batch[col].replace([np.inf, -np.inf], np.nan)

            # Write or append to file
            if first_batch:
                # First batch: create new file
                gdf_batch.to_file(output_file, driver="GPKG")
                first_batch = False
                print(f"Created output file: {output_file}")
            else:
                # Subsequent batches: append to existing file
                gdf_batch.to_file(output_file, driver="GPKG", mode='a')
                print(f"Appended to output file: {output_file}")

            total_contours_processed += len(batch_contour_results)

            # CLEAR batch results from memory
            del batch_contour_results
            del contour_data
            del gdf_batch
            gc.collect()

            print(f"Batch {batch_num + 1} written and cleared from memory")
        else:
            print(f"No contours found in batch {batch_num + 1}")

    print("\n" + "="*60)
    print("All batches completed!")
    print("="*60)

    # Update original dataframe with error flags and save back to the ORIGINAL peaks file
    print(f"\nUpdating original peaks file with error information...")

    # Initialize error columns
    df['error_with_contours'] = False
    df['error_no_contours'] = False

    # Mark peaks with contour nesting errors
    if all_error_peak_ids:
        print(f"  Marking {len(all_error_peak_ids)} peaks with contour nesting errors...")
        if 'peak_id' in df.columns:
            for peak_id in all_error_peak_ids:
                df.loc[df['peak_id'] == peak_id, 'error_with_contours'] = True
        else:
            for peak_id in all_error_peak_ids:
                if 0 <= peak_id - 1 < len(df):
                    df.loc[df.index[peak_id - 1], 'error_with_contours'] = True

    # Mark peaks that produced no contours (extraction failures)
    if all_nocontour_peak_ids:
        print(f"  Marking {len(all_nocontour_peak_ids)} peaks that produced no contours...")
        if 'peak_id' in df.columns:
            for peak_id in all_nocontour_peak_ids:
                df.loc[df['peak_id'] == peak_id, 'error_no_contours'] = True
        else:
            for peak_id in all_nocontour_peak_ids:
                if 0 <= peak_id - 1 < len(df):
                    df.loc[df.index[peak_id - 1], 'error_no_contours'] = True

    # Save updated peaks file back to the ORIGINAL location (overwrite)
    print(f"  Saving updated peaks file to: {gpkg_file}")
    df.to_file(gpkg_file, driver="GPKG")
    print(f"  Successfully updated original peaks file with error attributes")

    # Load final results and print summary statistics
    if os.path.exists(output_file):
        print(f"\nLoading final results for summary statistics...")
        gdf_contours = gpd.read_file(output_file)

        print(f"\nSummary Statistics:")
        print(f"  Total peaks processed: {len(df)}")
        print(f"  Total contours found: {len(gdf_contours)}")
        print(f"  Features with contours: {gdf_contours['feature_id'].nunique()}")

        print(f"\nContours by level:")
        for pct in sorted(gdf_contours['prominence_percentage'].unique(), reverse=True):
            count = len(gdf_contours[gdf_contours['prominence_percentage'] == pct])
            print(f"  {pct:3.0f}%: {count:4d} contours")

        print(f"\nAverage metrics:")
        print(f"  Circularity: {gdf_contours['circularity_percent'].mean():.1f}%")
        print(f"  Mean slope: {gdf_contours['mean_slope_deg'].mean():.1f}°")
        print(f"  Min slope: {gdf_contours['min_slope_deg'].mean():.1f}°")
        print(f"  Max slope: {gdf_contours['max_slope_deg'].mean():.1f}°")

        # Window size usage statistics
        window_stats = gdf_contours['window_size_used'].value_counts().sort_index()
        print(f"\nWindow size usage:")
        for window_size, count in window_stats.items():
            percentage = (count / len(gdf_contours)) * 100
            print(f"  {window_size:4.1f}°: {count:4d} contours ({percentage:5.1f}%)")
        print(f"  Average window size used: {gdf_contours['window_size_used'].mean():.1f}°")

        print(f"\nFinal output saved to: {output_file}")

        # Report error statistics if any
        if all_error_peak_ids or all_nocontour_peak_ids:
            print(f"\n{'='*60}")
            print(f"ERROR SUMMARY:")
            print(f"{'='*60}")

            if all_error_peak_ids:
                print(f"\nContour Nesting Errors (error_with_contours = True):")
                print(f"  Total peaks: {len(all_error_peak_ids)}")
                print(f"  Percentage: {(len(all_error_peak_ids) / len(df)) * 100:.2f}%")
                print(f"  Description: Invalid contour nesting detected (contours discarded)")

            if all_nocontour_peak_ids:
                print(f"\nNo Contours Extracted (error_no_contours = True):")
                print(f"  Total peaks: {len(all_nocontour_peak_ids)}")
                print(f"  Percentage: {(len(all_nocontour_peak_ids) / len(df)) * 100:.2f}%")
                print(f"  Description: Failed to extract any contours (invalid data, all NaN, or other errors)")

            total_errors = len(all_error_peak_ids) + len(all_nocontour_peak_ids)
            print(f"\nTotal peaks with errors: {total_errors}")
            print(f"Total percentage: {(total_errors / len(df)) * 100:.2f}%")
            print(f"Peaks with successful contours: {gdf_contours['feature_id'].nunique()}")
            print(f"\nError flags have been added to: {gpkg_file}")

            # Verify accounting
            peaks_with_contours = gdf_contours['feature_id'].nunique()
            peaks_with_errors = total_errors
            total_accounted = peaks_with_contours + peaks_with_errors
            print(f"\n{'='*60}")
            print(f"ACCOUNTING VERIFICATION:")
            print(f"{'='*60}")
            print(f"  Total peaks processed: {len(df)}")
            print(f"  Peaks with contours: {peaks_with_contours}")
            print(f"  Peaks with nesting errors: {len(all_error_peak_ids)}")
            print(f"  Peaks with no contours: {len(all_nocontour_peak_ids)}")
            print(f"  Total accounted: {total_accounted}")
            print(f"  Difference: {len(df) - total_accounted}")
            if len(df) == total_accounted:
                print(f"  ✓ All peaks accounted for!")
            else:
                print(f"  ✗ WARNING: Accounting mismatch!")
    else:
        print("No valid contours found")

    print("\nDone!")


if __name__ == "__main__":
    main()
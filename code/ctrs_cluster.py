#!/usr/bin/env python3
"""
Contour nesting analysis: Add nested_on_feature_id attribute to contour geopackage
Optimized with tile-based processing around analyzed peaks
"""

import argparse
import geopandas as gpd
import numpy as np
from tqdm import tqdm
from shapely.geometry import box, Polygon
try:
    from shapely import STRtree
except ImportError:
    from shapely.strtree import STRtree
import os


def normalize_longitude(lon):
    """
    Normalize longitude to [-180, 180] range
    """
    return ((lon + 180) % 360) - 180




def bbox_contains(parent_bbox, child_bbox):
    """
    Fast bounding box containment check.
    Returns True if child_bbox is completely contained within parent_bbox.
    Handles longitude normalization but not dateline crossing.
    """
    parent_min_lon, parent_min_lat, parent_max_lon, parent_max_lat = parent_bbox
    child_min_lon, child_min_lat, child_max_lon, child_max_lat = child_bbox

    # Check if child bbox is contained within parent bbox
    return (parent_min_lon <= child_min_lon and
            child_max_lon <= parent_max_lon and
            parent_min_lat <= child_min_lat and
            child_max_lat <= parent_max_lat)


def detect_contour_nesting(gdf_contours, tile_buffer_degrees=30.0):
    """
    Detect which features are nested within other feature structures.
    Uses spatial indexing with bounding box containment checks for efficiency.
    """
    print("Detecting feature nesting relationships using spatial indexing...")

    # Get unique features with their largest contour (base contour)
    print("Building feature index...")
    feature_data = {}
    feature_geometries = []
    feature_ids = []

    for idx, row in tqdm(gdf_contours.iterrows(), total=len(gdf_contours), desc="Processing contours"):
        feature_id = row['feature_id']
        area = row['area_sq_km']

        if feature_id not in feature_data:
            feature_data[feature_id] = {
                'largest_area': area,
                'largest_contour_idx': idx,
                'bbox': None,
                'bbox_geometry': None,
                'all_contour_indices': [idx]
            }
        else:
            feature_data[feature_id]['all_contour_indices'].append(idx)
            if area > feature_data[feature_id]['largest_area']:
                feature_data[feature_id]['largest_area'] = area
                feature_data[feature_id]['largest_contour_idx'] = idx

    # Calculate bounding boxes and create geometries for spatial index
    print("Building spatial index...")
    for feature_id, data in tqdm(feature_data.items(), desc="Computing bboxes"):
        largest_idx = data['largest_contour_idx']
        bounds = gdf_contours.loc[largest_idx, 'geometry'].bounds
        min_lon = normalize_longitude(bounds[0])
        max_lon = normalize_longitude(bounds[2])

        # Handle dateline crossing
        if min_lon > max_lon:
            data['bbox'] = bounds  # Keep original bounds for dateline crossing
            bbox_geom = box(bounds[0], bounds[1], bounds[2], bounds[3])
        else:
            data['bbox'] = (min_lon, bounds[1], max_lon, bounds[3])
            bbox_geom = box(min_lon, bounds[1], max_lon, bounds[3])

        data['bbox_geometry'] = bbox_geom
        feature_geometries.append(bbox_geom)
        feature_ids.append(feature_id)

    # Build spatial index
    spatial_index = STRtree(feature_geometries)

    # Sort features by area (largest first) for hierarchical nesting
    sorted_features = sorted(feature_data.items(), key=lambda x: x[1]['largest_area'], reverse=True)

    print(f"Checking nesting for {len(sorted_features)} features...")
    nesting_updates = []
    comparisons = 0
    bbox_matches = 0
    polygon_confirmations = 0

    # Check nesting using spatial index + bounding box containment
    for parent_id, parent_data in tqdm(sorted_features, desc="Finding nested features"):
        parent_bbox = parent_data['bbox']
        parent_bbox_geom = parent_data['bbox_geometry']

        # Query spatial index for potential overlaps
        potential_children_indices = spatial_index.query(parent_bbox_geom)

        # Check each potential child
        for child_idx in potential_children_indices:
            if child_idx < len(feature_ids):
                child_feature_id = feature_ids[child_idx]

                # Skip self and features with larger or equal area
                if (child_feature_id == parent_id or
                    feature_data[child_feature_id]['largest_area'] >= parent_data['largest_area']):
                    continue

                child_bbox = feature_data[child_feature_id]['bbox']
                comparisons += 1

                # Fast bounding box containment check
                if bbox_contains(parent_bbox, child_bbox):
                    bbox_matches += 1
                    # Confirm with actual polygon containment check
                    parent_geom = gdf_contours.loc[parent_data['largest_contour_idx'], 'geometry']
                    child_geom = gdf_contours.loc[feature_data[child_feature_id]['largest_contour_idx'], 'geometry']

                    # Handle longitude normalization for polygon containment
                    try:
                        # Normalize parent polygon coordinates
                        parent_coords = list(parent_geom.exterior.coords)
                        normalized_parent_coords = []
                        for lon, lat in parent_coords:
                            norm_lon = normalize_longitude(lon)
                            normalized_parent_coords.append((norm_lon, lat))

                        normalized_parent = Polygon(normalized_parent_coords)
                        if not normalized_parent.is_valid:
                            normalized_parent = normalized_parent.buffer(0)

                        # Normalize child polygon coordinates
                        child_coords = list(child_geom.exterior.coords)
                        normalized_child_coords = []
                        for lon, lat in child_coords:
                            norm_lon = normalize_longitude(lon)
                            normalized_child_coords.append((norm_lon, lat))

                        normalized_child = Polygon(normalized_child_coords)
                        if not normalized_child.is_valid:
                            normalized_child = normalized_child.buffer(0)

                        # Actual polygon containment check
                        if normalized_parent.contains(normalized_child):
                            polygon_contained = True
                        else:
                            polygon_contained = False

                    except Exception:
                        # Fallback to original polygons if normalization fails
                        try:
                            polygon_contained = parent_geom.contains(child_geom)
                        except Exception:
                            polygon_contained = False

                    # Only add nesting if polygon containment is confirmed
                    if polygon_contained:
                        polygon_confirmations += 1
                        # Add nesting relationship to all contours of child feature
                        for contour_idx in feature_data[child_feature_id]['all_contour_indices']:
                            current_nested = gdf_contours.at[contour_idx, 'nested_on_feature_id']
                            if current_nested:
                                new_nested = f"{current_nested},{parent_id}"
                            else:
                                new_nested = str(parent_id)

                            nesting_updates.append((contour_idx, new_nested))

    # Apply updates efficiently using bulk operations
    print(f"Applying {len(nesting_updates)} nesting relationships...")
    if nesting_updates:
        # Convert to dict for bulk update
        update_dict = dict(nesting_updates)
        indices = list(update_dict.keys())
        values = list(update_dict.values())

        # Bulk update using loc
        gdf_contours.loc[indices, 'nested_on_feature_id'] = values

    print(f"  Performed {comparisons:,} bounding box comparisons")
    print(f"  Found {bbox_matches:,} bbox containments")
    print(f"  Confirmed {polygon_confirmations:,} polygon containments")
    print(f"  Found {len(nesting_updates)} nesting relationships")
    if bbox_matches > 0:
        confirmation_rate = (polygon_confirmations / bbox_matches) * 100
        print(f"  Polygon confirmation rate: {confirmation_rate:.1f}%")

    return gdf_contours


def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description='Add nesting analysis to contour geopackage')

    parser.add_argument('--input_gpkg', type=str, required=True,
                       help='Path to the input contour GeoPackage file')
    parser.add_argument('--output_gpkg', type=str, default=None,
                       help='Path to the output GeoPackage file (default: overwrites input)')
    parser.add_argument('--tile_buffer_degrees', type=float, default=45.0,
                       help='Tile buffer size in degrees for optimization (default: 30.0)')

    return parser.parse_args()


def main():
    args = parse_arguments()

    input_gpkg = args.input_gpkg
    output_gpkg = args.output_gpkg or input_gpkg
    tile_buffer_degrees = args.tile_buffer_degrees

    print(f"Configuration:")
    print(f"  Input GPKG: {input_gpkg}")
    print(f"  Output GPKG: {output_gpkg}")
    print(f"  Tile buffer: {tile_buffer_degrees}°")

    # Validate input file
    if not os.path.exists(input_gpkg):
        raise FileNotFoundError(f"Input file not found: {input_gpkg}")

    # Load contour data
    print(f"\nLoading contour data from {input_gpkg}...")
    gdf_contours = gpd.read_file(input_gpkg)
    print(f"Loaded {len(gdf_contours)} contours for {gdf_contours['feature_id'].nunique()} features")

    # Validate required columns
    required_columns = ['feature_id', 'nested_on_feature_id', 'area_sq_km', 'geometry']
    missing_columns = [col for col in required_columns if col not in gdf_contours.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns: {missing_columns}")

    # Check if nesting analysis was already performed
    existing_nested = gdf_contours[gdf_contours['nested_on_feature_id'] != '']
    if len(existing_nested) > 0:
        print(f"Warning: {len(existing_nested)} contours already have nesting data. This will be overwritten.")
        # Reset existing nesting data
        gdf_contours['nested_on_feature_id'] = ''

    # Perform nesting analysis
    gdf_contours = detect_contour_nesting(gdf_contours, tile_buffer_degrees)

    # Convert data types for GeoPackage compatibility
    for col in gdf_contours.columns:
        if col != gdf_contours.geometry.name:
            dtype = gdf_contours[col].dtype
            if dtype.kind in ['i', 'u']:
                gdf_contours[col] = gdf_contours[col].astype('int64')
            elif dtype.kind == 'f':
                gdf_contours[col] = gdf_contours[col].astype('float64')
                gdf_contours[col] = gdf_contours[col].replace([np.inf, -np.inf], np.nan)

    # Sort contours by area for proper rendering order
    # When saved to GeoPackage, the FID will be assigned based on dataframe row order
    # Row 0 → FID 1 (largest), Row 1 → FID 2, ..., Last row → Last FID (smallest)
    # Standard rendering: FID 1 drawn first (bottom), last FID drawn last (top)
    print("\nSorting contours by area for rendering order...")

    gdf_contours = gdf_contours.sort_values('area_sq_km', ascending=False).reset_index(drop=True)

    # Verify sorting
    first_area = gdf_contours.iloc[0]['area_sq_km']
    last_area = gdf_contours.iloc[-1]['area_sq_km']
    print(f"  Sorted {len(gdf_contours)} contours by area:")
    print(f"    First (FID=1): {first_area:.2f} km² (largest)")
    print(f"    Last (FID={len(gdf_contours)}): {last_area:.2f} km² (smallest)")

    # Save updated data
    # The to_file() method will create FIDs sequentially based on row order
    print(f"\nSaving updated contours to {output_gpkg}...")
    gdf_contours.to_file(output_gpkg, driver="GPKG", compression="ZLIB")
    print(f"Saved {len(gdf_contours)} contours with FIDs sorted by area")

    # Summary statistics
    nested_contours = gdf_contours[gdf_contours['nested_on_feature_id'] != '']
    print(f"\nSummary Statistics:")
    print(f"  Total contours: {len(gdf_contours)}")
    print(f"  Contours with nesting: {len(nested_contours)}")
    print(f"  Features with contours: {gdf_contours['feature_id'].nunique()}")

    # Nesting by contour level
    print(f"\nNesting by contour level:")
    for pct in sorted(gdf_contours['prominence_percentage'].unique(), reverse=True):
        level_contours = gdf_contours[gdf_contours['prominence_percentage'] == pct]
        level_nested = level_contours[level_contours['nested_on_feature_id'] != '']
        print(f"  {pct:3.0f}%: {len(level_nested)}/{len(level_contours)} nested")

    print("Done!")


if __name__ == "__main__":
    main()
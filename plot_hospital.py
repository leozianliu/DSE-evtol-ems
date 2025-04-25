import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import csv
import os
from shapely.ops import transform
import pyproj
from functools import partial
import numpy as np
from matplotlib.patches import Circle

# Set environment variable for OSM reading
os.environ['OSM_USE_CUSTOM_INDEXING'] = 'NO'

# --- 1. Load Hospitals Data ---
hospitals = {}
name_lst = []
lat_lst = []
lon_lst = []

with open('barv_24h_hospitals.csv', newline='') as csvfile:
    reader = csv.reader(csvfile)
    next(reader)  # Skip header
    for row in reader:
        name_lst.append(row[0])
        lat_lst.append(float(row[1]))
        lon_lst.append(float(row[2]))

hospitals['Latitude'] = lat_lst
hospitals['Longitude'] = lon_lst

geometry = [Point(xy) for xy in zip(hospitals['Longitude'], hospitals['Latitude'])]
gdf_hospitals = gpd.GeoDataFrame(hospitals, geometry=geometry, crs="EPSG:4326")

# --- 2. Create Buffer Circles (in km) ---
def geodesic_point_buffer(lat, lon, km):
    # Azimuthal Equidistant Projection centered at the point
    local_proj = pyproj.Proj(proj='aeqd', lat_0=lat, lon_0=lon, datum='WGS84')
    wgs84 = pyproj.Proj('EPSG:4326')

    project_to_aeqd = partial(pyproj.transform, wgs84, local_proj)
    project_to_wgs84 = partial(pyproj.transform, local_proj, wgs84)

    # Transform point to local projection
    point = transform(project_to_aeqd, Point(lon, lat))

    # Buffer in meters, then transform back to lat/lon
    buffer = point.buffer(km * 1000)
    return transform(project_to_wgs84, buffer)


buffer_distance_km = 50
gdf_hospitals['buffer'] = gdf_hospitals.apply(
    lambda row: geodesic_point_buffer(row['Latitude'], row['Longitude'], buffer_distance_km),
    axis=1
)

# Convert buffers to GeoDataFrame
gdf_buffers = gpd.GeoDataFrame(gdf_hospitals, geometry='buffer', crs="EPSG:4326")

# --- 3. Load Boundaries ---
bavaria = gpd.read_file("bavaria.geojson")

region_files = [
    "Mittelfranken.osm", "Niederbayern.osm", "Oberbayern.osm",
    "Oberfranken.osm", "Oberpfalz.osm", "Schwaben.osm", "Unterfranken.osm"
]
regions = []

for file in region_files:
    for layer in ['multipolygons', 'lines']:
        try:
            region = gpd.read_file(file, layer=layer)
            if not region.empty:
                regions.append(region)
                break
        except:
            continue
    else:
        regions.append(gpd.GeoDataFrame())  # Empty if all fails

# --- 4. Extract Coordinates to NumPy Array ---
coords_array = np.array(list(zip(gdf_hospitals['Latitude'], gdf_hospitals['Longitude'])))
# print("Geodetic hospital coordinates (Lat, Lon):")
# print(coords_array)

# --- 5. Plot ---
fig, ax = plt.subplots(figsize=(12, 8))

# Plot Bavaria boundary
bavaria.plot(ax=ax, color='lightblue', edgecolor='black', alpha=0.3)

# Plot hospital points
gdf_hospitals.plot(ax=ax, marker='o', color='red', markersize=50)

# Option 1: Plot geodesic buffers (accurate, shapely-based)
gdf_buffers.plot(ax=ax, facecolor='purple', edgecolor='purple', linewidth=1, alpha=0.1)

# Customize plot
plt.title(f'24h Hospitals in Bavaria with {buffer_distance_km}km Radius')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True, alpha=0.3)

# Zoom to Bavaria bounds
ax.set_xlim(bavaria.total_bounds[[0, 2]])
ax.set_ylim(bavaria.total_bounds[[1, 3]])
ax.set_aspect('equal')

plt.show()

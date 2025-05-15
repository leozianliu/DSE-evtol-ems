import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import csv
import os
from shapely.ops import transform
import pyproj
from functools import partial
import numpy as np
import contextily as ctx

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
hospitals['Name'] = name_lst

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

buffer_distance_km = 60
gdf_hospitals['buffer'] = gdf_hospitals.apply(
    lambda row: geodesic_point_buffer(row['Latitude'], row['Longitude'], buffer_distance_km),
    axis=1
)

# Convert buffers to GeoDataFrame
gdf_buffers = gpd.GeoDataFrame(gdf_hospitals, geometry='buffer', crs="EPSG:4326")

# --- 3. Load Boundaries ---
bavaria = gpd.read_file("bavaria.geojson")

# --- 4. Create City Points ---
# Example cities with hospitals (you should replace with actual cities from your data)
# Format: {'City': [lat, lon]}
major_cities = {
# Major cities
    'Munich': (48.1351, 11.5820),
    'Nuremberg': (49.4521, 11.0767),
    'Augsburg': (48.3705, 10.8978),
    'Regensburg': (49.0134, 12.1016),
    'Würzburg': (49.7913, 9.9534),
    'Ingolstadt': (48.7665, 11.4258),
    'Fürth': (49.4771, 10.9887),
    'Erlangen': (49.5895, 11.0140),
    'Bayreuth': (49.9480, 11.5783),
    'Bamberg': (49.8917, 10.8850),
    'Aschaffenburg': (49.9738, 9.1526),
    'Kempten': (47.7345, 10.3146),
    'Rosenheim': (47.8564, 12.1284),
    'Landshut': (48.5393, 12.1507),
    'Passau': (48.5667, 13.4333),
    
    # Additional larger towns
    'Hof': (50.3167, 11.9167),
    'Schweinfurt': (50.0500, 10.2333),
    'Straubing': (48.8833, 12.5667),
    'Dachau': (48.2600, 11.4340),
    'Freising': (48.4000, 11.7500),
    'Coburg': (50.2667, 10.9667),
    'Memmingen': (47.9833, 10.1833),
    'Neu-Ulm': (48.4000, 10.0000),
    'Amberg': (49.4500, 11.8500),
    'Weiden': (49.6769, 12.1561),

    # Regional centers
    'Garmisch-Partenkirchen': (47.5000, 11.1000),
    'Berchtesgaden': (47.6333, 13.0000),
    'Bad Tölz': (47.7667, 11.5500),
    'Mühldorf': (48.2500, 12.5333),
    'Altötting': (48.2333, 12.6833),
    'Deggendorf': (48.8333, 12.9667),
    'Kelheim': (48.9167, 11.8667),
    'Eichstätt': (48.8833, 11.1833),
    'Donauwörth': (48.7167, 10.8000),
    'Ansbach': (49.3000, 10.5833),
    'Dinkelsbühl': (49.0667, 10.3167),
    'Rothenburg ob der Tauber': (49.3833, 10.1833),
    'Wunsiedel': (50.0333, 12.0167),
    'Pegnitz': (49.7500, 11.5333),
    'Kronach': (50.2500, 11.3333),
}

# Create GeoDataFrame for cities
city_data = {'City': [], 'Latitude': [], 'Longitude': []}
for city, coords in major_cities.items():
    city_data['City'].append(city)
    city_data['Latitude'].append(coords[0])
    city_data['Longitude'].append(coords[1])

geometry = [Point(xy) for xy in zip(city_data['Longitude'], city_data['Latitude'])]
gdf_cities = gpd.GeoDataFrame(city_data, geometry=geometry, crs="EPSG:4326")

# --- 5. Plot ---
fig, ax = plt.subplots(figsize=(12, 8))

# Plot Bavaria boundary
bavaria.plot(ax=ax, color='lightblue', edgecolor='black', alpha=0.3)

# Plot hospital points (larger and red)
gdf_hospitals.plot(ax=ax, marker='o', color='red', markersize=50, label='24h Hospitals')

# Plot city points (smaller and blue)
gdf_cities.plot(ax=ax, marker='o', color='blue', markersize=30, label='City Centers')

# Add city labels (offset slightly for better visibility)
for x, y, label in zip(gdf_cities.geometry.x, gdf_cities.geometry.y, gdf_cities['City']):
    ax.text(x + 0.02, y + 0.02, label, fontsize=9, ha='left', va='bottom')

# Plot geodesic buffers
gdf_buffers.plot(ax=ax, facecolor='purple', edgecolor='purple', linewidth=1, alpha=0.1)

# Add satellite map basemap
ctx.add_basemap(ax, crs=gdf_hospitals.crs.to_string(), source=ctx.providers.Esri.WorldImagery, alpha=0.7)

# Customize plot
plt.title(f'24h Hospitals and City Centers in Bavaria with {buffer_distance_km}km Radius')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True, alpha=0.3)
plt.legend()

# Zoom to Bavaria bounds
ax.set_xlim(bavaria.total_bounds[[0, 2]])
ax.set_ylim(bavaria.total_bounds[[1, 3]])
ax.set_aspect('equal')

plt.show()
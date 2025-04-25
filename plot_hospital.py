import geopandas as gpd
import matplotlib.pyplot as plt
from shapely.geometry import Point
import csv

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
gdf = gpd.GeoDataFrame(hospitals, geometry=geometry, crs="EPSG:4326")

bavaria = gpd.read_file("bavaria.geojson")  # Make sure this file exists

fig, ax = plt.subplots(figsize=(12, 8))

# Plot Bavaria boundary
bavaria.plot(ax=ax, color='lightblue', edgecolor='black')

# Plot hospital points
gdf.plot(ax=ax, marker='o', color='red', markersize=50)

# Set map style
plt.title('24h Hospitals in Bavaria')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.grid(True)

# Optional: Zoom to Bavaria bounds
ax.set_xlim(bavaria.total_bounds[[0, 2]])
ax.set_ylim(bavaria.total_bounds[[1, 3]])

plt.show()

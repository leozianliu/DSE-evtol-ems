import numpy as np
import pandas as pd

import numpy as np
import pandas as pd

import geopandas as gpd
import matplotlib.pyplot as plt

export = True
plot = True
normalize = True

def read_excel_to_numpy(file_path, sheet_name, column_range, row_range):
    """
    Reads a specific range of data from an Excel file and converts it to a NumPy array.

    Parameters:
    - file_path: str, path to the Excel file
    - sheet_name: str, name of the sheet to read from
    - column_range: str, range of columns to read (e.g., "A15:A2153")

    Returns:
    - numpy_array: np.ndarray, the data as a NumPy array
    """
    # Read the data from the specified range in the Excel file
    data = pd.read_excel(file_path, sheet_name=sheet_name, usecols=column_range, skiprows=row_range[0]-1, nrows=row_range[1]-row_range[0]+1)

    # Convert the DataFrame to a NumPy array
    numpy_array = data.to_numpy()

    return numpy_array

population_regions = read_excel_to_numpy('data/population_age_data.xlsx', 'Sheet 1', 'B', (15, 2153))
population_names = read_excel_to_numpy('data/population_age_data.xlsx', 'Sheet 1', 'A', (15, 2153))
healthcare_regions = read_excel_to_numpy('data/RYB2024_CH02_Health.xlsx', 'CH02M1', 'B', (2, 1199))

# Find indices where population_regions values are in healthcare_regions
matching_indices = np.where(np.isin(population_regions, healthcare_regions))[0]

population = read_excel_to_numpy('data/population_age_data.xlsx', 'Sheet 2', 'C', (15, 2153))
population = np.hstack((population, read_excel_to_numpy('data/population_age_data.xlsx', 'Sheet 3', 'C', (15, 2153))))
population = np.hstack((population, read_excel_to_numpy('data/population_age_data.xlsx', 'Sheet 4', 'C', (15, 2153))))
# Filter matching_indices to include only rows in population that contain integers
integer_indices = matching_indices[np.all(np.char.isnumeric(population[matching_indices].astype(str)), axis=1)]

matching_healthcare_indices = np.array([
    np.where(healthcare_regions == value)[0][0] if value in healthcare_regions else -1
    for value in population_regions[integer_indices]
])


healthcare = read_excel_to_numpy('data\RYB2024_CH02_Health.xlsx', 'CH02M1', 'C', (2, 1199))

# Annual emergencies per person by age group (<15, 15-64, 65+)
emergency_age = np.array([0.0293, 0.0717, 0.3626])
emergencies = np.sum(population[integer_indices] * emergency_age, axis=1, keepdims=True)
# Calculate unaccessible as a column vector
unaccessible = 1 - healthcare[matching_healthcare_indices].astype(float) * 0.01
# Calculate crit_ems with the desired shape (1145, 1)
crit_ems = (emergencies * unaccessible)[:, 0]

if normalize:
    crit_ems = crit_ems / np.sum(population[integer_indices], axis=1)

sorted_indices = np.argsort(crit_ems)
sorted_population_regions = population_names[integer_indices][sorted_indices]

# print(sorted_population_regions, crit_ems[sorted_indices])

if export:
    # Create a DataFrame with sorted_population_regions and crit_ems[sorted_indices]
    export_data = pd.DataFrame({
        'Region': sorted_population_regions.flatten(),
        'NUTS3 Code': population_regions[integer_indices][sorted_indices].flatten(),
        'Critical Emergencies': crit_ems[sorted_indices]
    })

    # Export the DataFrame to an Excel file
    export_data.to_excel('sorted_population_regions.xlsx', index=False)

    print("Data exported to 'sorted_population_regions.xlsx'")

if plot:
    # Load the shapefile or GeoJSON containing NUTS3 regions
    # Replace 'path_to_nuts3_shapefile' with the actual path to your shapefile or GeoJSON
    nuts3_map = gpd.read_file('data/NUTS_RG_20M_2024_3035.geojson')

    # Create a DataFrame with population_regions[integer_indices] and crit_ems[sorted_indices]
    region_data = pd.DataFrame({
        'region': population_regions[integer_indices].flatten(),
        'crit_ems': crit_ems
    })

    # Merge the shapefile data with your region data
    # Assuming the NUTS3 codes in the shapefile are stored in a column named 'NUTS_ID'
    merged_map = nuts3_map.merge(region_data, left_on='NUTS_ID', right_on='region', how='inner')

    # Plot the map
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))
    merged_map.plot(column='crit_ems', cmap='Reds', legend=True, ax=ax)

    # Add a title
    ax.set_title('NUTS3 Regions Colored by Annual Emergencies not reached in 15 Minutes')

    # Customize the legend to include the last 20 elements of sorted_population_regions
    legend_labels = [f"{sorted_population_regions[-x, 0]} - {round(crit_ems[sorted_indices][-x])}" for x in range(1, 31)]

    # Add the custom legend
    handles = [plt.Line2D([0], [0], color='red', lw=4, label=label) for label in legend_labels]
    ax.legend(handles=handles, title="Worst 30 Regions", loc='lower left', bbox_to_anchor=(1, 0))

    # Show the plot
    plt.tight_layout()
    plt.show()
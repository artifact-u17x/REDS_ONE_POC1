# Declaration: as this script produces non-critical information,
# it was vibe-coded less careful human checking than the analysis code.
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
from collections import defaultdict

# Read the CSV file
df = pd.read_csv('garden_coordinates_rural_urban_region.csv')


# Create Rural/Urban/Unknown category based on RUC21NM column
def categorize_location(ruc_name):
    if pd.isna(ruc_name):
        return 'Unknown'
    ruc_lower = str(ruc_name).lower()
    if 'rural' in ruc_lower:
        return 'Rural'
    elif 'urban' in ruc_lower:
        return 'Urban'
    else:
        return 'Unknown'


df['Location_Type'] = df['RUC21NM'].apply(categorize_location)

# Extract coordinates
x = df['epsg3857.x'].values
y = df['epsg3857.y'].values

# Snap to 5km grid
grid_size = 5000
x_snapped = np.round(x / grid_size) * grid_size
y_snapped = np.round(y / grid_size) * grid_size

# Check for duplicate snapped locations and handle them
snapped_coords = list(zip(x_snapped, y_snapped))
unique_coords = set(snapped_coords)

if len(snapped_coords) != len(unique_coords):
    num_duplicates = len(snapped_coords) - len(unique_coords)
    print(f"WARNING: Snapping resulted in {num_duplicates} duplicate location(s)")
    print(f"Original locations: {len(snapped_coords)}, Unique snapped locations: {len(unique_coords)}")

    # Group indices by coordinate
    coord_to_indices = defaultdict(list)
    for idx, coord in enumerate(snapped_coords):
        coord_to_indices[coord].append(idx)

    # Count overlaps
    overlapping = {coord: indices for coord, indices in coord_to_indices.items() if len(indices) > 1}
    print(f"Number of grid cells with multiple points: {len(overlapping)}")
    max_overlap = max(len(indices) for indices in coord_to_indices.values())
    print(f"Maximum points in a single grid cell: {max_overlap}")

    # Shift overlapping points
    points_shifted = 0
    for coord, indices in overlapping.items():
        if len(indices) == 2:
            # Shift first point north, second point south
            y_snapped[indices[0]] += grid_size
            y_snapped[indices[1]] -= grid_size
            points_shifted += 2
        elif len(indices) > 2:
            print(f"WARNING: {len(indices)} points at location {coord} - only handling pairs")
            # Still shift first two as a pair
            y_snapped[indices[0]] += grid_size
            y_snapped[indices[1]] -= grid_size
            points_shifted += 2

    print(f"Shifted {points_shifted} points to avoid overlaps")

# Create the plot - A4 quarter-page size (approximately 4" × 5.2")
fig, ax = plt.subplots(figsize=(4, 5.2))

# Download UK boundary from Natural Earth (more reliable source)
url = "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries.zip"
world = gpd.read_file(url)

# Filter for UK and transform to EPSG:3857
uk = world[world['ADMIN'] == 'United Kingdom']
uk = uk.to_crs('EPSG:3857')

# Plot UK boundary
uk.boundary.plot(ax=ax, edgecolor='black', linewidth=0.8)

# Set axis limits to UK bounds
ax.set_xlim(-909000, 200000)
ax.set_ylim(6420000, 7800000)

# Plot points by category with different markers
markers = {'Urban': 's', 'Rural': '^', 'Unknown': 'o'}  # square, triangle, circle
colors = {'Urban': 'red', 'Rural': 'green', 'Unknown': 'blue'}

for location_type in ['Urban', 'Rural', 'Unknown']:
    mask = df['Location_Type'] == location_type
    if mask.any():
        ax.scatter(x_snapped[mask], y_snapped[mask],
                   c=colors[location_type],
                   marker=markers[location_type],
                   s=30, alpha=0.7,
                   edgecolors='black', linewidth=0.5,
                   label=location_type, zorder=5)

# Add legend (top right)
ax.legend(loc='upper right', frameon=True, fancybox=False,
          edgecolor='black', facecolor='white')

# Remove axis ticks
ax.set_xticks([])
ax.set_yticks([])

# Set background to white
ax.set_facecolor('white')
fig.patch.set_facecolor('white')

# Remove the frame/spines for cleaner look
for spine in ax.spines.values():
    spine.set_visible(False)

# Adjust layout and save
plt.tight_layout()
plt.savefig('uk_gardens_map.png', dpi=300, bbox_inches='tight', facecolor='white')
print(f"Map saved as 'uk_gardens_map.png'")
print(f"Locations of {len(df)} garden included in model update snapped to {grid_size / 1000}km grid")
print(f"Urban: {sum(df['Location_Type'] == 'Urban')}, "
      f"Rural: {sum(df['Location_Type'] == 'Rural')}, "
      f"Unknown: {sum(df['Location_Type'] == 'Unknown')}")

# Optional: show the plot
plt.show()
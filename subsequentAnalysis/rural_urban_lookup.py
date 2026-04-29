# Declaration: as this script produces non-critical information,
# it was vibe-coded less careful human checking than the analysis code.

import pandas as pd
import requests
import time
from pyproj import Transformer

# Set up coordinate transformer: EPSG:3857 (Web Mercator) -> EPSG:4326 (WGS84 lat/long)
transformer = Transformer.from_crs("EPSG:3857", "EPSG:4326", always_xy=True)

def get_lsoa_code(lat, lon):
    time.sleep(0.1)  # 100ms delay to be nice to the API

    url = f"https://api.postcodes.io/postcodes?lon={lon}&lat={lat}"
    try:
        response = requests.get(url)

        if response.status_code != 200:
            return None

        data = response.json()

        if not data.get('result') or len(data['result']) == 0:
            return None

        result = data['result'][0]
        codes = result.get('codes', {})
        return codes.get('lsoa')

    except Exception as e:
        return None

# Load data
print("Loading garden coordinates...")
gardens = pd.read_csv('garden_coordinates.csv')

print(f"Processing all {len(gardens)} gardens\n")

print("Loading geography lookup files...")
# Load LSOA to OA lookup
# This file is not redistributed as I haven't checked rights to do so.
# It was originally obtained from:https://www.data.gov.uk/dataset/1676bb60-bef1-43dc-8c23-641554b066bd/postcode-to-oa-2021-to-lsoa-to-msoa-to-lad-february-2025-best-fit-lookup-in-the-uk
lsoa_oa_lookup = pd.read_csv('geography_databases/PCD_OA21_LSOA21_MSOA21_LAD_FEB25_UK_LU.csv', low_memory=False)

# Load RUC file
# This file is not redistributed as I haven't checked rights to do so.
# It was originally obtained from:
# https://geoportal.statistics.gov.uk/datasets/ons::rural-urban-classification-2021-of-output-areas-in-ew/about
ruc = pd.read_csv('geography_databases/Rural_Urban_Classification_(2021)_of_Output_Areas_in_EW.csv')

print(f"  LSOA-OA lookup: {len(lsoa_oa_lookup):,} rows")
print(f"  RUC file: {len(ruc):,} Output Areas")

# Convert coordinates from EPSG:3857 to lat/long (EPSG:4326)
print("\nConverting coordinates from Web Mercator to lat/long...")
lon_lat = transformer.transform(gardens['epsg3857.x'].values, gardens['epsg3857.y'].values)
gardens['longitude'] = lon_lat[0]
gardens['latitude'] = lon_lat[1]

# Get LSOA codes
print("\nFetching LSOA codes from postcodes.io API:")
lsoa_codes = []
for idx, row in gardens.iterrows():
    guid = row['garden.guid']
    lat = row['latitude']
    lon = row['longitude']

    if (idx + 1) % 10 == 0 or idx < 5:  # Show progress every 10 gardens
        print(f"  Processing garden {idx + 1}/{len(gardens)}...")

    lsoa = get_lsoa_code(lat, lon)
    lsoa_codes.append(lsoa)

gardens['lsoa_code'] = lsoa_codes
print(f"  Completed all {len(gardens)} gardens")

# Create LSOA to OA mapping (take the most common OA for each LSOA, or just pick first)
# Actually, for rural-urban classification, we want the majority classification per LSOA
print("\nMapping LSOA codes to OA codes...")
lsoa_to_oa = lsoa_oa_lookup[['lsoa21cd', 'oa21cd']].drop_duplicates()

# Join gardens with LSOA-OA lookup to get OA codes
gardens_with_oa = gardens.merge(
    lsoa_to_oa,
    left_on='lsoa_code',
    right_on='lsoa21cd',
    how='left'
)

print(f"  Matched {gardens_with_oa['oa21cd'].notna().sum()} out of {len(gardens)} gardens to OA codes")

# DIAGNOSTIC: Check sample of matched data
print(f"\n  Sample LSOA codes from API: {gardens['lsoa_code'].head(3).tolist()}")
print(f"  Sample OA codes after join: {gardens_with_oa['oa21cd'].head(3).tolist()}")

# For each LSOA, we might have multiple OAs. Let's get all possible RUC values and take the most common
print("\nJoining with RUC data...")

# Join with RUC to get classifications for all OA codes
gardens_with_ruc = gardens_with_oa.merge(
    ruc[['OA21CD', 'RUC21NM']],
    left_on='oa21cd',
    right_on='OA21CD',
    how='left'
)

print(f"  After RUC join: {len(gardens_with_ruc)} rows")
print(f"  Rows with RUC21NM: {gardens_with_ruc['RUC21NM'].notna().sum()}")
print(f"  Sample RUC values: {gardens_with_ruc['RUC21NM'].dropna().head(3).tolist()}")

# Check if RUC21NM column exists
if 'RUC21NM' not in gardens_with_ruc.columns:
    print("\nERROR: RUC21NM column not found after merge!")
    print(f"Available columns: {gardens_with_ruc.columns.tolist()}")
    exit(1)

# Group by garden and take the most common RUC classification per LSOA
# (since one LSOA can have multiple OAs)
def get_most_common_ruc(series):
    """Get most common non-null value from a series"""
    series_clean = series.dropna()
    if len(series_clean) == 0:
        return None
    return series_clean.mode()[0] if len(series_clean.mode()) > 0 else series_clean.iloc[0]

gardens_final = gardens_with_ruc.groupby('garden.guid').agg({
    'epsg3857.x': 'first',
    'epsg3857.y': 'first',
    'lsoa_code': 'first',
    'RUC21NM': get_most_common_ruc
}).reset_index()

# Save
output = gardens_final[['garden.guid', 'epsg3857.x', 'epsg3857.y', 'RUC21NM', 'lsoa_code']]

# Replace None with informative labels based on LSOA code
def categorize_missing(row):
    if pd.notna(row['RUC21NM']):
        return row['RUC21NM']
    elif pd.isna(row['lsoa_code']):
        return 'NO_POSTCODE_FOUND'
    elif str(row['lsoa_code']).startswith('S'):
        return 'SCOTLAND_NOT_IN_EW_RUC'
    elif str(row['lsoa_code']).startswith('N'):
        return 'NORTHERN_IRELAND_NOT_IN_EW_RUC'
    else:
        return 'NOT_MATCHED_IN_RUC'

output['RUC21NM'] = output.apply(categorize_missing, axis=1)
output = output[['garden.guid', 'epsg3857.x', 'epsg3857.y', 'RUC21NM']]

output.to_csv('garden_coordinates_rural_urban.csv', index=False)

print(f"\n✓ Done! Saved to garden_coordinates_rural_urban.csv")

# Create summary table
summary_lines = []
summary_lines.append("=" * 80)
summary_lines.append("SUMMARY")
summary_lines.append("=" * 80)

# Identify failed vs successful classifications
failed_labels = ['NO_POSTCODE_FOUND', 'SCOTLAND_NOT_IN_EW_RUC', 'NORTHERN_IRELAND_NOT_IN_EW_RUC', 'NOT_MATCHED_IN_RUC']
output['is_failed'] = output['RUC21NM'].isin(failed_labels)

total_gardens = len(output)
successful_gardens = (~output['is_failed']).sum()
failed_gardens = output['is_failed'].sum()

# Get counts by classification
counts = output['RUC21NM'].value_counts()

# Build summary table
summary_lines.append("")
summary_lines.append(f"{'Classification':<50} {'Count':>8} {'% of Success':>15}")
summary_lines.append(f"{'-'*50} {'-'*8} {'-'*15}")

for classification, count in counts.items():
    if classification not in failed_labels:
        pct = (count / successful_gardens * 100) if successful_gardens > 0 else 0
        summary_lines.append(f"{classification:<50} {count:>8} {pct:>14.1f}%")

summary_lines.append(f"{'-'*50} {'-'*8} {'-'*15}")
summary_lines.append(f"{'TOTAL SUCCESSFUL':<50} {successful_gardens:>8} {100.0:>14.1f}%")

if failed_gardens > 0:
    summary_lines.append("")
    summary_lines.append(f"{'FAILED CLASSIFICATIONS:':<50}")
    for classification, count in counts.items():
        if classification in failed_labels:
            pct = (count / total_gardens * 100)
            summary_lines.append(f"  {classification:<48} {count:>8} ({pct:.1f}% of all)")
    summary_lines.append(f"{'-'*50} {'-'*8}")
    summary_lines.append(f"{'TOTAL FAILED':<50} {failed_gardens:>8}")

summary_lines.append("")
summary_lines.append(f"{'GRAND TOTAL':<50} {total_gardens:>8}")
summary_lines.append("=" * 80)

# Print to console
print("\n" + "\n".join(summary_lines))

# Write to file
with open('rural_urban_summary.txt', 'w') as f:
    f.write("\n".join(summary_lines) + "\n")
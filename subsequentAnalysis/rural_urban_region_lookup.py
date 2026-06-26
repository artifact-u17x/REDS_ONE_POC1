# Declaration: this script produces output that is either non-critical (rural-urban distinction) or
# can easily be manually verified by map examination (region classification, and this verification
# was carried out. To save time, the script was therefore vibe-coded with less careful human checking
# than the analysis code.

import pandas as pd
import requests
import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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
# It was originally obtained from: https://www.data.gov.uk/dataset/1676bb60-bef1-43dc-8c23-641554b066bd/postcode-to-oa-2021-to-lsoa-to-msoa-to-lad-february-2025-best-fit-lookup-in-the-uk
lsoa_oa_lookup = pd.read_csv('geography_databases/PCD_OA21_LSOA21_MSOA21_LAD_FEB25_UK_LU.csv', low_memory=False)

# Load RUC file
# This file is not redistributed as I haven't checked rights to do so.
# It was originally obtained from:
# https://geoportal.statistics.gov.uk/datasets/ons::rural-urban-classification-2021-of-output-areas-in-ew/about
ruc = pd.read_csv('geography_databases/Rural_Urban_Classification_(2021)_of_Output_Areas_in_EW.csv')

# Load OA to Region lookup
# This file is not redistributed as I haven't checked rights to do so.
# It was originally obtained from:
# https://www.data.gov.uk/dataset/output-area-2021-to-parncp-to-lad-to-rgn-to-ctry-december-2024-best-fit-lookup-in-ew-v2
# Direct CSV download URL:
# https://open-geography-portalx-ons.hub.arcgis.com/api/download/v1/items/7507c0292db546ed83e4ba60f1115b1d/csv?layers=0
oa_rgn_lookup = pd.read_csv('geography_databases/OA21_to_RGN_CTRY_DEC24_EW.csv', low_memory=False)

print(f"  LSOA-OA lookup: {len(lsoa_oa_lookup):,} rows")
print(f"  RUC file: {len(ruc):,} Output Areas")
print(f"  OA-to-Region lookup: {len(oa_rgn_lookup):,} Output Areas")

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

print("\nJoining with Region data...")

# Identify the region name column (handles case variation in the downloaded file)
rgn_col_candidates = [c for c in oa_rgn_lookup.columns if c.upper().startswith('RGN') and c.upper().endswith('NM') and 'NMW' not in c.upper()]
if not rgn_col_candidates:
    print(f"\nERROR: Could not find region name column. Available columns: {oa_rgn_lookup.columns.tolist()}")
    exit(1)
rgn_col = rgn_col_candidates[0]
print(f"  Using region column: '{rgn_col}'")

# Identify the OA code column in the region lookup (handles case variation)
oa_col_candidates = [c for c in oa_rgn_lookup.columns if c.upper() == 'OA21CD']
if not oa_col_candidates:
    print(f"\nERROR: Could not find OA21CD column. Available columns: {oa_rgn_lookup.columns.tolist()}")
    exit(1)
oa_col_in_rgn = oa_col_candidates[0]

# Join with region lookup
gardens_with_rgn = gardens_with_ruc.merge(
    oa_rgn_lookup[[oa_col_in_rgn, rgn_col]],
    left_on='oa21cd',
    right_on=oa_col_in_rgn,
    how='left'
)
# Normalise the region column name
gardens_with_rgn = gardens_with_rgn.rename(columns={rgn_col: 'RGN24NM'})

print(f"  After Region join: {len(gardens_with_rgn)} rows")
print(f"  Rows with RGN24NM: {gardens_with_rgn['RGN24NM'].notna().sum()}")
print(f"  Sample Region values: {gardens_with_rgn['RGN24NM'].dropna().head(3).tolist()}")

# Group by garden and take the most common classification per LSOA
# (since one LSOA can have multiple OAs)
def get_most_common(series):
    """Get most common non-null value from a series"""
    series_clean = series.dropna()
    if len(series_clean) == 0:
        return None
    return series_clean.mode()[0] if len(series_clean.mode()) > 0 else series_clean.iloc[0]

gardens_final = gardens_with_rgn.groupby('garden.guid').agg({
    'epsg3857.x': 'first',
    'epsg3857.y': 'first',
    'lsoa_code': 'first',
    'RUC21NM': get_most_common,
    'RGN24NM': get_most_common,
}).reset_index()

# Replace None with informative labels based on LSOA code
def categorize_missing(row, col):
    if pd.notna(row[col]):
        return row[col]
    elif pd.isna(row['lsoa_code']):
        return 'NO_POSTCODE_FOUND'
    elif str(row['lsoa_code']).startswith('S'):
        return 'SCOTLAND_NOT_IN_EW_LOOKUP'
    elif str(row['lsoa_code']).startswith('N'):
        return 'NORTHERN_IRELAND_NOT_IN_EW_LOOKUP'
    else:
        return f'NOT_MATCHED_IN_{col}_LOOKUP'

gardens_final['RUC21NM'] = gardens_final.apply(lambda r: categorize_missing(r, 'RUC21NM'), axis=1)
gardens_final['RGN24NM'] = gardens_final.apply(lambda r: categorize_missing(r, 'RGN24NM'), axis=1)

# ── Assign region band ─────────────────────────────────────────────────────────

REGION_TO_BAND = {
    'South West':                       'South England',
    'South East':                       'South England',
    'London':                           'South England',
    'Wales':                            'Mid England and Wales',
    'West Midlands':                    'Mid England and Wales',
    'East Midlands':                    'Mid England and Wales',
    'East of England':                  'Mid England and Wales',
    'North West':                       'North England',
    'North East':                       'North England',
    'Yorkshire and The Humber':         'North England',
    'SCOTLAND_NOT_IN_EW_LOOKUP':        'Scotland and Northern Ireland',
    'NORTHERN_IRELAND_NOT_IN_EW_LOOKUP': 'Scotland and Northern Ireland',
}

def assign_band(rgn):
    return REGION_TO_BAND.get(rgn, None)

gardens_final['region_band'] = gardens_final['RGN24NM'].map(assign_band)

# ── Manual overrides for locations that failed automatic lookup ────────────────

MANUAL_BAND_OVERRIDES = {
    '4f340640-5011-422e-983f-4bb7f9b03deb': 'Scotland and Northern Ireland',  # Northern Ireland
    'f6a00867-44ec-47ec-8dc3-60ababbffd1e': 'South England',
    'ec7e9962-e65c-4f04-9049-a9321e380c55': 'South England',
    'c80e2b8f-7d5d-4c0b-809a-26209407e058': 'Scotland and Northern Ireland',  # Scotland
}

for guid, band in MANUAL_BAND_OVERRIDES.items():
    mask = gardens_final['garden.guid'] == guid
    if mask.sum() == 0:
        print(f"WARNING: Manual override GUID not found in data: {guid}")
    else:
        gardens_final.loc[mask, 'region_band'] = band

# Save
output = gardens_final[['garden.guid', 'epsg3857.x', 'epsg3857.y', 'RUC21NM', 'RGN24NM', 'region_band']]
output.to_csv('garden_coordinates_rural_urban_region.csv', index=False)

print(f"\n✓ Done! Saved to garden_coordinates_rural_urban_region.csv")

# ── Summary table ──────────────────────────────────────────────────────────────

failed_labels_ruc  = ['NO_POSTCODE_FOUND', 'SCOTLAND_NOT_IN_EW_LOOKUP', 'NORTHERN_IRELAND_NOT_IN_EW_LOOKUP', 'NOT_MATCHED_IN_RUC21NM_LOOKUP']
failed_labels_rgn  = ['NO_POSTCODE_FOUND', 'SCOTLAND_NOT_IN_EW_LOOKUP', 'NORTHERN_IRELAND_NOT_IN_EW_LOOKUP', 'NOT_MATCHED_IN_RGN24NM_LOOKUP']
failed_labels_band = []  # all locations should now have a band; none are treated as failures

def build_summary(df, col, failed_labels, label):
    lines = []
    lines.append("=" * 80)
    lines.append(f"SUMMARY: {label}")
    lines.append("=" * 80)

    df = df.copy()
    df['is_failed'] = df[col].isin(failed_labels)
    total = len(df)
    successful = (~df['is_failed']).sum()
    failed = df['is_failed'].sum()
    counts = df[col].value_counts()

    lines.append("")
    lines.append(f"{'Classification':<50} {'Count':>8} {'% of Success':>15}")
    lines.append(f"{'-'*50} {'-'*8} {'-'*15}")

    for classification, count in counts.items():
        if classification not in failed_labels:
            pct = (count / successful * 100) if successful > 0 else 0
            lines.append(f"{classification:<50} {count:>8} {pct:>14.1f}%")

    lines.append(f"{'-'*50} {'-'*8} {'-'*15}")
    lines.append(f"{'TOTAL SUCCESSFUL':<50} {successful:>8} {100.0:>14.1f}%")

    if failed > 0:
        lines.append("")
        lines.append(f"{'FAILED CLASSIFICATIONS:':<50}")
        for classification, count in counts.items():
            if classification in failed_labels:
                pct = (count / total * 100)
                lines.append(f"  {classification:<48} {count:>8} ({pct:.1f}% of all)")
        lines.append(f"{'-'*50} {'-'*8}")
        lines.append(f"{'TOTAL FAILED':<50} {failed:>8}")

    lines.append("")
    lines.append(f"{'GRAND TOTAL':<50} {total:>8}")
    lines.append("=" * 80)
    return lines

summary_lines = []
summary_lines += build_summary(output, 'RUC21NM', failed_labels_ruc, "Rural/Urban Classification (RUC 2021)")
summary_lines.append("")
summary_lines += build_summary(output, 'RGN24NM', failed_labels_rgn, "Region (December 2024)")
summary_lines.append("")
summary_lines += build_summary(output, 'region_band', failed_labels_band, "Region Band")

print("\n" + "\n".join(summary_lines))

with open('rural_urban_region_summary.txt', 'w') as f:
    f.write("\n".join(summary_lines) + "\n")

# ── Map: colour by region band, label unassigned locations ────────────────────

print("\nPlotting region band map...")

BAND_COLOURS = {
    'South England':                '#e06c4b',
    'Mid England and Wales':        '#4b9be0',
    'North England':                '#5ab56e',
    'Scotland and Northern Ireland': '#9b6bbf',
}
UNASSIGNED_COLOUR = 'red'

lon_lat = transformer.transform(output['epsg3857.x'].values, output['epsg3857.y'].values)
output = output.copy()
output['_lon'] = lon_lat[0]
output['_lat'] = lon_lat[1]

unassigned = output[output['region_band'].isna()]
assigned   = output[output['region_band'].notna()]

# Print unassigned identities to console
if len(unassigned) > 0:
    print(f"\nLocations with no region band assigned ({len(unassigned)}):")
    for _, row in unassigned.iterrows():
        print(f"  {row['garden.guid']}  (RGN24NM: {row['RGN24NM']})")
else:
    print("\nAll locations successfully assigned to a region band.")

fig, ax = plt.subplots(figsize=(7, 10))

for band, colour in BAND_COLOURS.items():
    grp = assigned[assigned['region_band'] == band]
    ax.scatter(grp['_lon'], grp['_lat'], s=10, alpha=0.7, linewidths=0,
               color=colour, label=band, zorder=2)

if len(unassigned) > 0:
    ax.scatter(unassigned['_lon'], unassigned['_lat'], s=30, color=UNASSIGNED_COLOUR,
               marker='x', linewidths=1.2, label='Unassigned', zorder=3)
    for _, row in unassigned.iterrows():
        ax.annotate(str(row['garden.guid']),
                    xy=(row['_lon'], row['_lat']),
                    xytext=(4, 4), textcoords='offset points',
                    fontsize=5, color=UNASSIGNED_COLOUR, zorder=4)

ax.set_aspect('equal')
ax.set_title('Gardens by region band', fontsize=11)
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
legend_handles = [mpatches.Patch(color=c, label=b) for b, c in BAND_COLOURS.items()]
if len(unassigned) > 0:
    legend_handles.append(plt.Line2D([0], [0], marker='x', color='w',
                                     markeredgecolor=UNASSIGNED_COLOUR,
                                     markersize=7, label='Unassigned'))
ax.legend(handles=legend_handles, fontsize=7, loc='lower left', framealpha=0.8)
plt.tight_layout()
plt.show()
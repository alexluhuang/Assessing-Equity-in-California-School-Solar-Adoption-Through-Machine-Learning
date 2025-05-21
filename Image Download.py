import os
import requests
import pandas as pd
import math
from dotenv import load_dotenv
import matplotlib.pyplot as plt
from PIL import Image
import geopandas as gpd

load_dotenv()
API_KEY = os.getenv("API_KEY")
BASE_URL = os.getenv("BASE_URL")

# Constants
EARTH_CIRCUMFERENCE_METERS = 40075017
TILE_SIZE_PIXELS = 400
ZOOM_LEVEL = 20

# Define the path to your Excel file
excel_file_path = r"C:\Users\huang\Documents\Deepsolar (Broken)\solar_models\Test Data\Midsized Suburb\Midsized Suburb Test Data.xlsx"

# Read Excel file
df = pd.read_excel(excel_file_path)

# Update these to match the exact column names
latitude_column = next((col for col in df.columns if str(col).strip().lower() == 'latitude'), None)
longitude_column = next((col for col in df.columns if str(col).strip().lower() == 'longitude'), None)
fed_id_column = next((col for col in df.columns if str(col).strip().lower() == 'fed id'), None)
locale_column = next((col for col in df.columns if str(col).strip().lower() == 'locale'), None)

# Define the target directory
target_directory = r'C:\Users\huang\Documents\Deepsolar (Broken)\solar_models\data_folder'
os.makedirs(target_directory, exist_ok=True)


def calculate_meters_per_pixel(latitude):
    meters_per_pixel_at_equator = EARTH_CIRCUMFERENCE_METERS / (2 ** ZOOM_LEVEL * TILE_SIZE_PIXELS)
    return meters_per_pixel_at_equator / math.cos(math.radians(latitude))


def calculate_offsets(lat, meters_per_pixel):
    lat_offset = TILE_SIZE_PIXELS * meters_per_pixel / 111111
    lon_offset = TILE_SIZE_PIXELS * meters_per_pixel / (111111 * abs(math.cos(math.radians(lat))))
    return lat_offset, lon_offset


def fetch_satellite_image(center, zoom, size, filename):
    params = {
        'center': center.strip(),  # Stripping spaces to ensure correctness
        'zoom': zoom,
        'size': size,
        'maptype': 'satellite',
        'key': API_KEY,
    }
    response = requests.get(BASE_URL, params=params)
    if response.status_code == 200:
        with open(filename, 'wb') as file:
            file.write(response.content)
        print(f"Image saved: {filename}")
    else:
        print(f"Error fetching image for '{center}': {response.status_code} - {response.text}")


def get_grid_size(locale):
    locale = locale.strip()
    if locale == "11 - City, Large":
        return 5
    elif locale == "12 - City, Midsize":
        return 7
    elif locale == "13 - City, Small":
        return 7
    elif locale == "21 - Suburban, Large":
        return 7
    elif locale == "22 - Suburban, Midsize":
        return 3
    else:
        return 7


def display_grid_images(fed_id, grid_size):
    """Arrange and display the downloaded images in a grid using matplotlib without borders."""
    fig, ax = plt.subplots(grid_size, grid_size, figsize=(10, 10))

    # Flip rows correctly: top row should be index (1,1)
    for i in range(grid_size):
        for j in range(grid_size):
            # Reverse the row indexing for display (bottom row on top)
            filename = os.path.join(target_directory, f"{fed_id}_{grid_size - i}_{j + 1}.jpg")
            if os.path.exists(filename):
                img = Image.open(filename)
                ax[i, j].imshow(img)
                ax[i, j].axis('off')
            else:
                ax[i, j].set_facecolor('gray')
                ax[i, j].text(0.5, 0.5, "No Image", ha='center', va='center')

    # Invert the order of the rows to fix the visualization
    ax = ax[::-1]

    # Removing all padding and spacing between images
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.margins(0, 0)
    plt.show()


# Main Processing Loop (Consistent File Naming and Display)
for idx, row in df.iterrows():
    lat = row[latitude_column]
    lon = row[longitude_column]
    locale = row[locale_column]
    grid_size = get_grid_size(locale)

    # If fed_id_column does not exist or is empty, create an ID from row index
    fed_id = row[fed_id_column] if fed_id_column in row and not pd.isna(row[fed_id_column]) else f"ID_{idx}"

    # Adjust meters per pixel for the specific latitude
    meters_per_pixel = calculate_meters_per_pixel(lat)

    # Calculate offsets for latitude and longitude
    lat_offset, lon_offset = calculate_offsets(lat, meters_per_pixel)

    # Generate images in the grid with corrected file naming
    half_size = grid_size // 2
    for i in range(-half_size, half_size + 1):
        for j in range(-half_size, half_size + 1):
            adj_lat = lat + (i * lat_offset)
            adj_lon = lon + (j * lon_offset)
            # Flip the row ordering for proper top-left indexing during saving
            filename = os.path.join(target_directory, f"{fed_id}_{i + half_size + 1}_{j + half_size + 1}.jpg")
            location_str = f"{adj_lat},{adj_lon}"
            print(f"Fetching image for {fed_id}: Lat = {adj_lat}, Lon = {adj_lon}")
            fetch_satellite_image(location_str, ZOOM_LEVEL, '400x400', filename)

    # Display the corrected grid of images for verification
    display_grid_images(fed_id, grid_size)

import os  # We'll handle file paths, directories, and environment variables with this module.
import requests  # Allows us to send HTTP requests, which we'll use to fetch satellite images.
import pandas as pd  # Handy for reading and manipulating tabular data, like Excel files.
import math  # Used for various math operations, including trigonometric functions.
from dotenv import load_dotenv  # Loads key-value pairs from a .env file into environment variables.
import geopandas as gpd  # Extends pandas to work with geospatial data (GIS shapes).
from shapely.geometry import Polygon, Point  # Geometry objects we need for bounding boxes and location points.
import pyproj  # A library for handling spatial references and coordinate transformations.
import hashlib
import hmac
import base64
import urllib.parse as urlparse

# We tell pyproj where to find its projection files.
os.environ["PROJ_LIB"] = r"C:\Users\huang\anaconda3\envs\Deepsolar\Library\share\proj"

# Just to confirm that pyproj is indeed seeing the correct data directory:
print("PyProj data directory:", pyproj.datadir.get_data_dir())

# 1. Load your .gdb (geodatabase) which contains the school polygons
gdb_path = r"C:\Users\huang\Documents\Deepsolar (Broken)\solar_models\CSCD_2021.gdb"
layer_name = "Schools_Current_Stacked"

# Read the specified layer from the geodatabase as a GeoDataFrame
school_parcels_gdf = gpd.read_file(gdb_path, layer=layer_name)
print("Original parcels CRS:", school_parcels_gdf.crs)

# If GeoPandas didn't recognize the coordinate reference system, set it manually:
if school_parcels_gdf.crs is None:
    # We assume EPSG:3310 (NAD83 / California (Teale) Albers) for these polygons
    school_parcels_gdf.crs = "EPSG:3310"

# Reproject polygons from EPSG:3310 to EPSG:4326 (lat/lon)
school_parcels_gdf_4326 = school_parcels_gdf.to_crs(epsg=4326)
print("Reprojected parcels CRS:", school_parcels_gdf_4326.crs)

# Load environment variables from .env (should have API_KEY, BASE_URL, etc.)
load_dotenv()
API_KEY = os.getenv("API_KEY")
BASE_URL = os.getenv("BASE_URL")
SIGNATURE = os.getenv("SIGNATURE")

def sign_url(input_url=None, secret=None):
    """ Sign a request URL with a URL signing secret.
      Usage:
      from urlsigner import sign_url
      signed_url = sign_url(input_url=my_url, secret=SECRET)
      Args:
      input_url - The URL to sign
      secret    - Your URL signing secret
      Returns:
      The signed request URL
  """

    if not input_url or not secret:
        raise Exception("Both input_url and secret are required")

    url = urlparse.urlparse(input_url)

    # We only need to sign the path+query part of the string
    url_to_sign = url.path + "?" + url.query

    # Decode the private key into its binary format
    # We need to decode the URL-encoded private key
    decoded_key = base64.urlsafe_b64decode(secret)

    # Create a signature using the private key and the URL-encoded
    # string using HMAC SHA1. This signature will be binary.
    signature = hmac.new(decoded_key, str.encode(url_to_sign), hashlib.sha1)

    # Encode the binary signature into base64 for use within a URL
    encoded_signature = base64.urlsafe_b64encode(signature.digest())

    original_url = url.scheme + "://" + url.netloc + url.path + "?" + url.query

    # Return signed URL
    return original_url + "&signature=" + encoded_signature.decode()

# Constants we’ll use when requesting static map tiles
EARTH_CIRCUMFERENCE_METERS = 40075017  # Earth's approximate circumference
TILE_SIZE_PIXELS = 400                # We'll request 400x400 tiles
ZOOM_LEVEL = 20                       # High zoom level to see details

# Path to the Excel file that contains the school's lat/lon info
excel_file_path = r"C:\Users\huang\Documents\Deepsolar (Broken)\Test Data.xlsx"
df = pd.read_excel(excel_file_path)

# Identify relevant columns by name, ignoring capitalization
latitude_column = next((c for c in df.columns if c.strip().lower() == 'latitude'), None)
longitude_column = next((c for c in df.columns if c.strip().lower() == 'longitude'), None)
cds_code_column = next((c for c in df.columns if c.strip().lower() == 'cds code'), None)
locale_column = next((c for c in df.columns if c.strip().lower() == 'locale'), None)

# Directory where we’ll store downloaded satellite images
target_directory = r"C:\Users\huang\Documents\Deepsolar (Broken)\solar_models\data_folder\Test Batch\Division then None"
os.makedirs(target_directory, exist_ok=True)  # Create it if it doesn't exist

# Helper Functions
def calculate_meters_per_pixel(latitude):
    """
    Calculates how many meters correspond to one pixel at zoom=20 for a given latitude.
    The closer to the poles, the larger the 'cos(lat)' effect.
    """
    # Base meters/pixel at the equator for zoom 20 with our tile size
    meters_per_pixel_at_equator = EARTH_CIRCUMFERENCE_METERS / (2 ** ZOOM_LEVEL * TILE_SIZE_PIXELS)
    return meters_per_pixel_at_equator / math.cos(math.radians(latitude))

def fetch_satellite_image(center, zoom, size, filename):
    """
    Downloads one tile from the Google Maps Static API, saving it to 'filename'.
    Incorporates URL signing to ensure the requests are properly authenticated.
    """
    # Step 1: Build the base (unsigned) URL string
    #         (Do not include 'signature=' in this step.)
    query_params = (
        f"center={center.strip()}"
        f"&zoom={zoom}"
        f"&size={size}"
        f"&maptype=satellite"
        f"&key={API_KEY}"
    )

    # e.g., https://maps.googleapis.com/maps/api/staticmap?center=...&zoom=...&...
    unsigned_url = BASE_URL + "?" + query_params

    # Step 2: Sign the URL using your sign_url function
    signed_url = sign_url(unsigned_url, SIGNATURE)  # <--- your secret key!

    # Step 3: Fetch the signed URL
    response = requests.get(signed_url)
    if response.status_code == 200:
        with open(filename, 'wb') as file:
            file.write(response.content)
        print(f"Image saved: {filename}")
    else:
        print(f"Error fetching image for '{center}': {response.status_code} - {response.text}")

def create_tile_polygon(center_lat, center_lon, meters_per_pixel):
    """
    Builds the bounding box polygon for a 400x400 tile at zoom=20
    by calculating how many degrees each pixel covers.
    """
    # Convert meters to degrees. 1 degree ~111111 meters in lat dimension.
    deg_per_pixel_lat = meters_per_pixel / 111111
    deg_per_pixel_lon = meters_per_pixel / (111111 * abs(math.cos(math.radians(center_lat))))

    # Half the tile means half its dimension in degrees
    half_tile_offset_lat = (TILE_SIZE_PIXELS / 2) * deg_per_pixel_lat
    half_tile_offset_lon = (TILE_SIZE_PIXELS / 2) * deg_per_pixel_lon

    # Compute min/max lat/lon bounds of the tile
    lat_min = center_lat - half_tile_offset_lat
    lat_max = center_lat + half_tile_offset_lat
    lon_min = center_lon - half_tile_offset_lon
    lon_max = center_lon + half_tile_offset_lon

    # Return a shapely Polygon that outlines the tile
    coords = [
        (lon_min, lat_min),
        (lon_min, lat_max),
        (lon_max, lat_max),
        (lon_max, lat_min),
        (lon_min, lat_min)
    ]
    return Polygon(coords)

def get_polygon_by_cds_code(cds_code, gdf_4326):
    """
    Looks up a polygon in the GeoDataFrame by matching the 'CDSCode' field.
    Returns a single geometry (via unary_union) if found, or None if it doesn't exist.
    """
    matching = gdf_4326[gdf_4326["CDSCode"] == cds_code]
    if matching.empty:
        return None
    return matching.unary_union

def get_school_parcel_geometry(lat, lon, gdf_4326):
    """
    Uses the point-in-polygon check to find which polygon(s) contain a given lat/lon.
    Returns the combined polygon geometry (unary_union) or None if no match is found.
    """
    point_geom = Point(lon, lat)  # remember: x=lon, y=lat
    matching_rows = gdf_4326[gdf_4326.geometry.contains(point_geom)]
    if matching_rows.empty:
        return None
    return matching_rows.unary_union

# Main Processing Loop
for idx, row in df.iterrows():
    # Grab latitude & longitude from the Excel row
    lat_from_excel = row[latitude_column]
    lon_from_excel = row[longitude_column]

    # If there's a valid CDS Code, use it; otherwise just label it with index
    if cds_code_column in df.columns and not pd.isna(row[cds_code_column]):
        cds_code = str(row[cds_code_column])
    else:
        cds_code = f"CDS_{idx}"

    # First, see if we have a polygon in the GDB that matches this CDS code
    polygon_by_code = get_polygon_by_cds_code(cds_code, school_parcels_gdf_4326)
    if polygon_by_code is not None:
        print(f"CDS code {cds_code} found in gdb by matching code. Using bounding-box approach.")
        school_polygon = polygon_by_code
    else:
        # If not, we try a point-in-polygon approach based on lat/lon
        school_polygon = get_school_parcel_geometry(lat_from_excel, lon_from_excel, school_parcels_gdf_4326)
        if school_polygon is None:
            print(f"No polygon found for {cds_code} by code or point. We'll do a 7x7 fallback.")
            school_polygon = None

    if school_polygon is not None:
        # --------------------------------------------
        # BOUNDING BOX approach
        # --------------------------------------------
        minx, miny, maxx, maxy = school_polygon.bounds

        # We'll take the vertical midpoint to get an approximate latitude
        # for the meters-per-pixel calculation
        center_lat_for_calc = (miny + maxy) / 2
        meters_per_pixel = calculate_meters_per_pixel(center_lat_for_calc)

        # Convert from meters to degrees for lat and lon
        deg_per_pixel_lat = meters_per_pixel / 111111
        deg_per_pixel_lon = meters_per_pixel / (111111) * abs(math.cos(math.radians(center_lat_for_calc)))

        # Total width/height of one tile in degrees
        tile_width_deg = deg_per_pixel_lon * TILE_SIZE_PIXELS
        tile_height_deg = deg_per_pixel_lat * TILE_SIZE_PIXELS

        # Expand the bounding box slightly so we don't cut off edges
        left_bound = minx - (tile_width_deg / 2)
        right_bound = maxx + (tile_width_deg / 2)
        bottom_bound = miny - (tile_height_deg / 2)
        top_bound = maxy + (tile_height_deg / 2)

        print(f"DEBUG: {cds_code} bounding box approach - "
              f"({left_bound}, {bottom_bound}) to ({right_bound}, {top_bound})")

        # We'll move in rows (y-direction) and columns (x-direction) across the bounding box
        row_index = 0
        current_y = bottom_bound

        while current_y < top_bound:
            col_index = 0
            current_x = left_bound
            while current_x < right_bound:
                # Find center of this tile in degrees
                center_lon_tile = current_x + (tile_width_deg / 2)
                center_lat_tile = current_y + (tile_height_deg / 2)

                # Build a polygon for this tile and see if it intersects the school polygon
                tile_poly = create_tile_polygon(center_lat_tile, center_lon_tile, meters_per_pixel)
                if tile_poly.intersects(school_polygon):
                    # Intersects, so let's fetch it
                    filename = os.path.join(
                        target_directory,
                        f"{cds_code}_r{row_index}_c{col_index}.jpg"
                    )
                    center_str = f"{center_lat_tile},{center_lon_tile}"
                    print(f"Fetching tile row={row_index}, col={col_index} for {cds_code}: {center_str}")
                    fetch_satellite_image(center_str, ZOOM_LEVEL, '400x400', filename)
                else:
                    print(f"SKIP tile row={row_index}, col={col_index} at "
                          f"{center_lat_tile},{center_lon_tile} (no intersection)")

                col_index += 1
                current_x += tile_width_deg  # Move one tile width in the x-direction
            row_index += 1
            current_y += tile_height_deg  # Move one tile height in the y-direction

    else:
        # If no polygon was found, use a 7×7 fallback grid.
        print(f"CDS code {cds_code} not in GDB polygon. Using 7x7 fallback around lat={lat_from_excel}, lon={lon_from_excel}.")

        grid_size = 7
        half_size = grid_size // 2

        meters_per_pixel = calculate_meters_per_pixel(lat_from_excel)
        lat_offset = TILE_SIZE_PIXELS * meters_per_pixel / 111111
        lon_offset = TILE_SIZE_PIXELS * meters_per_pixel / (
            111111) #* abs(math.cos(math.radians(lat_from_excel))

        print(f"DEBUG: 7x7 fallback approach for {cds_code} - offsets: lat={lat_offset}, lon={lon_offset}")

        # We'll loop in both directions from -half_size to +half_size
        for i in range(-half_size, half_size + 1):
            row_index = i + half_size
            for j in range(-half_size, half_size + 1):
                col_index = j + half_size

                # Adjust lat/lon for the current tile in the fallback grid
                adj_lat = lat_from_excel + (i * lat_offset)
                adj_lon = lon_from_excel + (j * lon_offset)

                # Create a filename for the image
                filename = os.path.join(
                    target_directory,
                    f"{cds_code}_r{row_index}_c{col_index}.jpg"
                )
                center_str = f"{adj_lat},{adj_lon}"

                # Fetch the tile
                print(f"Fetching tile row={row_index}, col={col_index} for {cds_code}: {center_str}")
                fetch_satellite_image(center_str, ZOOM_LEVEL, '400x400', filename)

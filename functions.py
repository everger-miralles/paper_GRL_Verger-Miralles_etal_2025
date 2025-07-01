# functions

import pandas as pd
import numpy  as np
import xarray as xr
from scipy.signal import butter, filtfilt

def process_and_concat_duacs(data_array, lon_min, lon_max, lat_min, lat_max, df_duacs=None):
  """
  Processes a single data array and concatenates it to an existing dataframe (if provided).

  """
  # Spatial subsetting using vectorized operations
  data_array = data_array.where((data_array['longitude'].compute() < lon_max) & (data_array['longitude'].compute() > lon_min) & 
                              (data_array['latitude'].compute() < lat_max) & (data_array['latitude'].compute() > lat_min), drop=True)

  # Convert to dataframe and select desired variables
  df = data_array.to_dataframe().reset_index()
  df = df[['time', 'longitude', 'latitude', 'adt']]

  # Concatenate if df_duacs is provided
  if df_duacs is None:
      return df
  else:
      return pd.concat([df_duacs, df], ignore_index=True)
  
# def process_and_concat_swot(data_array, lon_min, lon_max, lat_min, lat_max, df_swot=None):
#   """
#   Processes a single data array and concatenates it to an existing dataframe (if provided).

#   """
#   # Remove num_nadir dimension
#   data_array = data_array.drop_dims('num_nadir')
#   # Spatial subsetting using vectorized operations
#   data_array = data_array.where((data_array['longitude'].compute() < lon_max) & (data_array['longitude'].compute() > lon_min) & 
#                               (data_array['latitude'].compute() < lat_max) & (data_array['latitude'].compute() > lat_min), drop=True)

#   # Convert to dataframe and select desired variables
#   df = data_array.to_dataframe().reset_index()
#   df = df[['time', 'longitude', 'latitude', 'ssha_noiseless', 'mdt']]

#   # Concatenate if df_duacs is provided
#   if df_swot is None:
#       return df
#   else:
#       return pd.concat([df_swot, df], ignore_index=True)
  
def process_and_concat_swot(data_array, lon_min, lon_max, lat_min, lat_max, df_swot=None):
  """
  Processes a single data array and concatenates it to an existing dataframe (if provided).

  """
  # Remove num_nadir dimension
#   data_array = data_array.drop_dims('num_nadir')
  # Spatial subsetting using vectorized operations
  data_array = data_array.where((data_array['longitude'].compute() < lon_max) & (data_array['longitude'].compute() > lon_min) & 
                              (data_array['latitude'].compute() < lat_max) & (data_array['latitude'].compute() > lat_min), drop=True)

  # Convert to dataframe and select desired variables
  df = data_array.to_dataframe().reset_index()
  df = df[['time', 'longitude', 'latitude', 'ssha_filtered', 'mdt', 'ugos_filtered', 'vgos_filtered']]

  # Concatenate if df_duacs is provided
  if df_swot is None:
      return df
  else:
      return pd.concat([df_swot, df], ignore_index=True)
  
# Define a low-pass Butterworth filter
def low_pass_filter(data, cutoff, fs, order=4):
    nyquist = 0.5 * fs
    normal_cutoff = cutoff / nyquist  # cutoff = 1 / [km]
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y

    # Digital filter critical frequencies must be 0 < Wn < 1


def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great-circle distance between two points on the Earth's surface.
    
    Parameters:
    lat1, lon1 -- latitude and longitude of point 1 (in decimal degrees)
    lat2, lon2 -- latitude and longitude of point 2 (in decimal degrees)
    
    Returns:
    distance -- distance between the two points (in kilometers)
    """
    R = 6371.0  # Earth’s radius in kilometers
    
    # Convert latitude and longitude from degrees to radians
    lat1, lon1, lat2, lon2 = map(np.radians, [lat1, lon1, lat2, lon2])
    
    # Haversine formula
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = np.sin(dlat / 2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0)**2
    c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
    
    distance = R * c
    return distance


def compute_rossby_number_2d(zeta, latitude):
    """
    Compute the Rossby number Ro = ζ / f in 2D.
    
    Parameters:
    zeta: 2D array of relative vorticity (s⁻¹)
    latitude: Single latitude value (degrees)

    Returns:
    Ro: 2D array of Rossby number
    """
    omega = 7.2921e-5  # Earth's rotation rate (rad/s)
    f = 2 * omega * np.sin(np.radians(latitude))  # Compute Coriolis parameter

    if f == 0:
        f = np.nan  # Avoid division by zero at the equator

    Ro = zeta / f
    return Ro


def compute_vorticity_2d(u_cross, y):
    """
    Compute vertical relative vorticity (ζ) in 2D using finite differences.
    
    Parameters:
    u_cross: 2D array (depth x section distance) of cross-section velocity (m/s)
    y: 1D array of distances along the section (m)

    Returns:
    zeta: 2D array (depth x section distance) of relative vorticity (s⁻¹)
    """
    zeta = -np.gradient(u_cross, y, axis=1)  # Compute -du/dy along the section
    return zeta

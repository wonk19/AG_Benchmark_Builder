"""
GOCI-II Data Reader Module

Provides functions to read and handle GOCI-II Level-2 NetCDF files.
"""

import netCDF4
import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional, Tuple


@dataclass
class GOCI2Data:
    """Container for GOCI-II data"""
    rhoc: Dict[str, np.ndarray]
    rrs: Dict[str, np.ndarray]
    latitude: np.ndarray
    longitude: np.ndarray
    flag: np.ndarray
    dimensions: Tuple[int, int]
    
    def get_band(self, band_type: str, wavelength: int) -> Optional[np.ndarray]:
        """
        Get a specific band by type and wavelength.
        
        Parameters:
        -----------
        band_type : str
            'rhoc' or 'rrs'
        wavelength : int
            Wavelength in nm (e.g., 443, 490, 555, etc.)
        
        Returns:
        --------
        np.ndarray or None
            Band data or None if not found
        """
        if band_type.lower() == 'rhoc':
            key = f'RhoC_{wavelength}'
            return self.rhoc.get(key)
        elif band_type.lower() == 'rrs':
            key = f'Rrs_{wavelength}'
            return self.rrs.get(key)
        return None
    
    def get_rgb(self, band_type: str = 'rhoc', 
                red: int = 660, green: int = 555, blue: int = 443) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get RGB bands for visualization.
        
        Parameters:
        -----------
        band_type : str
            'rhoc' or 'rrs'
        red, green, blue : int
            Wavelengths for RGB channels
        
        Returns:
        --------
        Tuple[np.ndarray, np.ndarray, np.ndarray]
            Red, green, blue bands
        """
        r = self.get_band(band_type, red)
        g = self.get_band(band_type, green)
        b = self.get_band(band_type, blue)
        return r, g, b


def read_goci2(filepath: str) -> GOCI2Data:
    """
    Read GOCI-II Level-2 NetCDF file.
    
    Parameters:
    -----------
    filepath : str
        Path to GOCI-II NetCDF file
    
    Returns:
    --------
    GOCI2Data
        Data container with all bands and geolocation
    """
    ds = netCDF4.Dataset(filepath, 'r')
    
    dimensions = (ds.dimensions['number_of_lines'].size, 
                  ds.dimensions['pixels_per_line'].size)
    
    rhoc_group = ds.groups['geophysical_data'].groups['RhoC']
    rrs_group = ds.groups['geophysical_data'].groups['Rrs']
    nav_group = ds.groups['navigation_data']
    
    rhoc = {}
    for var_name in rhoc_group.variables.keys():
        rhoc[var_name] = rhoc_group.variables[var_name][:]
    
    rrs = {}
    for var_name in rrs_group.variables.keys():
        rrs[var_name] = rrs_group.variables[var_name][:]
    
    latitude = nav_group.variables['latitude'][:]
    longitude = nav_group.variables['longitude'][:]
    
    flag = ds.groups['geophysical_data'].variables['flag'][:]
    
    ds.close()
    
    return GOCI2Data(
        rhoc=rhoc,
        rrs=rrs,
        latitude=latitude,
        longitude=longitude,
        flag=flag,
        dimensions=dimensions
    )


# Flag handling functions
GOCI2_FLAG_MEANINGS = {
    0: 'COASTLINE',
    1: 'LAND',
    2: 'CLOUD',
    3: 'HIGH_GLINT',
    4: 'CLOUD_SHADOW',
    5: 'NEGATIVE_RRS',
    6: 'TURBID_WATER',
    7: 'COCCOLITHOPHORE',
    16: 'AC_FAIL',
}


def get_flag_bit(flag_array: np.ndarray, bit: int) -> np.ndarray:
    """Extract a specific bit from the flag array"""
    return ((flag_array >> bit) & 1).astype(np.uint8)


def get_flag_by_name(flag_array: np.ndarray, flag_name: str) -> np.ndarray:
    """
    Extract flag by name.
    
    Parameters:
    -----------
    flag_array : np.ndarray
        2D array of flag values
    flag_name : str
        Name of the flag (e.g., 'LAND', 'CLOUD')
    
    Returns:
    --------
    np.ndarray
        Binary array (1 where flag is set, 0 otherwise)
    """
    flag_name = flag_name.upper()
    
    for bit, name in GOCI2_FLAG_MEANINGS.items():
        if name == flag_name:
            return get_flag_bit(flag_array, bit)
    
    raise ValueError(f"Unknown flag name: {flag_name}")


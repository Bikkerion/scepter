"""
tleforger.py

This module generates artificial TLEs (Two-Line Element sets) for satellites.

Author: boris.sorokin <mralin@protonmail.com>
Date: 16-01-2025
"""
import numpy as np
from typing import List
from astropy.time import Time
from astropy.constants import GM_earth, R_earth

def _compute_tle_checksum(tle_line: str) -> int:
    """
    Compute TLE checksum for a single line.
    
    Sums all digits (ignoring '.'), and for each '-', adds +1, then takes mod 10.
    
    Parameters:
    tle_line (str): The TLE line for which to compute the checksum.
    
    Returns:
    int: The computed checksum.
    """
    return sum(int(ch) if ch.isdigit() else 1 if ch == '-' else 0 for ch in tle_line) % 10

def forge_tle_single(
    sat_name: str,
    altitude_m: float,
    eccentricity: float,
    inclination_deg: float,
    raan_deg: float,
    argp_deg: float,
    anomaly_deg: float,
    start_time: Time
) -> str:
    """
    Forge a minimal TLE from a simplified set of inputs.
    
    Returns a 3-line string: line0 (name), line1, line2
    
    Parameters:
    sat_name (str): Satellite name.
    altitude_m (float): Altitude above Earth's surface in meters.
    eccentricity (float): Orbital eccentricity.
    inclination_deg (float): Inclination in degrees.
    raan_deg (float): Right Ascension of the Ascending Node in degrees.
    argp_deg (float): Argument of Perigee in degrees.
    anomaly_deg (float): Mean Anomaly in degrees.
    start_time (Time): Start time as an astropy Time object.
    
    Returns:
    str: Three-line TLE string.
    """
    def format_leading_zero(value: float) -> str:
        """
        Format the second derivative of mean motion for TLE line 1 (columns 45–52).

        Returns an 8-character string:
        [sign or space][5-digit mantissa][sign exponent][1-digit exponent]

        Examples:
        0.0 -> " 00000-0"
        1.2345e-5 -> " 12345-5"
        -2.3e-3 -> "-02300-3"
        """
        # If it's effectively zero, return " 00000-0"
        if abs(value) < 1e-99:  # treat as zero
            return " 00000-0"
        
        sign_mantissa = '-' if value < 0 else ' '
        vabs = abs(value)

        # Step 1: Convert to scientific notation. We want 1 <= vabs < 10 if possible.
        exponent = 0
        if vabs >= 1.0:
            while vabs >= 10.0 and exponent < 9:
                vabs /= 10.0
                exponent += 1
        else:  # vabs < 1
            while vabs < 1.0 and exponent > -9:
                vabs *= 10.0
                exponent -= 1

        # Step 2: We only keep 5 digits in the mantissa
        # e.g. 1.2345678 -> 12346 -> "12346"
        mantissa_int = round(vabs * 10000)
        if mantissa_int >= 100000:
            # e.g. 9.9995 -> 100000
            # push exponent up by 1
            mantissa_int //= 10
            exponent += 1
            if exponent > 9:
                # saturate if needed, though extremely rare
                exponent = 9

        mantissa_str = f"{mantissa_int:05d}"  # exactly 5 digits

        # Step 3: sign and digit for exponent
        sign_exp = '-' if exponent < 0 else '+'
        exponent_str = str(abs(exponent))

        # Combine into 8 chars
        return f"{sign_mantissa}{mantissa_str}{sign_exp}{exponent_str}"
    
    # 1. Basic TLE identifiers
    sat_number = 0           # Fake NORAD ID
    classification = "U"
    int_desg = "25001A"      # Fake International Designator
    
    # 2. Epoch calculation
    year = start_time.datetime.year
    year_short = year % 100
    start_of_year = Time(start_time.datetime.replace(month=1, day=1, hour=0, minute=0, second=0, microsecond=0))
    day_of_year = (start_time - start_of_year).to_value('day') + 1  # Day count starts at 1
    epoch_str = f"{year_short:02d}{day_of_year:012.8f}"  # Leading space for epoch
    
    # 3. Mean motion calculation (circular orbit assumption)
    a_m = R_earth.value + altitude_m
    mean_motion_rad_s = np.sqrt(GM_earth.value / a_m**3)
    mean_motion_rev_day = mean_motion_rad_s * (86400.0 / (2.0 * np.pi))
    
    # 4. Zero placeholders
    mm_dot = 0.0
    mm_dot_string = f"{mm_dot:+.8f}"[:1]+f"{mm_dot:+.8f}"[2:]
    mm_ddot = 0.0
    bstar = 0.00
    ephemeris_type = 0
    element_set_number = 1
    rev_number = 1  # Dummy revolution number
    
    # Build line 1 (no checksum yet)
    line1 = (f"1 {sat_number:05d}{classification} {int_desg:8} {epoch_str:14} {mm_dot_string} {format_leading_zero(mm_ddot)} {format_leading_zero(bstar)} {ephemeris_type} {element_set_number:4d}")
    # f"1 {sat_number:05d}{classification} {int_desg:8} {epoch_str:14}{mm_dot:+11.8f}{mm_ddot:+11.8f}{bstar_str:>9}{ephemeris_type:1d}{element_set_number:4d}"
    
    # Ensure line1 is 68 characters before adding checksum
    if len(line1) != 68:
        raise ValueError("Line 1 is not 68 characters long before adding checksum.")
    
    # Build line 2 (no checksum yet)
    ecc_str = f"{eccentricity:.7f}"[2:]  # Remove '0.' to get 7 digits
    line2 = (
        "2 {:05d} {:8.4f} {:8.4f} {:7} {:8.4f} {:8.4f} {:11.8f}{:05d}"
        .format(
            sat_number, inclination_deg, raan_deg, ecc_str,
            argp_deg, anomaly_deg, mean_motion_rev_day, rev_number
        )
    )
    
    # Ensure line2 is 68 characters before adding checksum
    if len(line2) != 68:
        print(f"'{line2}'")
        print(f"Length is {len(line2)}")
        raise ValueError("Line 2 is not 68 characters long before adding checksum.")
    
    # 5. Compute and append checksums
    line1 += str(_compute_tle_checksum(line1))
    line2 += str(_compute_tle_checksum(line2))
    
    # 6. Construct the full TLE string
    tle = f"{sat_name}\n{line1}\n{line2}"
    
    return tle

def forge_tle_belt(
    num_sats_per_plane: int,
    plane_count: int,
    altitude_km: float,
    eccentricity: float,
    inclination_deg: float,
    argp_deg: float,
    start_time: Time
) -> List[str]:
    """
    Generate a list of TLEs for a belt of satellites.
    
    Parameters:
    num_sats_per_plane (int): Number of satellites per plane.
    plane_count (int): Number of planes.
    altitude_km (float): Altitude of the belt in kilometers.
    eccentricity (float): Eccentricity of the orbits.
    inclination_deg (float): Inclination of the orbits in degrees.
    argp_deg (float): Argument of perigee in degrees.
    start_time (Time): Start time as an astropy Time object.
    
    Returns:
    List[str]: List of TLE strings.
    """
    tle_list = []
    altitude_m = altitude_km * 1000  # Convert altitude to meters
    step = 360.0 / num_sats_per_plane
    for plane_idx in range(plane_count):
        raan_deg_plane = plane_idx * 10.0
        for satellite_idx in range(num_sats_per_plane):
            anomaly_deg_sat = satellite_idx * step
            sat_name = f"SystemC_Belt_1_Plane_{plane_idx+1}_Satellite_{satellite_idx+1}"
            tle_str = forge_tle_single(
                sat_name, altitude_m, eccentricity,
                inclination_deg, raan_deg_plane, argp_deg,
                anomaly_deg_sat, start_time
            )
            tle_list.append(tle_str)
    return tle_list
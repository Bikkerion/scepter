from astropy import time
from datetime import datetime
import numpy as np
from scepter import tleforger

# ========================
# MAIN SCRIPT STARTS HERE
# Chapter I. Calculations
# ========================

# --------------------------------------------------------
# 1) Time range
# --------------------------------------------------------
start_time = time.Time(datetime(2025, 1, 1, 0, 0, 0))
timestep = 10 # steps in seconds
td = time.TimeDelta(np.arange(0, .005*3600*24, timestep), format='sec')  # x days above steps
times = start_time + td
n_times = len(times)  # total frames

# --------------------------------------------------------
# 2) Satellite batch creation
# --------------------------------------------------------
sat_name='Test'
altitude_m=400000
eccentricity=0
inclination_deg=90
raan_deg=0
argp_deg=0
anomaly_deg=0
tle=tleforger.forge_tle_single(sat_name, altitude_m, eccentricity, inclination_deg, raan_deg, argp_deg, anomaly_deg, start_time)
print(tle)
# num_sats_per_plane = 40
# plane_count = 18
# altitude_km = 1200
# eccentricity = 0.0
# inclination_deg = 87.9
# argp_deg = 0.0

# tle_list=tleforger.forge_tle_belt(num_sats_per_plane, plane_count, altitude_km, eccentricity, inclination_deg, argp_deg, start_time)
# print(tle_list)
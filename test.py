from astropy import time
from datetime import datetime
import numpy as np
import cysgp4
from scepter import tleforger
from astropy import units as u

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
# 2) RAS Observatory definition
# --------------------------------------------------------
latitude = -30.712777 * u.deg
longitude = 21.443611 * u.deg
elevation = 1052.0 * u.m
SKAO=cysgp4.PyObserver(longitude.value,latitude.value,elevation.to(u.km).value)

# --------------------------------------------------------
# 3) Satellite batch creation
# --------------------------------------------------------
belt_name='SystemC_Belt_1'
num_sats_per_plane = 40
plane_count = 18
altitude_km = 1200
eccentricity = 0.0
inclination_deg = 87.9
argp_deg = 0.0

tle_list=tleforger.forge_tle_belt(belt_name=belt_name, num_sats_per_plane=num_sats_per_plane, plane_count=plane_count, altitude_m=altitude_km*1000, eccentricity=eccentricity, inclination_deg=inclination_deg, argp_deg=argp_deg)

# --------------------------------------------------------
# 4) Processing positions
# --------------------------------------------------------

results = cysgp4.propagate_many(
    mjds=times.mjd,
    tles=tle_list,
    observers=SKAO,
    do_eci_pos=True,
    do_eci_vel=False,
    do_geo=True,
    do_topo=True,
    do_obs_pos=False,
    do_sat_azel=False,
    do_sat_rotmat=False,
    sat_frame='zxy',
    on_error='raise',
    method='dwarner'
)
from astropy import time
from datetime import datetime
import numpy as np

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
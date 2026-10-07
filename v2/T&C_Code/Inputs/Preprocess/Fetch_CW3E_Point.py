"""Fetch_CW3E_Point.py

The CW3E 1-km, 1-hourly meteorological forcing (Pan 2025, doi
10.5281/zenodo.14714512; downscaled from NLDAS-2/HRRR, Stage IV and PRISM)
for the conus2 grid cell holding one point, from HydroFrame's hf_hydrodata
service, written as the point time series Build_TC_Forcing_From_NLDAS2
reads: one text file per variable, a header of key=value lines, then
"yyyy-mm-ddTHH:MM:SSZ,value" lines.

Needs hf_hydrodata and a HydroFrame account whose API PIN is registered on
this machine (hf.register_api_pin, run by the account holder); the PIN never
passes through this script.

usage:
  python Fetch_CW3E_Point.py LAT LON START END OUT_DIR
    START, END  UTC dates (yyyy-mm-dd); END is exclusive
    OUT_DIR     an empty or new folder

The service labels each hour by its start (the record starts at START);
the time convention of the values -- instantaneous states, precipitation
over the hour -- is checked against the NLDAS-2 forcing before building
(the label of an hour here is the service's, written as given).
"""
import datetime as dt
import os
import sys

import numpy as np
import hf_hydrodata as hf

# CW3E name in hf_hydrodata, the builder's name, its units
VARIABLES = [
    ('air_temp', 'Tair', 'K'),
    ('specific_humidity', 'Qair', 'kg/kg'),
    ('atmospheric_pressure', 'PSurf', 'Pa'),
    ('east_windspeed', 'Wind_E', 'm/s'),
    ('north_windspeed', 'Wind_N', 'm/s'),
    ('downward_shortwave', 'SWdown', 'W/m2'),
    ('precipitation', 'Rainf', 'kg/m2'),
]


def year_windows(start, end):
    """[start, end) in calendar-year pieces, so no request is large."""
    t = start
    while t < end:
        nxt = min(dt.datetime(t.year + 1, 1, 1), end)
        yield t, nxt
        t = nxt


def fetch(variable, x, y, start, end):
    values = []
    for a, b in year_windows(start, end):
        data = hf.get_gridded_data({
            'dataset': 'CW3E', 'grid': 'conus2', 'variable': variable,
            'temporal_resolution': 'hourly',
            'start_time': a.strftime('%Y-%m-%d %H:%M:%S'),
            'end_time': b.strftime('%Y-%m-%d %H:%M:%S'),
            'grid_bounds': [x, y, x + 1, y + 1]})
        data = np.asarray(data, dtype=float).reshape(-1)
        hours = int((b - a).total_seconds() // 3600)
        if data.size != hours:
            raise RuntimeError('%s %s-%s: %d values for %d hours'
                               % (variable, a, b, data.size, hours))
        values.append(data)
        print('  %s %s: %d hours' % (variable, a.year, hours), flush=True)
    return np.concatenate(values)


def main(argv):
    if len(argv) != 6:
        sys.exit(__doc__)
    lat, lon = float(argv[1]), float(argv[2])
    start = dt.datetime.strptime(argv[3], '%Y-%m-%d')
    end = dt.datetime.strptime(argv[4], '%Y-%m-%d')
    out_dir = argv[5]
    if end <= start:
        sys.exit('END must follow START.')
    os.makedirs(out_dir, exist_ok=True)
    if os.listdir(out_dir):
        sys.exit('%s is not empty.' % out_dir)
    xf, yf = hf.from_latlon('conus2', lat, lon)
    x, y = int(np.floor(xf)), int(np.floor(yf))
    clat, clon = hf.to_latlon('conus2', x + 0.5, y + 0.5)
    print('conus2 cell x=%d y=%d, centre %.5f %.5f' % (x, y, clat, clon), flush=True)
    hours = int((end - start).total_seconds() // 3600)
    stamps = [(start + dt.timedelta(hours=h)).strftime('%Y-%m-%dT%H:%M:%SZ')
              for h in range(hours)]
    for variable, name, units in VARIABLES:
        v = fetch(variable, x, y, start, end)
        if not np.all(np.isfinite(v)):
            raise RuntimeError('%s holds non-finite values' % variable)
        path = os.path.join(out_dir, 'CW3E_%s.csv' % name)
        with open(path, 'w', newline='\n') as f:
            f.write('variable=%s\nunits=%s\nlat=%.6f\nlon=%.6f\n' % (name, units, clat, clon))
            f.write('source=CW3E 1-km 1-hourly forcing (hf_hydrodata dataset CW3E, '
                    'conus2 x=%d y=%d), variable %s\n' % (x, y, variable))
            for s, value in zip(stamps, v):
                f.write('%s,%.7g\n' % (s, value))
        print('wrote %s' % path, flush=True)


if __name__ == '__main__':
    main(sys.argv)

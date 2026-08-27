#!/usr/bin/env python3
"""Build a small synthetic netCDF climate file for testing cli_ncdf_meas.

The data encodes its own date so a misplaced record is visible in the output:

    pcp(t) = 0.01 * (year - 1900) + 0.001 * day_of_year        [mm]

Values are uniform in space, so the nearest-grid-cell search cannot influence
the result and any HRU's annual precipitation equals the annual total below.
Carrying the year matters: a day-of-year-only encoding would hide a whole-year
shift between two leap years, which is exactly how the original bug survived
against the real 20crv3-era5 file (1960 and 1980 are both leap years).

Deliberate omissions and defects, to exercise the missing-data path:
  * `wnd` is not written at all         -> reader must fall back to the generator
  * `hmd` is NaN for 1980-06-01..30     -> NaN must reach the generator, not 0.0
  * `slr` is -9999.0 for 1981-02-01..28 -> negative fill must reach the generator

Usage:
    python3 test/make_ncdf_fixture.py OUT.nc [--start 1978-01-01] [--end 1982-12-31]
"""

import argparse
import datetime as dt

import netCDF4
import numpy as np

REF = dt.date(1900, 1, 1)


def parse_date(s):
    return dt.date(*[int(p) for p in s.split("-")])


def day_of_year(d):
    return d.timetuple().tm_yday


def pcp_value(d):
    return 0.01 * (d.year - 1900) + 0.001 * day_of_year(d)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--start", default="1978-01-01")
    ap.add_argument("--end", default="1982-12-31")
    args = ap.parse_args()

    start, end = parse_date(args.start), parse_date(args.end)
    dates = [start + dt.timedelta(days=i) for i in range((end - start).days + 1)]

    # same grid as 20crv3_era5_1960_2021.nc
    lat = 40.0 - 0.5 * np.arange(22) - 0.25
    lon = -94.0 + 0.5 * np.arange(11) + 0.25
    nt, nlat, nlon = len(dates), len(lat), len(lon)

    ds = netCDF4.Dataset(args.out, "w", format="NETCDF4")
    ds.createDimension("time", nt)
    ds.createDimension("lat", nlat)
    ds.createDimension("lon", nlon)

    v = ds.createVariable("time", "f8", ("time",))
    v.units = "days since 1900-01-01"
    v.calendar = "proleptic_gregorian"
    v[:] = np.array([(d - REF).days for d in dates], dtype="f8")

    v = ds.createVariable("lat", "f8", ("lat",))
    v.units = "degrees_north"
    v[:] = lat

    v = ds.createVariable("lon", "f8", ("lon",))
    v.units = "degrees_east"
    v[:] = lon

    def field(name, values, units, fill=np.nan):
        var = ds.createVariable(name, "f4", ("time", "lat", "lon"), fill_value=fill)
        var.units = units
        # same value in every cell. netCDF4 1.7.4 emits a numpy 2.5
        # DeprecationWarning from inside its own __setitem__ here; harmless.
        var[:] = values[:, None, None] * np.ones((1, nlat, nlon), dtype="f4")

    pcp = np.array([pcp_value(d) for d in dates], dtype="f4")
    field("pcp", pcp, "mm")
    field("tmax", np.full(nt, 25.0, dtype="f4"), "degC")
    field("tmin", np.full(nt, 15.0, dtype="f4"), "degC")

    hmd = np.full(nt, 0.6, dtype="f4")
    hmd[[i for i, d in enumerate(dates) if d.year == 1980 and d.month == 6]] = np.nan
    field("hmd", hmd, "mm/mm")

    slr = np.full(nt, 18.0, dtype="f4")
    slr[[i for i, d in enumerate(dates) if d.year == 1981 and d.month == 2]] = -9999.0
    field("slr", slr, "MJ m-2 day-1")

    # wnd intentionally absent

    ds.close()

    print(f"wrote {args.out}: {nt} days, {args.start}..{args.end}, grid {nlat}x{nlon}")
    print("expected annual precipitation totals (mm), per station and per HRU:")
    for year in sorted({d.year for d in dates}):
        total = sum(pcp_value(d) for d in dates if d.year == year)
        print(f"  {year}: {total:.3f}")

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

With --atmodep, also writes a monthly atmospheric deposition file on the same
grid, again uniform in space, with each value carrying its own running month
index k (k=1 is January of the start year):

    nh4_rf(k)  = 0.10 + 0.01 * k    [mg/L]
    no3_rf(k)  = 0.20 + 0.01 * k    [mg/L]
    nh4_dry(k) = 0.30 + 0.01 * k    [kg/ha/month]
    no3_dry(k) = 0.40 + 0.01 * k    [kg/ha/month]

The 0.01 step is chosen against the output format: hru_nb_*.txt prints these
through f17.3, so a 0.001 step would put the month index in the last printed
digit, where float32 rounding can swallow an off-by-one. At 0.01 a wrong month
is ten times the print resolution.

Usage:
    python3 test/make_ncdf_fixture.py OUT.nc [--start 1978-01-01] [--end 1982-12-31]
                                             [--atmodep DEP.nc]
"""

import argparse
import datetime as dt
import warnings

import netCDF4
import numpy as np

# netCDF4 1.7.4 sets .shape on a numpy array inside its own Variable.__setitem__
# (_netCDF4.pyx:5633), which numpy 2.5 deprecates. Nothing this script passes in
# can avoid it and it fires once per variable written, so silence just that one.
warnings.filterwarnings(
    "ignore", category=DeprecationWarning,
    message="Setting the shape on a NumPy array")

REF = dt.date(1900, 1, 1)


def parse_date(s):
    return dt.date(*[int(p) for p in s.split("-")])


def day_of_year(d):
    return d.timetuple().tm_yday


def pcp_value(d):
    return 0.01 * (d.year - 1900) + 0.001 * day_of_year(d)


# base value, units; the running month index k is added on at 0.01 per month
DEP_VARS = [("nh4_rf", 0.10, "mg L-1"),
            ("no3_rf", 0.20, "mg L-1"),
            ("nh4_dry", 0.30, "kg month-1 ha-1"),
            ("no3_dry", 0.40, "kg month-1 ha-1")]


def dep_value(base, k):
    return base + 0.01 * k


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--start", default="1978-01-01")
    ap.add_argument("--end", default="1982-12-31")
    ap.add_argument("--atmodep", default=None,
                    help="also write a monthly deposition fixture to this path")
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
        # same value in every cell
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

    # Monthly deposition fixture on the same grid. Uniform in space, so the
    # nearest-cell search cannot change the answer, and carrying the running
    # month index k means a wrong atmodep_cont%ts lands on a visibly wrong
    # number instead of a plausible one.
    if args.atmodep:
        months = [(y, m) for y in range(start.year, end.year + 1)
                  for m in range(1, 13)]
        epoch = dt.date(start.year, 1, 1)
        offsets = [(dt.date(y, m, 1) - epoch).days for (y, m) in months]

        dep = netCDF4.Dataset(args.atmodep, "w", format="NETCDF4")
        dep.createDimension("time", len(months))
        dep.createDimension("lat", nlat)
        dep.createDimension("lon", nlon)

        dv = dep.createVariable("time", "i8", ("time",))
        dv.units = f"days since {start.year}-01-01 00:00:00"
        dv.calendar = "proleptic_gregorian"
        dv[:] = offsets

        dv = dep.createVariable("lat", "f8", ("lat",))
        dv.units = "degrees_north"
        dv[:] = lat

        dv = dep.createVariable("lon", "f8", ("lon",))
        dv.units = "degrees_east"
        dv[:] = lon

        k = np.arange(1, len(months) + 1)
        for name, base, units in DEP_VARS:
            dv = dep.createVariable(name, "f4", ("time", "lat", "lon"))
            dv.units = units
            vals = dep_value(base, k).astype("f4")
            dv[:] = vals[:, None, None] * np.ones((1, nlat, nlon), dtype="f4")

        dep.close()

        print(f"wrote {args.atmodep}: {len(months)} months, "
              f"{start.year}-01 .. {end.year}-12, grid {nlat}x{nlon}")
        print("expected monthly deposition, per station and per HRU:")
        print("  " + "  ".join(f"{h:>9}" for h in
                               ["month", "k"] + [n for n, _, _ in DEP_VARS]))
        for i, (y, m) in enumerate(months, start=1):
            row = "  ".join(f"{dep_value(b, i):9.3f}" for _, b, _ in DEP_VARS)
            print(f"  {y}-{m:02d}      {i:9d}  {row}")
        print("hru_nb_mon.txt check, with P_mon from hru_wb_mon.txt field 8:")
        print("  no3atmo = 0.01 * no3_rf(k) * P_mon + no3_dry(k)")
        print("  nh4atmo = 0.01 * nh4_rf(k) * P_mon + nh4_dry(k)")

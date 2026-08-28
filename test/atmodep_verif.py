#!/usr/bin/env python3
"""Verify SWAT+ netCDF atmospheric deposition against the source netCDF file.

Picks HRUs at random, follows each one to its weather station and from there to
the deposition grid cell the reader would have used, then checks every month of
monthly HRU output against a value computed independently from atmodep.nc.

Why monthly rather than daily. nut_nrain (src/nut_nrain.f90:49-52) computes, for
each day of a month with N days,

    atmo_day = 0.01 * rf(k) * precip_day + dry(k) / N

and hru_output.f90:42 accumulates the monthly total as a plain sum with no
month-end division, so the /N cancels over a whole month and the dry term
becomes a pure additive constant. No leap-year dependence, and a wrong record
index shifts a month by a visible amount.

The record index k. cli_ncdf_read_atmodep sets atmodep_cont%first = 0 and seeds
ts with atmodep_start_index(), after which cli_atmodep_time_control increments it
once per month, so for a monthly record

    k = 12 * (yr - yr_init) + (mo - mo_init) + 1

with yr_init/mo_init decoded from the deposition file's own time axis. k is an
index into the whole record, not into the simulated period. nut_nrain zeroes the
deposition when k falls outside 1..ntime.

THE LAST DAY OF EACH MONTH USES THE NEXT MONTH'S RECORD. In time_control.f90,
time%end_mo is set at line 182 and cli_atmodep_time_control increments ts at
line 230, but hru_control -> nut_nrain does not run until `command` at line 249.
So on day N the counter has already moved on and nut_nrain reads record k+1:

    atmo_mon = 0.01 * rf(k)   * (P_mon - p_last) + dry(k)   * (N-1)/N
             + 0.01 * rf(k+1) * p_last           + dry(k+1) / N

This is upstream SWAT+ behaviour, present on the text atmodep.cli path too
(cli_atmodep_time_control.f90 is untouched by the netCDF work), so this script
reproduces it rather than flagging it. It matters only when a month's last day
is wet AND the deposition jumps between consecutive months -- with real ISIMIP
data that reaches 0.08 kg/ha, well outside any rounding.

Getting p_last therefore means sampling the climate netCDF, which this script
does at the station's own grid cell. That sampling is self-checked: the daily
precip summed over each month must reproduce P_mon from hru_wb_mon.txt, and a
month where it does not is reported as UNVERIFIED rather than quietly passed.
Note the climate reader uses a 2-D Euclidean nearest cell (cli_ncdf_meas.f90:288)
while the deposition reader uses a separable per-axis nearest -- the same answer
on a regular grid, but they are mirrored separately here rather than shared.

Usage:
    python3 test/atmodep_verif.py path/to/TxtInOut [-n 2] [--seed 7]
    python3 test/atmodep_verif.py path/to/TxtInOut --hru hru0001 --hru hru0042

Requires monthly output for both objects in print.prt:
    hru_wb   ... monthly = y      (precip, field 8)
    hru_nb   ... monthly = y      (no3atmo/nh4atmo, fields 21/22)

Assumes plaps = 0 in parameters.bsn, so HRU precip equals station precip; a
nonzero lapse rate shows up as a failed precip self-check, not a silent error.
"""

import argparse
import calendar
import datetime as dt
import pathlib
import random
import re
import sys

import netCDF4
import numpy as np

# Column indices in the monthly output files, 0-based.
WB_PRECIP = 7                    # hru_wb_mon.txt  field 8
NB_NO3ATMO, NB_NH4ATMO = 20, 21  # hru_nb_mon.txt  fields 21, 22
OUT_NAME, OUT_YR, OUT_MO = 6, 3, 1

# hru.con, 0-based: 0 id, 1 name, 2 gis_id, 3 area, 4 lat, 5 lon, 6 elev,
#                   7 hru, 8 wst
CON_NAME, CON_WST = 1, 8

DEP_VARS = ("nh4_rf", "no3_rf", "nh4_dry", "no3_dry")
DEP_SKIP = -99.0   # an .ncw factor of -99 means "no deposition on this station"
MISSING = -99.0    # SWAT+'s "not measured, use the weather generator" sentinel


def die(msg):
    sys.exit(f"atmodep_verif: {msg}")


def data_lines(path, skip=3):
    """Output-file rows, dropping the header block and any short trailing line."""
    with open(path) as fh:
        for line in fh.readlines()[skip:]:
            f = line.split()
            if f and f[0].lstrip("-").isdigit():
                yield f


def read_cio(txtinout):
    """netcdf.ncw (climate slot 1), the deposition file (slot 9) and pcp_path."""
    cio = (txtinout / "file.cio").read_text().splitlines()
    line = next((l for l in cio if l.startswith("climate")), None)
    if line is None:
        die("no 'climate' line in file.cio")
    f = line.split()
    if len(f) < 10:
        die(f"climate line has {len(f)} fields, need at least 10")
    if not f[1].endswith(".ncw"):
        die(f"climate slot 1 is {f[1]!r}, not a .ncw file -- not a netCDF run")
    if f[9] in ("null", "0"):
        die("file.cio climate slot 9 is 'null': deposition is not configured")

    # Only pcp_path is read by cli_ncdf_meas; the other four are parsed and
    # discarded, so all climate variables live in the file it names.
    pcp = next((l.split()[1] for l in cio
                if l.startswith("pcp_path") and len(l.split()) > 1), None)
    if pcp is None or pcp == "null":
        die("no usable pcp_path line in file.cio")
    return f[1], f[9], pcp


def read_ncw(path):
    """Station name -> {lat, lon, factors}, keyed off the header, not a count."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    if len(lines) < 3:
        die(f"{path.name} has no station rows")
    head = lines[1].split()

    def col(name):
        if name not in head:
            die(f"{path.name} has no {name!r} column; header is {head}")
        return head.index(name)

    ilat, ilon, ipcp = col("latitude"), col("longitude"), col("pcp")
    ifac = {v: col(v) for v in DEP_VARS}

    stations = {}
    for row in lines[2:]:
        f = row.split()
        if len(f) <= max([ilat, ilon, ipcp] + list(ifac.values())):
            die(f"{path.name}: station row has {len(f)} fields, header has "
                f"{len(head)} -- column count mismatch")
        try:
            stations[f[0]] = dict(
                lat=np.float32(f[ilat]), lon=np.float32(f[ilon]),
                pcp_factor=np.float32(f[ipcp]),
                fac={v: np.float32(f[i]) for v, i in ifac.items()})
        except ValueError as e:
            die(f"{path.name}: non-numeric field in row for {f[0]!r}: {e}")
    return stations


def read_hru_con(path):
    """HRU name -> weather station name."""
    lines = [l for l in path.read_text().splitlines() if l.strip()]
    out = {}
    for row in lines[2:]:
        f = row.split()
        if len(f) > CON_WST:
            out[f[CON_NAME]] = f[CON_WST]
    if not out:
        die(f"{path.name}: no HRU rows found")
    return out


def parse_time_units(units):
    """Epoch date out of a CF 'days since YYYY-MM-DD ...' units string."""
    m = re.search(r"(\d{1,4})-(\d{1,2})-(\d{1,2})", units)
    if not m:
        die(f"cannot parse a date out of time units {units!r}")
    if "days since" not in units.replace("  ", " ").lower():
        print(f"  ! warning: time units are {units!r}; the reader assumes days")
    return dt.date(*(int(g) for g in m.groups()))


def classify_timestep(offsets):
    """Mirrors cli_ncdf_atmodep_util::classify_timestep."""
    if len(offsets) <= 1:
        return "aa"
    span = (float(offsets[-1]) - float(offsets[0])) / (len(offsets) - 1)
    if 20.0 <= span <= 40.0:
        return "mo"
    if 350.0 <= span <= 380.0:
        return "yr"
    return "??"


def nearest_cell(vals, target):
    """Mirrors cli_ncdf_atmodep_util::nearest_cell: 1-D nearest, first wins ties.

    Kept as an explicit loop with a strict '<' rather than argmin, so a tie
    breaks the same way the Fortran does.
    """
    best_i, best_d = 0, abs(np.float32(vals[0]) - target)
    for i in range(1, len(vals)):
        d = abs(np.float32(vals[i]) - target)
        if d < best_d:
            best_i, best_d = i, d
    return best_i, best_d


def dep_value(raw, factor):
    """Mirrors cli_ncdf_read_atmodep::dep_value."""
    if factor <= DEP_SKIP + 1.0:
        return 0.0
    if np.isnan(raw):
        return 0.0
    if raw < 0.0:
        return 0.0
    return float(raw) * float(factor)


def nc_value(raw, factor, clamp_low=0.0, clamp_high=None):
    """Mirrors cli_ncdf_meas::nc_value. Returns -99 for 'missing, use the wgn'."""
    if np.isnan(raw):
        return MISSING
    if raw <= -97.0:
        return MISSING
    v = float(raw) * float(factor)
    if clamp_low is not None and v < clamp_low:
        v = clamp_low
    if clamp_high is not None and v > clamp_high:
        v = clamp_high
    return v


def nearest_cell_2d(lat_vals, lon_vals, tgt_lat, tgt_lon):
    """Mirrors cli_ncdf_meas.f90:286-292: 2-D Euclidean nearest, lat-major scan.

    Not the same routine as the deposition reader's separable search. On a
    regular lat/lon grid both land on the same cell, but they are mirrored
    separately so that a future divergence is caught rather than assumed away.
    """
    best = (None, None, float("inf"))
    for a in range(len(lat_vals)):
        dlat = float(lat_vals[a]) - float(tgt_lat)
        for o in range(len(lon_vals)):
            dlon = float(lon_vals[o]) - float(tgt_lon)
            d = (dlat * dlat + dlon * dlon) ** 0.5
            if d < best[2]:
                best = (a, o, d)
    return best


class ClimateSampler:
    """Daily station precipitation, read from the climate netCDF as SWAT+ does.

    cli_ncdf_meas places each record on the calendar date its own time value
    decodes to (populate_timeseries_data), then scales it by the station's pcp
    factor through nc_value. Reproduced here so the last day of a month, which
    the deposition counter has already moved past, can be priced correctly.
    """

    def __init__(self, path):
        self.ds = netCDF4.Dataset(path)
        for v in ("pcp", "lat", "lon", "time"):
            if v not in self.ds.variables:
                die(f"{path.name} has no {v!r} variable")
        self.lat = np.asarray(self.ds["lat"][:], dtype="f4")
        self.lon = np.asarray(self.ds["lon"][:], dtype="f4")
        epoch = parse_time_units(self.ds["time"].units)
        offsets = np.asarray(self.ds["time"][:])
        self.date_index = {epoch + dt.timedelta(days=int(o)): i
                           for i, o in enumerate(offsets)}
        self._cache = {}

    def cell(self, st):
        return nearest_cell_2d(self.lat, self.lon, st["lat"], st["lon"])

    def series(self, st):
        """Whole-record precip at this station's cell, scaled and clamped."""
        a, o, _ = self.cell(st)
        if (a, o) not in self._cache:
            raw = np.asarray(self.ds["pcp"][:, a, o], dtype="f4")
            self._cache[(a, o)] = np.array(
                [nc_value(r, st["pcp_factor"]) for r in raw])
        return self._cache[(a, o)]

    def month(self, st, yr, mo):
        """(daily precip list, n_missing) for one month; None where uncovered."""
        s = self.series(st)
        vals, missing = [], 0
        for day in range(1, calendar.monthrange(yr, mo)[1] + 1):
            i = self.date_index.get(dt.date(yr, mo, day))
            if i is None or s[i] == MISSING:
                missing += 1
                vals.append(None)
            else:
                vals.append(float(s[i]))
        return vals, missing

    def close(self):
        self.ds.close()


def main():
    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__)
    ap.add_argument("txtinout", type=pathlib.Path)
    ap.add_argument("-n", type=int, default=2,
                    help="number of random HRUs to check (default 2)")
    ap.add_argument("--hru", action="append", default=[],
                    help="check this HRU by name; repeatable, disables -n")
    ap.add_argument("--seed", type=int, default=None,
                    help="seed for the random pick, for a reproducible run")
    ap.add_argument("--tol", type=float, default=2e-3,
                    help="max allowed |obs-exp| in kg/ha (default 0.002)")
    args = ap.parse_args()

    d = args.txtinout
    if not d.is_dir():
        die(f"{d} is not a directory")

    ncw_name, dep_name, pcp_name = read_cio(d)
    stations = read_ncw(d / ncw_name)
    hru2wst = read_hru_con(d / "hru.con")
    if not (d / pcp_name).exists():
        die(f"climate netCDF {pcp_name} not found (needed for the last day of "
            "each month)")
    clim = ClimateSampler(d / pcp_name)

    for f in ("hru_wb_mon.txt", "hru_nb_mon.txt"):
        if not (d / f).exists():
            die(f"{f} not found. Set the monthly column to 'y' for hru_wb and "
                "hru_nb in print.prt and rerun the model.")

    precip, atmo = {}, {}
    for f in data_lines(d / "hru_wb_mon.txt"):
        precip[(f[OUT_NAME], int(f[OUT_YR]), int(f[OUT_MO]))] = float(f[WB_PRECIP])
    for f in data_lines(d / "hru_nb_mon.txt"):
        atmo[(f[OUT_NAME], int(f[OUT_YR]), int(f[OUT_MO]))] = (
            float(f[NB_NO3ATMO]), float(f[NB_NH4ATMO]))
    if not atmo:
        die("hru_nb_mon.txt has no data rows")

    ds = netCDF4.Dataset(d / dep_name)
    lat_vals = np.asarray(ds["lat"][:], dtype="f4")
    lon_vals = np.asarray(ds["lon"][:], dtype="f4")
    offsets = np.asarray(ds["time"][:])
    ntime = len(offsets)

    ts = classify_timestep(offsets)
    if ts != "mo":
        die(f"deposition record is {ts!r}; this script checks monthly records "
            "only (the reader also supports 'aa' and 'yr')")

    epoch = parse_time_units(ds["time"].units)
    first = epoch + dt.timedelta(days=int(offsets[0]))
    yr_init, mo_init = first.year, first.month

    print(f"deposition file : {dep_name}")
    print(f"  record        : {ntime} monthly steps from "
          f"{yr_init:04d}-{mo_init:02d}")
    print(f"  grid          : {len(lat_vals)} lat x {len(lon_vals)} lon")
    print(f"  time units    : {ds['time'].units}")

    # Only HRUs that actually appear in both output files are checkable.
    have = sorted({k[0] for k in atmo} & {k[0] for k in precip} & set(hru2wst))
    if not have:
        die("no HRU appears in hru_nb_mon.txt, hru_wb_mon.txt and hru.con alike")

    if args.hru:
        missing = [h for h in args.hru if h not in have]
        if missing:
            die(f"not checkable (absent from output or hru.con): {missing}")
        chosen = args.hru
    else:
        seed = args.seed if args.seed is not None else random.randrange(2**31)
        random.seed(seed)
        chosen = random.sample(have, min(args.n, len(have)))
        print(f"  random seed   : {seed}  (pass --seed {seed} to repeat)")

    # Cache the per-station series so two HRUs on one station cost one read.
    series_cache = {}

    def station_series(wst_name):
        if wst_name in series_cache:
            return series_cache[wst_name]
        st = stations.get(wst_name)
        if st is None:
            die(f"weather station {wst_name!r} is in hru.con but not in {ncw_name}")
        ilat, dlat = nearest_cell(lat_vals, st["lat"])
        ilon, dlon = nearest_cell(lon_vals, st["lon"])
        vals = {v: np.asarray(ds[v][:, ilat, ilon], dtype="f4") for v in DEP_VARS}
        series_cache[wst_name] = (st, ilat, ilon, dlat, dlon, vals)
        return series_cache[wst_name]

    worst_overall, failures, checked, unverified = 0.0, 0, 0, 0

    for hru in chosen:
        wst_name = hru2wst[hru]
        st, ilat, ilon, dlat, dlon, vals = station_series(wst_name)

        print(f"\n=== {hru}  station {wst_name} ===")
        print(f"  station at    : lat {float(st['lat']):.4f}, "
              f"lon {float(st['lon']):.4f}")
        print(f"  nearest cell  : [lat {ilat}] {float(lat_vals[ilat]):.4f}, "
              f"[lon {ilon}] {float(lon_vals[ilon]):.4f} "
              f"(off by {float(dlat):.4f}, {float(dlon):.4f} deg)")
        facs = "  ".join(f"{v}={float(st['fac'][v]):g}" for v in DEP_VARS)
        print(f"  .ncw factors  : {facs}")
        skipped = [v for v in DEP_VARS if st["fac"][v] <= DEP_SKIP + 1.0]
        if skipped:
            print(f"  ! {skipped} set to -99: the reader forces these to zero")

        ca, co, cd = clim.cell(st)
        print(f"  climate cell  : [lat {ca}] {float(clim.lat[ca]):.4f}, "
              f"[lon {co}] {float(clim.lon[co]):.4f} "
              f"(pcp factor {float(st['pcp_factor']):g})")

        months = sorted(m for m in atmo if m[0] == hru)
        print(f"\n  {'month':>8} {'k':>5} {'P_mon':>8} {'p_last':>7} "
              f"{'no3 obs':>9} {'no3 exp':>9} {'diff':>7} "
              f"{'nh4 obs':>9} {'nh4 exp':>9} {'diff':>7}")

        worst = 0.0
        for key in months:
            _, yr, mo = key
            k = 12 * (yr - yr_init) + (mo - mo_init) + 1
            n_days = calendar.monthrange(yr, mo)[1]
            p_mon = precip[key]
            obs_no3, obs_nh4 = atmo[key]

            # Self-check the climate sampling: our daily precip must add up to
            # the precip the model reported, or p_last is not to be trusted.
            daily, n_missing = clim.month(st, yr, mo)
            p_last = daily[-1]
            summed = sum(v for v in daily if v is not None)
            bad_precip = (n_missing > 0
                          or abs(summed - p_mon) > max(0.05, 2e-4 * p_mon))

            def term(kk, wet_precip, dry_share):
                """One record's contribution; zero outside the record, as
                nut_nrain's `ist > 0 .and. ist <= num` guard leaves it."""
                if not 1 <= kk <= ntime:
                    return 0.0, 0.0
                i = kk - 1
                no3 = (0.01 * dep_value(vals["no3_rf"][i], st["fac"]["no3_rf"])
                       * wet_precip
                       + dep_value(vals["no3_dry"][i], st["fac"]["no3_dry"])
                       * dry_share)
                nh4 = (0.01 * dep_value(vals["nh4_rf"][i], st["fac"]["nh4_rf"])
                       * wet_precip
                       + dep_value(vals["nh4_dry"][i], st["fac"]["nh4_dry"])
                       * dry_share)
                return no3, nh4

            if bad_precip:
                note = f"  <- UNVERIFIED ({n_missing} generated days)" \
                    if n_missing else \
                    f"  <- UNVERIFIED (precip sum {summed:.2f} != {p_mon:.2f})"
                unverified += 1
                checked += 1
                pl = -1.0 if p_last is None else p_last
                print(f"  {yr}-{mo:02d} {k:5d} {p_mon:8.2f} {pl:7.2f} "
                      f"{obs_no3:9.3f} {'':>9} {'':>7} "
                      f"{obs_nh4:9.3f} {'':>9} {'':>7}{note}")
                continue

            # Days 1..N-1 use record k; day N uses k+1, because the counter is
            # advanced before hru_control runs. See the module docstring.
            a3, a4 = term(k, p_mon - p_last, (n_days - 1) / n_days)
            b3, b4 = term(k + 1, p_last, 1.0 / n_days)
            exp_no3, exp_nh4 = a3 + b3, a4 + b4

            d3, d4 = obs_no3 - exp_no3, obs_nh4 - exp_nh4
            worst = max(worst, abs(d3), abs(d4))
            checked += 1
            note = "" if 1 <= k <= ntime else "  <- k outside record"
            print(f"  {yr}-{mo:02d} {k:5d} {p_mon:8.2f} {p_last:7.2f} "
                  f"{obs_no3:9.3f} {exp_no3:9.3f} {d3:7.3f} "
                  f"{obs_nh4:9.3f} {exp_nh4:9.3f} {d4:7.3f}{note}")

        ok = worst <= args.tol
        failures += 0 if ok else 1
        worst_overall = max(worst_overall, worst)
        print(f"  worst |obs-exp| = {worst:.4f} kg/ha  "
              f"{'PASS' if ok else 'FAIL'}")

    ds.close()
    clim.close()

    print(f"\n{checked} month-values checked across {len(chosen)} HRU(s); "
          f"worst |obs-exp| = {worst_overall:.4f} kg/ha (tolerance {args.tol})")
    if unverified:
        print(f"  {unverified} month(s) UNVERIFIED: the climate sampling did "
              "not reproduce the reported precip, so p_last is unknown there.")
    if failures:
        print(f"FAIL: {failures} of {len(chosen)} HRU(s) outside tolerance")
        print("  Same offset every month      -> the record index k")
        print("  Varies month to month        -> the grid cell or the factors")
        print("  Only months with a wet 31st  -> the last-day rule in the "
              "docstring has changed")
        return 1
    print("PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())

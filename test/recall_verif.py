"""
Check SWAT+ recall output against the orgmin.nc it was fed.

A point source's yearly recall output should equal that year's sum of the rows the netCDF
holds for it. Run against a scenario directory that has already been simulated:

    python test/recall_verif.py data/mississippi_downstream/rec_netcdf

recall_yr.txt is laid out as a title line, a line of column names, a line of units, then
the data. Its flo units read m^3/s but the values are the summed m3 the input carried, so
the comparison is a plain sum with no conversion.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

TOL = 1e-3      # relative; float32 on both sides

# Output column -> netCDF variable. The rest of the hyd_output fields are spelled the same
# in both places; these four are not.
RENAMES = {"san": "sand", "sil": "silt", "cla": "clay", "grv": "gravel"}

CHECK = ["flo", "sed", "orgn", "sedp", "no3", "solp", "chla", "nh3", "no2", "cbod", "dox"]


def main(scen_dir):
    scen = Path(scen_dir)

    nc_FN = scen / "orgmin.nc"
    out_FN = scen / "recall_yr.txt"
    for f in (nc_FN, out_FN):
        if not f.exists():
            print(f"missing {f}")
            return 1

    ds = xr.open_dataset(nc_FN).set_index(row=["name", "time"])

    lines = out_FN.read_text().splitlines()
    cols = lines[1].split()
    df = pd.read_csv(out_FN, skiprows=3, sep=r"\s+", header=None,
                     names=cols, usecols=range(len(cols)))

    names = sorted(set(np.asarray(ds.name.values).tolist()))
    print(f"{len(names)} record(s) in {nc_FN.name}, {len(df)} row(s) in {out_FN.name}\n")

    failures = 0
    for name in names:
        rec = ds.sel(name=name)
        rows = df[df["name"] == name]

        if rows.empty:
            print(f"  {name}: no rows in {out_FN.name}")
            failures += 1
            continue

        bad_cols = set()
        for _, row in rows.iterrows():
            yr = int(row["yr"])
            for col in CHECK:
                if col not in df.columns:
                    continue
                want = float(rec.sel(time=str(yr))[RENAMES.get(col, col)].sum())
                have = float(row[col])
                if abs(have - want) > TOL * max(1.0, abs(want)):
                    if col not in bad_cols:
                        print(f"  {name} {yr} {col}: out={have:.6E} nc={want:.6E}")
                    bad_cols.add(col)

        if bad_cols:
            failures += 1
            print(f"  {name}: {len(bad_cols)} column(s) disagree: {sorted(bad_cols)}")
        else:
            print(f"  {name}: {len(rows)} year(s) match across {len(CHECK)} columns")

    print("\nrecall output matches orgmin.nc" if not failures
          else f"\n{failures} record(s) disagree")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))

#!/usr/bin/env bash
# Run the SWAT+ netCDF reader against a synthetic fixture and report hru0001's
# annual precipitation, which the fixture makes analytically predictable.
#
# Also checks hru0001's monthly atmospheric deposition against the deposition
# fixture, unless NCDF_DEP=0, which instead strips the deposition wiring to
# exercise the backward-compatible 12-column path.
#
# Usage: test/run_ncdf_fixture.sh <scratch-dir> [start-yr] [end-yr]
#        NCDF_DEP=0 test/run_ncdf_fixture.sh <scratch-dir>
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$REPO/data/mississippi_downstream/TxtInOut"
DST="${1:?usage: run_ncdf_fixture.sh <scratch-dir> [start-yr] [end-yr]}"
YR0="${2:-1980}"
YR1="${3:-1981}"
DEP="${NCDF_DEP:-1}"

# The executable name embeds `git describe --tags`, so it changes with every
# commit. Take the most recently built match rather than hard-coding a name.
EXE="$(ls -t "$REPO"/build/debug_nc/swatplus-*-Dbg-nc 2>/dev/null | head -1)"
[ -n "$EXE" ] || { echo "no swatplus-*-Dbg-nc in $REPO/build/debug_nc" >&2; exit 1; }

# ifx-built binaries need the Intel runtime on the library path
export LD_LIBRARY_PATH="/home/jopato/intel/oneapi/compiler/2026.1/lib:${LD_LIBRARY_PATH:-}"

echo "exe:     $EXE"
echo "scratch: $DST"

rm -rf "$DST"
mkdir -p "$DST"
cp -r "$SRC"/. "$DST"/
rm -f "$DST"/20crv3_era5_1960_2021.nc "$DST"/atmodep.nc

if [ "$DEP" = 1 ]; then
  python3 "$REPO/test/make_ncdf_fixture.py" "$DST/fixture.nc" \
          --atmodep "$DST/atmodep.nc"
else
  python3 "$REPO/test/make_ncdf_fixture.py" "$DST/fixture.nc"
fi

# file.cio carries the netCDF path on its pcp_path line (and tmp/slr/hmd/wnd,
# which cli_ncdf_meas ignores -- only pcp_path is ever read)
sed -i 's#20crv3_era5_1960_2021\.nc#fixture.nc#' "$DST/file.cio"

# The deposition check reads hru_nb_mon.txt and hru_wb_mon.txt, but print.prt
# ships with the monthly column off, so turn it on for those two objects.
#
# With NCDF_DEP=1 also assert the deposition wiring is in place -- a netcdf.ncw
# without the four columns, or a file.cio without atmodep.nc in slot 9, runs
# fine and silently skips deposition, which would look like a passing test.
# With NCDF_DEP=0 do the opposite and strip both, since the shipped project now
# has deposition configured; that is the backward-compatible 12-column path.
python3 - "$DST" "$DEP" <<'PY'
import re, sys, pathlib
d = pathlib.Path(sys.argv[1])
dep = sys.argv[2] == "1"
COLS = ("nh4_rf", "no3_rf", "nh4_dry", "no3_dry")

ncwp = d / "netcdf.ncw"
ncw = ncwp.read_text().splitlines()
cio_p = d / "file.cio"
cio = cio_p.read_text().splitlines()
ic = next(i for i, l in enumerate(cio) if l.startswith("climate"))

if dep:
    missing = [c for c in COLS if c not in ncw[1]]
    if missing:
        sys.exit(f"netcdf.ncw is missing deposition columns: {missing}")
    if cio[ic].split()[9] != "atmodep.nc":
        sys.exit(f"file.cio slot 9 is {cio[ic].split()[9]!r}, "
                 "expected 'atmodep.nc'")
else:
    # Drop the last four fields of the header and of every station row, then
    # blank slot 9. Both files are whitespace-delimited and read list-directed.
    keep = len(re.findall(r"\S+", ncw[1])) - len(COLS)
    for i in range(1, len(ncw)):
        if ncw[i].strip():
            ncw[i] = "".join(f"{t:<14}" for t in ncw[i].split()[:keep])
    ncwp.write_text("\n".join(ncw) + "\n")

    f = cio[ic].split()
    f[9] = "null"
    cio[ic] = "".join(f"{x:<20}" for x in f)
    cio_p.write_text("\n".join(cio) + "\n")
    print("stripped deposition columns from netcdf.ncw and file.cio slot 9")

prt = d / "print.prt"
lines = prt.read_text().splitlines()
for i, line in enumerate(lines):
    if re.match(r"(hru_wb|hru_nb)\s", line):
        tok = list(re.finditer(r"\S+", line))
        s, e = tok[2].span()          # 0 name, 1 daily, 2 monthly
        lines[i] = line[:s] + "y".ljust(e - s) + line[e:]
prt.write_text("\n".join(lines) + "\n")
print("print.prt: hru_wb and hru_nb monthly output enabled")
PY

cat > "$DST/time.sim" <<EOF
time.sim written by the netCDF fixture driver
day_start     yrc_start     day_end       yrc_end       step
0             $YR0          0             $YR1          0
EOF

cd "$DST"
"$EXE" 2>&1 | tee run.log

echo
echo "=== hru0001 annual precipitation (yr, precip mm) ==="
# hru_wb_yr.txt columns: 1 jday, 2 mon, 3 day, 4 yr, 5 unit, 6 gis_id, 7 name, 8 precip
awk 'NR>3 && $7=="hru0001" {printf "  %s  %s\n", $4, $8}' hru_wb_yr.txt

if [ "$DEP" != 1 ]; then
  echo
  echo "=== deposition disabled (NCDF_DEP=0) ==="
  if grep -qi "deposition" run.log; then
    echo "  FAIL: reader ran anyway"; exit 1
  fi
  echo "  PASS: no deposition reader output, run completed"
  exit 0
fi

echo
echo "=== hru0001 monthly atmospheric deposition (kg/ha) ==="
# Summing the daily nut_nrain expression over the N days of a month makes the
# dry term's /const cancel exactly, so the monthly total is
#     atmo = 0.01 * rf(k) * P_mon + dry(k)
# with k the running month index of the deposition record (k=1 is 1978-01).
# No leap-year dependence, and the dry term is a pure additive constant that
# steps 0.01 per month, so a frozen or offset atmodep_cont%ts is unmissable.
python3 - "$YR0" <<'PY'
import sys

HRU, DEP_EPOCH_YR = "hru0001", 1978
yr0 = int(sys.argv[1])


def rows(path, cols):
    out = {}
    for line in open(path).readlines()[3:]:
        f = line.split()
        if len(f) > max(cols) and f[6] == HRU:
            out[(int(f[3]), int(f[1]))] = [float(f[c]) for c in cols]
    return out


# hru_wb_mon.txt: 8th field (index 7) is precip
# hru_nb_mon.txt: 21st and 22nd fields (index 20, 21) are no3atmo, nh4atmo
precip = rows("hru_wb_mon.txt", [7])
atmo = rows("hru_nb_mon.txt", [20, 21])

hdr = ("  month      k   P_mon      no3atmo             nh4atmo\n"
       "                            obs      exp   diff    obs      exp   diff")
print(hdr)
worst = 0.0
for (yr, mo) in sorted(atmo):
    k = (yr - DEP_EPOCH_YR) * 12 + mo
    p = precip[(yr, mo)][0]
    obs_no3, obs_nh4 = atmo[(yr, mo)]
    exp_no3 = 0.01 * (0.20 + 0.01 * k) * p + (0.40 + 0.01 * k)
    exp_nh4 = 0.01 * (0.10 + 0.01 * k) * p + (0.30 + 0.01 * k)
    d3, d4 = obs_no3 - exp_no3, obs_nh4 - exp_nh4
    worst = max(worst, abs(d3), abs(d4))
    print(f"  {yr}-{mo:02d} {k:6d} {p:7.2f} "
          f"{obs_no3:8.3f} {exp_no3:8.3f} {d3:6.3f} "
          f"{obs_nh4:8.3f} {exp_nh4:8.3f} {d4:6.3f}")

# f17.3 output rounds at the third decimal, so anything under 0.001 is print
# resolution; a wrong month index would show up as a multiple of 0.01.
print(f"\nworst absolute difference: {worst:.4f}")
print("PASS" if worst < 1e-3 else "FAIL")
PY

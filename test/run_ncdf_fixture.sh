#!/usr/bin/env bash
# Run the SWAT+ netCDF reader against a synthetic fixture and report hru0001's
# annual precipitation, which the fixture makes analytically predictable.
#
# Usage: test/run_ncdf_fixture.sh <scratch-dir> [start-yr] [end-yr]
set -euo pipefail

REPO="$(cd "$(dirname "$0")/.." && pwd)"
SRC="$REPO/data/mississippi_downstream/TxtInOut"
DST="${1:?usage: run_ncdf_fixture.sh <scratch-dir> [start-yr] [end-yr]}"
YR0="${2:-1980}"
YR1="${3:-1981}"

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
rm -f "$DST"/20crv3_era5_1960_2021.nc

python3 "$REPO/test/make_ncdf_fixture.py" "$DST/fixture.nc"

# file.cio carries the netCDF path on its pcp_path line (and tmp/slr/hmd/wnd,
# which cli_ncdf_meas ignores -- only pcp_path is ever read)
sed -i 's#20crv3_era5_1960_2021\.nc#fixture.nc#' "$DST/file.cio"

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

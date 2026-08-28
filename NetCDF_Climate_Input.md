# NetCDF Climate Input Support for SWAT+

This document describes how to use NetCDF files as climate input for SWAT+ simulations. NetCDF climate input allows you to use gridded climate data (such as reanalysis or climate model output) directly in SWAT+ without converting to individual station files.

## Overview

Instead of using traditional text-based climate files (`.pcp`, `.tmp`, `.slr`, `.hmd`, `.wnd`), SWAT+ can read climate data directly from a NetCDF4 file. The model extracts data for each weather station by finding the nearest grid cell in the NetCDF file.

## Requirements

### Build Requirements

NetCDF support must be enabled at compile time:

```bash
cmake -DENABLE_NETCDF=ON -DNETCDF_ROOT=/path/to/netcdf ..
make
```

If SWAT+ was not built with NetCDF support, you will see an error message:
```
! Error: NetCDF support is not enabled in this build.
```

### NetCDF File Requirements


if you have a project you want to test with, you can convert automatically using https://github.com/celray/swat2netcdf

the C++ version is very recommended as it converts faster even though the Python Version is parallelized, check the release section of the repo or compile yourself.

run to check options:
```bash
swat2netcdf -h
```


The NetCDF file must contain:

**Dimensions:**
- `time` - Number of time steps (daily data)
- `lat` - Number of latitude points
- `lon` - Number of longitude points

**Coordinate Variables:**
- `lat` - Latitude values (degrees)
- `lon` - Longitude values (degrees)
- `time` - Time values with a `units` attribute of the form `"days since YYYY-MM-DD"`

> Only `days since` works. The reader ignores the unit token before `since`, so
> an `"hours since ..."` axis is read as if the numbers were whole days, and the
> dates come out wrong with no warning. Convert the axis to daily first.

**Climate Variables** (at least one required):
| Variable | Description | Units |
|----------|-------------|-------|
| `pcp` | Precipitation | mm/day |
| `tmin` | Minimum temperature | degrees C |
| `tmax` | Maximum temperature | degrees C |
| `slr` | Solar radiation | MJ/m2/day |
| `hmd` | Relative humidity | fraction (0-1) |
| `wnd` | Wind speed | m/s |

All six variables have to be in the same file. Only the `pcp_path` line of
`file.cio` is ever read; `tmp_path`, `slr_path`, `hmd_path` and `wnd_path` get
parsed and then thrown away.

Missing data goes to the weather generator. Three things count as "not
measured" for a given station-day:

- a variable absent from the file, which is named in a startup message
- a negative fill value such as `-9999`
- `NaN`, which is what most reanalysis products use for `_FillValue`

Days the time axis does not cover go to the generator too. That includes a
partial first or last year and any gap in the middle, so they read as missing
rather than as zero.

## Configuration Files

### 1. The `.ncw` File (Weather Station Definitions)

The `.ncw` file defines the weather stations and how to apply data from the NetCDF grid.

> Name the file exactly `netcdf.ncw`. Detection is a match on that string, not
> a check of the `.ncw` extension. Call it `mygrid.ncw` and SWAT+ quietly goes
> back to text climate input without telling you.

Format:

```
netcdf.ncw: written by [author]
name                 wgn        latitude     longitude     elevation        pcp       tmin       tmax        slr        hmd       wnd        pet
station1         wgn_name1         -17.25         29.25           811        1.0        1.0        1.0        1.0        1.0        1.0      null
station2         wgn_name2         -17.25         29.75          1156        1.0        1.0        1.0        1.0        1.0        1.0      null
```

**Column Descriptions:**

| Column | Description |
|--------|-------------|
| `name` | Station identifier (used internally) |
| `wgn` | Weather generator station name (from `weather-wgn.cli`) |
| `latitude` | Station latitude in decimal degrees |
| `longitude` | Station longitude in decimal degrees |
| `elevation` | Station elevation in meters |
| `pcp` | Precipitation scale factor |
| `tmin` | Minimum temperature scale factor |
| `tmax` | Maximum temperature scale factor |
| `slr` | Solar radiation scale factor |
| `hmd` | Humidity scale factor |
| `wnd` | Wind speed scale factor |
| `pet` | PET specification, read but never used in NetCDF mode (see below) |

**Scale Factor Values:**
- `1.0` - Use the NetCDF data as-is
- Other positive values - Multiply NetCDF data by this factor (e.g., `1.1` increases values by 10%)
- `-99` - Generate this variable with the weather generator and ignore the NetCDF data
- `null` - Variable not used or handled separately

Measured PET is not available in NetCDF mode. `cli_petmeas` never runs, so the
`pet` column gets read into the station record and then sits there unused. Leave
it `null`.

Atmospheric deposition still uses the traditional `atmodep.cli` file. Its station
names need to match the `name` column of `netcdf.ncw`.

### 2. The `file.cio` Configuration

Modify your `file.cio` to enable NetCDF climate input:

**Climate Line:**
Change the first entry in the `climate` line to point to your `.ncw` file:

```
climate           netcdf.ncw        weather-wgn.cli   null              null              null              null              null              null              null
```

**Path Lines:**
Add path lines specifying the NetCDF file:

```
pcp_path          climate_data.nc4
tmp_path          climate_data.nc4
slr_path          climate_data.nc4
hmd_path          climate_data.nc4
wnd_path          climate_data.nc4
```

Only `pcp_path` is read. The other four are parsed and discarded, so all six
climate variables have to be in the one file `pcp_path` names. Point the other
lines at that same path anyway, so nobody reading the file later assumes
otherwise.

## Complete Example

### Example `file.cio` (relevant sections):

```
file.cio: written by swat-netcdf converter
simulation        time.sim          print.prt         null              object.cnt        null
basin             codes.bsn         parameters.bsn
climate           netcdf.ncw        weather-wgn.cli   null              null              null              null              null              null              null
connect           hru.con           null              rout_unit.con     ...
...
pcp_path          save.nc4
tmp_path          save.nc4
slr_path          save.nc4
hmd_path          save.nc4
wnd_path          save.nc4
```

### Example `netcdf.ncw`:

```
netcdf.ncw: written by Celray James
name                 wgn        latitude     longitude     elevation        pcp       tmin       tmax        slr        hmd       wnd        pet
s17250s29250e   173s294e          -17.25         29.25           811        1.0        1.0        1.0        1.0        1.0        1.0      null
s17250s29750e   173s297e          -17.25         29.75          1156        1.0        1.0        1.0        1.0        1.0        1.0      null
s17750s29250e   176s294e          -17.75         29.25           929        1.0        1.0        1.0        1.0        1.0        1.0      null
```

## How It Works

1. **Detection**: SWAT+ enters NetCDF mode when the climate file named in `file.cio` is exactly `netcdf.ncw`

2. **Station Reading**: The `.ncw` file is read to get station locations and scale factors

3. **NetCDF Loading**: The NetCDF file specified in `pcp_path` is opened and all climate data is loaded into memory

4. **Spatial Matching**: For each station defined in `.ncw`, the model finds the closest grid cell in the NetCDF file using Euclidean distance in latitude/longitude space

5. **Date Placement**: The first and last values of `time` are decoded through the `units` attribute to get the record's calendar extent. Each record then goes on the date its own `time` value gives. Records earlier than the simulation start year are skipped. The simulation's start date does not shift the data

6. **Time Series Population**: Daily time series are extracted from the nearest grid cell for each station

7. **Missing Data**: An absent variable, a negative fill value and `NaN` all become SWAT+'s `-99` "not measured" sentinel, which sends that station-day to the weather generator

8. **Scaling**: Scale factors from the `.ncw` file are applied to each measured value. Missing values are never scaled

9. **Bounds Checking**: Measured values are clamped, precipitation, radiation and wind to >= 0 and humidity to 0-1. Temperature is not clamped

## Troubleshooting

### Common Errors

**"No NetCDF file specified in pcp path"**
- Ensure `pcp_path` is defined in `file.cio`

**"NetCDF file does not exist"**
- Check that the path in `pcp_path` is correct and the file exists

**"Cannot find 'time' dimension"**
- Your NetCDF file must have dimensions named exactly `time`, `lat`, and `lon`

**"NetCDF support is not enabled"**
- Rebuild SWAT+ with `-DENABLE_NETCDF=ON`

**Climate data looks shifted by a whole number of years**
- Fixed by the date placement change above. On an older build the reader ignored
  the file's start date and put record 1 on 1 January of the simulation start year

**Simulation runs with zero rainfall despite a "pcp will use wgn" message**
- Fixed by the missing value change above. On an older build the `-99` sentinel
  was clamped to `0`, so the generator never kicked in

**`forrtl: error (65): floating invalid` while reading climate**
- Fixed by the missing value change. On an older build a `NaN` fill value hit an
  ordered comparison under `-fpe0` and killed the run

### Tips

- The NetCDF file path in `file.cio` is relative to the simulation directory
- The climate variables have to be dimensioned `(time, lat, lon)`. The reader
  does not check, and a differently ordered file is read wrong with no warning
- Time units have to be `"days since YYYY-MM-DD"`, per the note under Requirements
- The startup log prints `netcdf climate record: <start> to <end>`. That is the
  extent the reader decoded from your time axis, and the fastest way to confirm
  it read the dates the way you meant

## Creating NetCDF Climate Files

NetCDF climate files can be created from various sources:

- **ERA5/ERA-Interim**: Download from Copernicus Climate Data Store
- **CHIRPS**: Precipitation data from Climate Hazards Group
- **Custom data**: Use Python (xarray, netCDF4) or R (ncdf4) to create NetCDF files

Example Python code to create a compatible NetCDF file:

```python
import xarray as xr
import numpy as np
import pandas as pd

# Create sample data
times = pd.date_range('2000-01-01', periods=365, freq='D')
lats = np.arange(-22, -17, 0.5)
lons = np.arange(29, 38, 0.5)

# Create dataset
ds = xr.Dataset({
    'pcp': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 10),
    'tmax': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 15 + 20),
    'tmin': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 10 + 10),
    'slr': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 20 + 5),
    'hmd': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 0.5 + 0.3),
    'wnd': (['time', 'lat', 'lon'], np.random.rand(365, len(lats), len(lons)) * 3 + 1),
},
coords={
    'time': times,
    'lat': lats,
    'lon': lons,
})

# Set time encoding
ds.time.encoding['units'] = 'days since 1970-01-01'

# Save to NetCDF4
ds.to_netcdf('climate_data.nc4', format='NETCDF4')
```

## References

- SWAT+ Documentation: https://swatplus.gitbook.io/
- NetCDF Documentation: https://www.unidata.ucar.edu/software/netcdf/
- CF Conventions: https://cfconventions.org/

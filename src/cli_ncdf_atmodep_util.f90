      module cli_ncdf_atmodep_util

      !! Helpers for the netCDF atmospheric deposition reader (src/cli_ncdf_read_atmodep.f90).

      implicit none
      private

      public :: classify_timestep, atmodep_start_index, nearest_cell

      contains

      pure function classify_timestep (day_offsets, n) result (ts)

      !! SWAT+ timestep code for a CF time axis given in days.
      !! Returns "aa" for a single record, "mo" for roughly monthly spacing and "yr" for roughly annual.
      !! Anything else is "??", which later should be rejected to avoid silently passing other timesteps.

        integer, intent(in) :: n
        real, intent(in) :: day_offsets(n)
        character(len=2) :: ts
        real :: span

        if (n <= 1) then
          ts = "aa"
          return
        end if

        span = (day_offsets(n) - day_offsets(1)) / real(n - 1)

        if (span >= 20. .and. span <= 40.) then
          ts = "mo"
        else if (span >= 350. .and. span <= 380.) then
          ts = "yr"
        else
          ts = "??"
        end if

      end function classify_timestep

      pure integer function atmodep_start_index (timestep, yr_init, mo_init, yr_start, mo_start)
        !! 1-based index into the deposition series for the first simulated
        !! month (or year). A result < 1 or > num means the simulation period
        !! falls outside the record; the caller warns and lets nut_nrain's
        !! own range check drop the deposition to zero.

        character(len=*), intent(in) :: timestep
        integer, intent(in) :: yr_init, mo_init, yr_start, mo_start

        select case (trim(timestep))
        case ("mo")
          atmodep_start_index = 12 * (yr_start - yr_init) + (mo_start - mo_init) + 1
        case ("yr")
          atmodep_start_index = yr_start - yr_init + 1
        case default
          atmodep_start_index = 1
        end select

      end function atmodep_start_index


      pure subroutine nearest_cell (lat_vals, nlat, lon_vals, nlon,tgt_lat, tgt_lon, ilat, ilon)

        !! Indices of the grid cell closest to a target point.

        integer, intent(in) :: nlat, nlon
        real, intent(in) :: lat_vals(nlat), lon_vals(nlon)
        real, intent(in) :: tgt_lat, tgt_lon
        integer, intent(out) :: ilat, ilon

        integer :: k
        real :: best, d

        ilat = 1
        best = abs(lat_vals(1) - tgt_lat)
        do k = 2, nlat
          d = abs(lat_vals(k) - tgt_lat)
          if (d < best) then
            best = d
            ilat = k
          end if
        end do

        ilon = 1
        best = abs(lon_vals(1) - tgt_lon)
        do k = 2, nlon
          d = abs(lon_vals(k) - tgt_lon)
          if (d < best) then
            best = d
            ilon = k
          end if
        end do

      end subroutine nearest_cell

      end module cli_ncdf_atmodep_util

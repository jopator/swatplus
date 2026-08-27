      module cli_ncdf_date

      !! Calendar helpers for the netCDF climate reader (src/cli_ncdf_meas.f90).
      !! No netCDF or SWAT+ dependency, so it compiles without -DENABLE_NETCDF.
      !! Do NOT wrap this file in #ifdef USE_NETCDF.
      !!
      !! Dates convert through a Julian Day Number (a running count of days from
      !! a fixed origin), so "N days since a reference date" and "day of year"
      !! are both plain arithmetic. Proleptic Gregorian calendar, matching the
      !! climate files in use. Algorithm: Fliegel & Van Flandern.

      implicit none
      private

      public :: is_leap_year, ymd_to_jdn, jdn_to_ymd
      public :: day_of_year, add_days_to_date, parse_time_units

      contains

      pure logical function is_leap_year (yr)
        integer, intent(in) :: yr

        is_leap_year = (mod(yr, 4) == 0 .and. mod(yr, 100) /= 0) .or. mod(yr, 400) == 0

      end function is_leap_year

      !! Julian Day Number of a Gregorian calendar date
      pure integer function ymd_to_jdn (yr, mo, dy)
        integer, intent(in) :: yr, mo, dy
        integer :: a, y, m

        a = (14 - mo) / 12
        y = yr + 4800 - a
        m = mo + 12 * a - 3
        ymd_to_jdn = dy + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045

      end function ymd_to_jdn

      !! Gregorian calendar date of a Julian Day Number
      pure subroutine jdn_to_ymd (jdn, yr, mo, dy)
        integer, intent(in) :: jdn
        integer, intent(out) :: yr, mo, dy
        integer :: a, b, c, d, e, m

        a = jdn + 32044
        b = (4 * a + 3) / 146097
        c = a - (146097 * b) / 4
        d = (4 * c + 3) / 1461
        e = c - (1461 * d) / 4
        m = (5 * e + 2) / 153
        dy = e - (153 * m + 2) / 5 + 1
        mo = m + 3 - 12 * (m / 10)
        yr = 100 * b + d - 4800 + m / 10

      end subroutine jdn_to_ymd

      !! day of year, 1..366 -- the first index of the SWAT+ %ts arrays
      pure integer function day_of_year (yr, mo, dy)
        integer, intent(in) :: yr, mo, dy

        day_of_year = ymd_to_jdn(yr, mo, dy) - ymd_to_jdn(yr, 1, 1) + 1

      end function day_of_year

      !! date ndays after yr0-mo0-dy0
      pure subroutine add_days_to_date (yr0, mo0, dy0, ndays, yr, mo, dy)
        integer, intent(in) :: yr0, mo0, dy0, ndays
        integer, intent(out) :: yr, mo, dy

        call jdn_to_ymd (ymd_to_jdn(yr0, mo0, dy0) + ndays, yr, mo, dy)

      end subroutine add_days_to_date

      !! Reference date out of a CF "units" attribute, e.g. "days since 1900-01-01"
      !! or "days since 1970-01-01 00:00:00". -> NOTE" Always days !!
      
      pure subroutine parse_time_units (units_str, ref_yr, ref_mo, ref_dy)
        character(len=*), intent(in) :: units_str
        integer, intent(out) :: ref_yr, ref_mo, ref_dy

        integer :: pos1, pos2, ios
        character(len=20) :: date_part

        ref_yr = 1970
        ref_mo = 1
        ref_dy = 1

        pos1 = index(units_str, "since")
        if (pos1 > 0) then
          pos1 = pos1 + 5

          !! skip whitespace after "since"
          do while (pos1 <= len(units_str))
            if (units_str(pos1:pos1) /= ' ') exit
            pos1 = pos1 + 1
          end do

          !! the date runs to the next blank, or to the end of the string
          pos2 = index(units_str(pos1:), ' ')
          if (pos2 == 0) then
            pos2 = len(units_str) + 1
          else
            pos2 = pos1 + pos2 - 1
          end if

          date_part = units_str(pos1:pos2-1)

          !! YYYY-MM-DD
          read (date_part(1:4), *, iostat=ios) ref_yr
          if (ios == 0 .and. len_trim(date_part) >= 7) then
            read (date_part(6:7), *, iostat=ios) ref_mo
          end if
          if (ios == 0 .and. len_trim(date_part) >= 10) then
            read (date_part(9:10), *, iostat=ios) ref_dy
          end if
        end if

        !! reject anything implausible and fall back to the CF-common default
        if (ref_yr < 1000 .or. ref_yr > 5000 .or. ref_mo < 1 .or. ref_mo > 12 &
            .or. ref_dy < 1 .or. ref_dy > 31) then
          ref_yr = 1970
          ref_mo = 1
          ref_dy = 1
        end if

      end subroutine parse_time_units

      end module cli_ncdf_date

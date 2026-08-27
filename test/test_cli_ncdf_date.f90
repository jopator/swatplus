program test_cli_ncdf_date
    use cli_ncdf_date
    implicit none

    integer :: nfail = 0
    integer :: jdn, yr, mo, dy, k

    !! Known values. A round trip stays self-consistent even if every date is
    !! offset by the same amount, so at least one absolute anchor is required.
    call check_int ("jdn 2000-01-01", ymd_to_jdn(2000, 1, 1), 2451545)
    call check_int ("days 1900->1960",                                        &
                    ymd_to_jdn(1960, 1, 1) - ymd_to_jdn(1900, 1, 1), 21914)

    !! The exact decode 20crv3_era5_1960_2021.nc needs: "days since 1900-01-01"
    !! with time(1) = 21914, which must come out as 1960-01-01.
    call add_days_to_date (1900, 1, 1, 21914, yr, mo, dy)
    call check_int ("era5 first yr", yr, 1960)
    call check_int ("era5 first mo", mo, 1)
    call check_int ("era5 first dy", dy, 1)

    !! Round trip every 97 days over ~210 years, crossing 1900 (not a leap
    !! year), 2000 (a leap year) and 2100 (not).
    do k = 0, 800
        jdn = ymd_to_jdn(1890, 1, 1) + 97 * k
        call jdn_to_ymd (jdn, yr, mo, dy)
        call check_int ("jdn round trip", ymd_to_jdn(yr, mo, dy), jdn)
    end do

    !! Year boundaries, which is what Task 4 relies on day_of_year for
    call check_int ("doy 1980-12-31", day_of_year(1980, 12, 31), 366)
    call check_int ("doy 1981-12-31", day_of_year(1981, 12, 31), 365)
    call check_int ("doy 1900-12-31", day_of_year(1900, 12, 31), 365)

    if (nfail > 0) then
        write (*,'(a,i0,a)') "FAILED: ", nfail, " assertion(s)"
        error stop 1
    end if
    write (*,*) "all cli_ncdf_date tests passed"

contains

    subroutine check_int (label, got, want)
        character(len=*), intent(in) :: label
        integer, intent(in) :: got, want
        if (got /= want) then
            write (*,'(a,a,a,i0,a,i0)') "  FAIL ", label, ": got ", got, " want ", want
            nfail = nfail + 1
        end if
    end subroutine check_int

end program test_cli_ncdf_date

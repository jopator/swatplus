program test_cli_ncdf_atmodep_util
    use cli_ncdf_atmodep_util
    implicit none

    integer :: nfail = 0
    integer :: ilat, ilon
    real :: monthly(24), yearly(5), single(1), odd(3)
    real :: lats(4), lons(3)

    !! A monthly axis in "days since": irregular 28-31 day steps.
    !! These are the real offsets of 1960-01-01 .. 1961-12-01.
    monthly = [   0.,  31.,  60.,  91., 121., 152., 182., 213., 244., 274.,   &
                305., 335., 366., 397., 425., 456., 486., 517., 547., 578.,   &
                609., 639., 670., 700. ]
    call check_str ("monthly axis", classify_timestep(monthly, 24), "mo")

    yearly = [ 0., 365., 730., 1095., 1461. ]
    call check_str ("yearly axis", classify_timestep(yearly, 5), "yr")

    single = [ 0. ]
    call check_str ("single record", classify_timestep(single, 1), "aa")

    !! Daily spacing is not supported and must be reported, not guessed at.
    odd = [ 0., 1., 2. ]
    call check_str ("daily axis rejected", classify_timestep(odd, 3), "??")

    !! Start index. Record starts 1960-01, simulation starts 1980-01 ->
    !! January 1980 is month 241 of the series.
    call check_int ("mo start 1980", atmodep_start_index("mo", 1960, 1, 1980, 1), 241)
    !! The case the upstream text formula gets wrong: a record that does not
    !! start in January. Record starts 1960-06, simulation starts 1960-06 -> 1.
    call check_int ("mo start mid-year", atmodep_start_index("mo", 1960, 6, 1960, 6), 1)
    call check_int ("mo start offset", atmodep_start_index("mo", 1960, 6, 1961, 1), 8)
    call check_int ("yr start", atmodep_start_index("yr", 1960, 1, 1980, 1), 21)
    call check_int ("aa start", atmodep_start_index("aa", 1960, 1, 1980, 1), 1)
    !! A simulation starting before the record yields an index < 1, which the
    !! caller must detect rather than clamp.
    call check_int ("mo start before record", atmodep_start_index("mo", 1980, 1, 1960, 1), -239)

    !! Nearest cell, on a descending latitude axis like ISIMIP's.
    lats = [ 40.25, 39.75, 39.25, 38.75 ]
    lons = [ -94.25, -93.75, -93.25 ]
    call nearest_cell (lats, 4, lons, 3, 39.70, -93.70, ilat, ilon)
    call check_int ("nearest lat", ilat, 2)
    call check_int ("nearest lon", ilon, 2)
    !! Exactly on a cell centre
    call nearest_cell (lats, 4, lons, 3, 38.75, -94.25, ilat, ilon)
    call check_int ("exact lat", ilat, 4)
    call check_int ("exact lon", ilon, 1)
    !! Outside the grid clamps to the edge cell rather than failing
    call nearest_cell (lats, 4, lons, 3, 99.0, 0.0, ilat, ilon)
    call check_int ("clamp lat", ilat, 1)
    call check_int ("clamp lon", ilon, 3)

    if (nfail > 0) then
        write (*,'(a,i0,a)') "FAILED: ", nfail, " assertion(s)"
        error stop 1
    end if
    write (*,*) "all cli_ncdf_atmodep_util tests passed"

contains

    subroutine check_int (label, got, want)
        character(len=*), intent(in) :: label
        integer, intent(in) :: got, want
        if (got /= want) then
            write (*,'(a,a,a,i0,a,i0)') "FAIL ", label, ": got ", got, " want ", want
            nfail = nfail + 1
        end if
    end subroutine check_int

    subroutine check_str (label, got, want)
        character(len=*), intent(in) :: label, got, want
        if (trim(got) /= trim(want)) then
            write (*,'(a,a,a,a,a,a)') "FAIL ", label, ": got '", trim(got), "' want '", trim(want)
            nfail = nfail + 1
        end if
    end subroutine check_str

end program test_cli_ncdf_atmodep_util

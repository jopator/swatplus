#ifdef USE_NETCDF
subroutine cli_ncdf_meas
    
    ! This subroutine reads NetCDF climate data for SWAT+ simulation
    ! Using NetCDF C interface instead of Fortran interface

    use climate_module
    use maximum_data_module
    use time_module
    use input_file_module
    use iso_c_binding
    use cli_ncdf_date, only: day_of_year, add_days_to_date, parse_time_units
    use ieee_arithmetic, only: ieee_is_nan

    implicit none
    
    ! NetCDF C constants
    integer(c_int), parameter :: NC_NOERR = 0
    integer(c_int), parameter :: NC_NOWRITE = 0
    integer(c_int), parameter :: NC_FLOAT = 5
    integer(c_int), parameter :: NC_CHAR = 2
    integer(c_int), parameter :: NC_MAX_NAME = 256
    integer(c_int), parameter :: NC_MAX_VAR_DIMS = 32

    ! Sentinel used to disable a clamp in nc_value
    real, parameter :: NO_LIMIT = 1.e30

    ! NetCDF C function interfaces
    interface
        function nc_open_c(path, mode, ncidp) bind(c, name='nc_open')
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int), value, intent(in) :: mode
            integer(c_int), intent(out) :: ncidp
            integer(c_int) :: nc_open_c
        end function
        
        function nc_close_c(ncid) bind(c, name='nc_close')
            import :: c_int
            integer(c_int), value, intent(in) :: ncid
            integer(c_int) :: nc_close_c
        end function
        
        function nc_inq_dimid_c(ncid, name, dimidp) bind(c, name='nc_inq_dimid')
            import :: c_int, c_char
            integer(c_int), value, intent(in) :: ncid
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: dimidp
            integer(c_int) :: nc_inq_dimid_c
        end function
        
        function nc_inq_dimlen_c(ncid, dimid, lenp) bind(c, name='nc_inq_dimlen')
            import :: c_int, c_size_t
            integer(c_int), value, intent(in) :: ncid
            integer(c_int), value, intent(in) :: dimid
            integer(c_size_t), intent(out) :: lenp
            integer(c_int) :: nc_inq_dimlen_c
        end function
        
        function nc_inq_varid_c(ncid, name, varidp) bind(c, name='nc_inq_varid')
            import :: c_int, c_char
            integer(c_int), value, intent(in) :: ncid
            character(kind=c_char), intent(in) :: name(*)
            integer(c_int), intent(out) :: varidp
            integer(c_int) :: nc_inq_varid_c
        end function
        
        function nc_get_var_float_c(ncid, varid, ip) bind(c, name='nc_get_var_float')
            import :: c_int, c_float
            integer(c_int), value, intent(in) :: ncid
            integer(c_int), value, intent(in) :: varid
            real(c_float), intent(out) :: ip(*)
            integer(c_int) :: nc_get_var_float_c
        end function
        
        function nc_get_att_text_c(ncid, varid, name, ip) bind(c, name='nc_get_att_text')
            import :: c_int, c_char
            integer(c_int), value, intent(in) :: ncid
            integer(c_int), value, intent(in) :: varid
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(out) :: ip(*)
            integer(c_int) :: nc_get_att_text_c
        end function
        
        function nc_inq_attlen_c(ncid, varid, name, lenp) bind(c, name='nc_inq_attlen')
            import :: c_int, c_char, c_size_t
            integer(c_int), value, intent(in) :: ncid
            integer(c_int), value, intent(in) :: varid
            character(kind=c_char), intent(in) :: name(*)
            integer(c_size_t), intent(out) :: lenp
            integer(c_int) :: nc_inq_attlen_c
        end function
        
        function nc_strerror_c(ncerr) bind(c, name='nc_strerror')
            import :: c_int, c_ptr
            integer(c_int), value, intent(in) :: ncerr
            type(c_ptr) :: nc_strerror_c
        end function
    end interface
    
    ! NetCDF variables
    integer(c_int) :: ncid, varid, dimid
    integer(c_int) :: status
    character(len=256) :: ncdf_file
    character(len=257, kind=c_char) :: ncdf_file_c
    
    ! NetCDF dimensions
    integer(c_int) :: time_dimid, lat_dimid, lon_dimid
    integer(c_size_t) :: ntime_c, nlat_c, nlon_c
    integer :: ntime, nlat, nlon
    
    ! NetCDF variable IDs
    integer(c_int) :: time_varid, lat_varid, lon_varid
    integer(c_int) :: pcp_varid, tmin_varid, tmax_varid, slr_varid, hmd_varid, wnd_varid
    
    ! Arrays for reading ALL data
    real(c_float), dimension(:), allocatable :: time_vals, lat_vals, lon_vals
    real(c_float), dimension(:,:,:), allocatable :: pcp_data, tmin_data, tmax_data
    real(c_float), dimension(:,:,:), allocatable :: slr_data, hmd_data, wnd_data
    
    ! Variables for finding closest grid point
    integer :: target_lat_idx, target_lon_idx
    real :: min_dist, dist
    integer :: ilat, ilon
    
    ! Variables for date calculation
    integer :: days_since_ref
    ! Calendar extent of the netCDF file itself, decoded from its time axis
    integer :: nc_start_yr, nc_start_mo, nc_start_dy
    integer :: nc_end_yr, nc_end_mo, nc_end_dy
    ! first_yr = calendar year stored in row 1 of the %ts arrays
    ! nbyr_nc  = number of rows those arrays need
    integer :: first_yr, nbyr_nc
    integer :: ref_year, ref_month, ref_day
    character(len=256) :: time_units
    character(len=257, kind=c_char) :: time_units_c
    integer(c_size_t) :: att_len
    
    ! Loop counters and diagnostics
    integer :: itime, iyear, iday, i, iwst
    logical :: exists
    
    ! Helper function to convert C string pointer to Fortran string
    interface
        function c_strlen(str) bind(c, name='strlen')
            import :: c_ptr, c_size_t
            type(c_ptr), value, intent(in) :: str
            integer(c_size_t) :: c_strlen
        end function
    end interface
    
    write (*,*) "reading data using netcdf C interface"
    write (9003,*) "reading data using netcdf C interface"

    ! Get NetCDF filename based on precipitation path from file.cio
    if (in_path_pcp%pcp == "null" .or. trim(in_path_pcp%pcp) == " ") then
        write(*,*) "! error: No NetCDF file specified in pcp path in 'file.cio'"
        write (9003,*) "! error: No NetCDF file specified in pcp path in 'file.cio'"
        stop
    else
        ncdf_file = TRIM(ADJUSTL(in_path_pcp%pcp))
        inquire(file=trim(ncdf_file), exist=exists)
        if (.not. exists) then
            write(*,*) "! error: NetCDF file does not exist at ", trim(ncdf_file)
            write (9003,*) "! error: NetCDF file does not exist at ", trim(ncdf_file)
            stop
        end if
    endif
    
    ! Convert filename to C string
    call f_to_c_string(ncdf_file, ncdf_file_c)
    
    ! Open NetCDF file
    status = nc_open_c(ncdf_file_c, NC_NOWRITE, ncid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot open NetCDF file: ", trim(ncdf_file)
        write (*,*) "! NetCDF Error code: ", status
        write (9003,*) "! error: Cannot open NetCDF file: ", trim(ncdf_file)
        stop
    endif
    
    ! Get dimensions for gridded data
    status = nc_inq_dimid_c(ncid, "time" // c_null_char, time_dimid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot find 'time' dimension, code:", status
        write (9003,*) "! error: Cannot find 'time' dimension, code:", status
        stop
    endif
    
    status = nc_inq_dimid_c(ncid, "lat" // c_null_char, lat_dimid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot find 'lat' dimension, code:", status
        write (9003,*) "! error: Cannot find 'lat' dimension, code:", status
        stop
    endif
    
    status = nc_inq_dimid_c(ncid, "lon" // c_null_char, lon_dimid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot find 'lon' dimension, code:", status
        write (9003,*) "! error: Cannot find 'lon' dimension, code:", status
        stop
    endif
    
    ! Get dimension sizes
    status = nc_inq_dimlen_c(ncid, time_dimid, ntime_c)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot get time dimension length, code:", status
        stop
    endif
    ntime = int(ntime_c)
    
    status = nc_inq_dimlen_c(ncid, lat_dimid, nlat_c)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot get lat dimension length, code:", status
        stop
    endif
    nlat = int(nlat_c)
    
    status = nc_inq_dimlen_c(ncid, lon_dimid, nlon_c)
    if (status /= NC_NOERR) then
        write (*,*) "! error: Cannot get lon dimension length, code:", status
        stop
    endif
    nlon = int(nlon_c)
    
    ! Allocate arrays - storage order is (lon, lat, time)
    allocate(time_vals(ntime), lat_vals(nlat), lon_vals(nlon))
    allocate(pcp_data(nlon, nlat, ntime), tmin_data(nlon, nlat, ntime), tmax_data(nlon, nlat, ntime))
    allocate(slr_data(nlon, nlat, ntime), hmd_data(nlon, nlat, ntime), wnd_data(nlon, nlat, ntime))
    
    ! Read coordinate data
    status = nc_inq_varid_c(ncid, "lat" // c_null_char, lat_varid)
    if (status == NC_NOERR) then
        status = nc_get_var_float_c(ncid, lat_varid, lat_vals)
        if (status /= NC_NOERR) then
            write (*,*) "! error reading lat data, code: ", status
            write (9003,*) "! error reading lat data, code: ", status
            stop
        endif
    else
        write (*,*) "WARNING: Cannot find lat variable, code:", status
    endif
    
    status = nc_inq_varid_c(ncid, "lon" // c_null_char, lon_varid)
    if (status == NC_NOERR) then
        status = nc_get_var_float_c(ncid, lon_varid, lon_vals)
        if (status /= NC_NOERR) then
            write (*,*) "! error reading lon data, code: ", status
            write (9003,*) "! error reading lon data, code: ", status
            stop
        endif
    else
        write (*,*) "! WARNING: Cannot find lon variable, code:", status
    endif
    
    ! Read time data to get the date for first time step
    status = nc_inq_varid_c(ncid, "time" // c_null_char, time_varid)
    if (status == NC_NOERR) then
        ! Try to read time units attribute
        status = nc_inq_attlen_c(ncid, time_varid, "units" // c_null_char, att_len)
        if (status == NC_NOERR .and. att_len > 0) then
            ! Allocate space for the attribute (including null terminator)
            if (att_len < 256) then
                time_units_c = repeat(c_null_char, 257)
                status = nc_get_att_text_c(ncid, time_varid, "units" // c_null_char, time_units_c)
                if (status == NC_NOERR) then
                    call c_to_f_string(time_units_c, time_units)
                else
                    time_units = "days since 1970-01-01 00:00:00"
                    write (*,*) "! warning: Cannot read time units attribute, using default"
                    write (9003,*) "! warning: Cannot read time units attribute, using default"
                endif
            else
                time_units = "days since 1970-01-01 00:00:00"
                write (*,*) "! warning: Time units attribute too long, using default"
            endif
        else
            time_units = "days since 1970-01-01 00:00:00"
            write (*,*) "! warning: Cannot find time units attribute, using default"
            write (9003,*) "! warning: Cannot find time units attribute, using default"
        endif
        
        ! Parse reference date from time units string
        call parse_time_units(time_units, ref_year, ref_month, ref_day)
        
        status = nc_get_var_float_c(ncid, time_varid, time_vals)
        if (status /= NC_NOERR) then
            write (*,*) "! error reading time data, code: ", status
            write (9003,*) "! error reading time data, code: ", status
            stop
        endif

        if (ntime < 1) then
            write (*,*) "! error: NetCDF time dimension is empty"
            write (9003,*) "! error: NetCDF time dimension is empty"
            stop
        endif

        ! Decode the first and last records of the file's own time axis.

        days_since_ref = int(time_vals(1))
        call add_days_to_date(ref_year, ref_month, ref_day, days_since_ref,   &
                              nc_start_yr, nc_start_mo, nc_start_dy)
        days_since_ref = int(time_vals(ntime))
        call add_days_to_date(ref_year, ref_month, ref_day, days_since_ref,   &
                              nc_end_yr, nc_end_mo, nc_end_dy)

        ! Row 1 of the %ts arrays holds the first year at or after the
        ! simulation start; earlier records are skipped when populating.

        first_yr = max(nc_start_yr, time%yrc)
        nbyr_nc  = nc_end_yr - first_yr + 1

        write (*,'(a,i4.4,a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2)')               &
            " netcdf climate record: ", nc_start_yr, "-", nc_start_mo, "-",   &
            nc_start_dy, " to ", nc_end_yr, "-", nc_end_mo, "-", nc_end_dy
        write (9003,'(a,i4.4,a,i2.2,a,i2.2,a,i4.4,a,i2.2,a,i2.2)')            &
            " netcdf climate record: ", nc_start_yr, "-", nc_start_mo, "-",   &
            nc_start_dy, " to ", nc_end_yr, "-", nc_end_mo, "-", nc_end_dy

        if (nbyr_nc < 1) then
            write (*,*) "! warning: NetCDF climate record does not reach the simulation period;"
            write (*,*) "!          the weather generator will be used throughout"
            write (9003,*) "! warning: NetCDF climate record ends before the simulation starts"
            nbyr_nc = 1
        endif
    else
        write (*,*) "WARNING: No time variable found in NetCDF, code:", status
        write (9003,*) "WARNING: No time variable found in NetCDF, code:", status
        stop
    endif
    
    ! Read all climate variables
    call read_climate_variable("pcp", pcp_varid, pcp_data, "precipitation")
    call read_climate_variable("tmin", tmin_varid, tmin_data, "minimum temperature")
    call read_climate_variable("tmax", tmax_varid, tmax_data, "maximum temperature")
    call read_climate_variable("slr", slr_varid, slr_data, "solar radiation")
    call read_climate_variable("hmd", hmd_varid, hmd_data, "humidity")
    call read_climate_variable("wnd", wnd_varid, wnd_data, "wind speed")

    ! allocate and populate climate arrays
    allocate (pcp(0:db_mx%wst))
    allocate (pcp_n(db_mx%wst))
    allocate (tmp(0:db_mx%wst))
    allocate (tmp_n(db_mx%wst))  
    allocate (slr(0:db_mx%wst))
    allocate (slr_n(db_mx%wst))
    allocate (hmd(0:db_mx%wst))
    allocate (hmd_n(db_mx%wst))
    allocate (wnd(0:db_mx%wst))
    allocate (wnd_n(db_mx%wst))
    db_mx%pcpfiles = db_mx%wst
    db_mx%tmpfiles = db_mx%wst  
    db_mx%slrfiles = db_mx%wst
    db_mx%rhfiles = db_mx%wst
    db_mx%wndfiles = db_mx%wst
    
    ! Set the climate file indices for each station
    do i = 1, db_mx%wst
        wst(i)%wco%pgage = i
        wst(i)%wco%tgage = i
        wst(i)%wco%sgage = i
        wst(i)%wco%hgage = i
        wst(i)%wco%wgage = i
    end do
    
    ! Populate station metadata and time series data
    do iwst = 1, db_mx%wst
        
        ! Find closest grid point to this station
        min_dist = 999999.0
        target_lat_idx = 1
        target_lon_idx = 1
        
        do ilat = 1, nlat
            do ilon = 1, nlon
                dist = sqrt((lat_vals(ilat) - wst(iwst)%lat)**2 + (lon_vals(ilon) - wst(iwst)%lon)**2)
                if (dist < min_dist) then
                    min_dist = dist
                    target_lat_idx = ilat
                    target_lon_idx = ilon
                endif
            end do
        end do
        
        ! Set station names and metadata
        call setup_station_metadata(iwst)
        
        ! Setup time series arrays
        call setup_timeseries_arrays(iwst)
        
        ! Populate time series data
        call populate_timeseries_data(iwst, target_lat_idx, target_lon_idx, ntime)
        
    end do
    
    ! Close NetCDF file
    status = nc_close_c(ncid)
    if (status /= NC_NOERR) then
        write (*,*) "Warning: Error closing NetCDF file, code:", status
    endif
    
    ! Clean up allocated arrays
    deallocate(time_vals, lat_vals, lon_vals)
    deallocate(pcp_data, tmin_data, tmax_data, slr_data, hmd_data, wnd_data)
    
    write(*,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"
    write(9003,'(A,I0,A)') " successfully populated time series for ", db_mx%wst, " stations"

    return

contains

    ! Helper subroutine to convert Fortran string to C string
    subroutine f_to_c_string(f_str, c_str)
        character(len=*), intent(in) :: f_str
        character(len=*, kind=c_char), intent(out) :: c_str
        integer :: i, f_len
        
        f_len = len_trim(f_str)
        do i = 1, f_len
            c_str(i:i) = f_str(i:i)
        end do
        c_str(f_len+1:f_len+1) = c_null_char
        
    end subroutine f_to_c_string
    
    ! Helper subroutine to convert C string to Fortran string
    subroutine c_to_f_string(c_str, f_str)
        character(len=*, kind=c_char), intent(in) :: c_str
        character(len=*), intent(out) :: f_str
        integer :: i, c_len
        
        ! Find the null terminator
        c_len = 0
        do i = 1, len(c_str)
            if (c_str(i:i) == c_null_char) then
                c_len = i - 1
                exit
            endif
        end do
        
        if (c_len == 0) c_len = len(c_str)
        
        ! Copy characters
        do i = 1, min(c_len, len(f_str))
            f_str(i:i) = c_str(i:i)
        end do
        
        ! Pad with spaces if necessary
        if (c_len < len(f_str)) then
            f_str(c_len+1:) = ' '
        endif
        
    end subroutine c_to_f_string

    ! Helper subroutine to read a climate variable
    subroutine read_climate_variable(var_name, var_id, var_data, description)
        character(len=*), intent(in) :: var_name, description
        integer(c_int), intent(out) :: var_id
        real(c_float), dimension(:,:,:), intent(out) :: var_data
        
        status = nc_inq_varid_c(ncid, trim(var_name) // c_null_char, var_id)
        if (status == NC_NOERR) then
            status = nc_get_var_float_c(ncid, var_id, var_data)
            if (status /= NC_NOERR) then
                write (*,*) "! error reading ", description, ", code: ", status
                write (9003,*) "! error reading ", description, ", code: ", status
                stop
            endif
        else
            write (*,*) trim(var_name), " will use wgn"
            write (9003,*) "WARNING: ", trim(var_name), " was not found in NetCDF, wgn will be used, code:", status
            ! -99. is SWAT+'s "no measured data" sentinel. nc_value() carries it
            ! through to the %ts arrays unscaled and unclamped, which is what
            ! makes climate_control's "<= -97." test fire and the generator run.
            var_data = -99.0
        endif
        
    end subroutine read_climate_variable

    ! Convert one raw netCDF value into a SWAT+ time series value.
    !
    ! SWAT+ treats any value <= -97. as "not measured -> goes to the weather
    ! generator". Three things must therefore survive untouched all the way
    ! into the %ts arrays as -99.:
    !   * the -99. that read_climate_variable writes for a variable that is
    !     absent from the file altogether
    !   * a negative _FillValue, e.g. the -9999. used by some source files
    !   * NaN, which is the _FillValue of e.g. the 20crv3-era5 files in use
    !
    ! Scaling or clamping any of those turns a gap into an actual value: before
    ! this function existed, a missing pcp was clamped from -99. to 0. and the
    ! whole simulation ran with zero rainfall instead of a generated one.
    !
    ! NaN is tested first and with ieee_is_nan, never with an ordered
    ! comparison: under the repo's -fpe0 / -ffpe-trap=invalid flags a NaN in a
    ! "<" test raises the invalid exception and aborts the run. That is not
    ! hypothetical -- it is what the old line 600 clamp did to a NaN humidity.
    real function nc_value(raw, factor, clamp_low, clamp_high)
        real, intent(in) :: raw, factor, clamp_low, clamp_high

        if (ieee_is_nan(raw)) then
            nc_value = -99.
        else if (raw <= -97.) then
            nc_value = -99.
        else
            nc_value = raw * factor
            if (nc_value < clamp_low) nc_value = clamp_low
            if (nc_value > clamp_high) nc_value = clamp_high
        end if

    end function nc_value

    ! Helper subroutine to setup station metadata
    subroutine setup_station_metadata(iwst)
        integer, intent(in) :: iwst
        
        ! Set station names as "filenames"
        pcp_n(iwst) = wst(iwst)%name
        tmp_n(iwst) = wst(iwst)%name
        slr_n(iwst) = wst(iwst)%name 
        hmd_n(iwst) = wst(iwst)%name
        wnd_n(iwst) = wst(iwst)%name
        
        ! Populate station metadata - use station coordinates (not NetCDF grid)
        pcp(iwst)%filename = wst(iwst)%name
        pcp(iwst)%lat = wst(iwst)%lat
        pcp(iwst)%long = wst(iwst)%lon
        pcp(iwst)%elev = wst(iwst)%elev
        
        ! Copy metadata to other climate arrays
        tmp(iwst)%filename = pcp(iwst)%filename
        tmp(iwst)%lat = pcp(iwst)%lat
        tmp(iwst)%long = pcp(iwst)%long
        tmp(iwst)%elev = pcp(iwst)%elev
        
        slr(iwst)%filename = pcp(iwst)%filename
        slr(iwst)%lat = pcp(iwst)%lat
        slr(iwst)%long = pcp(iwst)%long
        slr(iwst)%elev = pcp(iwst)%elev
        
        hmd(iwst)%filename = pcp(iwst)%filename
        hmd(iwst)%lat = pcp(iwst)%lat
        hmd(iwst)%long = pcp(iwst)%long
        hmd(iwst)%elev = pcp(iwst)%elev
        
        wnd(iwst)%filename = pcp(iwst)%filename
        wnd(iwst)%lat = pcp(iwst)%lat
        wnd(iwst)%long = pcp(iwst)%long
        wnd(iwst)%elev = pcp(iwst)%elev
        
    end subroutine setup_station_metadata
    
    ! Helper subroutine to setup time series arrays
    !
    ! The bounds follow the convention cli_pmeas.f90 establishes and that
    ! climate_control.f90 / cli_precip_control.f90 rely on:
    !   %start_yr / %start_day / %end_yr / %end_day  describe the FILE's extent and feed cli_bounds_check
    !   %yrs_start  simulation years elapsed before the record begins
    !   row r of %ts(day, r) holds calendar year first_yr + r - 1, where
    !   first_yr = max(nc_start_yr, time%yrc), because the run-time row index is time%yrs - %yrs_start

    subroutine setup_timeseries_arrays(iwst)
        integer, intent(in) :: iwst

        pcp(iwst)%nbyr = nbyr_nc
        tmp(iwst)%nbyr = nbyr_nc
        slr(iwst)%nbyr = nbyr_nc
        hmd(iwst)%nbyr = nbyr_nc
        wnd(iwst)%nbyr = nbyr_nc

        ! Set timestep (0 = daily)
        pcp(iwst)%tstep = 0
        tmp(iwst)%tstep = 0
        slr(iwst)%tstep = 0
        hmd(iwst)%tstep = 0
        wnd(iwst)%tstep = 0

        ! Initialize counters
        pcp(iwst)%days_gen = 0
        tmp(iwst)%days_gen = 0
        slr(iwst)%days_gen = 0
        hmd(iwst)%days_gen = 0
        wnd(iwst)%days_gen = 0

        ! Bounds taken from the netCDF file's own time axis
        pcp(iwst)%start_yr  = nc_start_yr
        pcp(iwst)%start_day = day_of_year(nc_start_yr, nc_start_mo, nc_start_dy)
        pcp(iwst)%end_yr    = nc_end_yr
        pcp(iwst)%end_day   = day_of_year(nc_end_yr, nc_end_mo, nc_end_dy)
        pcp(iwst)%yrs_start = max(0, nc_start_yr - time%yrc)

        ! Copy to other climate variables
        tmp(iwst)%start_yr = pcp(iwst)%start_yr
        tmp(iwst)%end_yr = pcp(iwst)%end_yr
        tmp(iwst)%start_day = pcp(iwst)%start_day
        tmp(iwst)%end_day = pcp(iwst)%end_day
        tmp(iwst)%yrs_start = pcp(iwst)%yrs_start

        slr(iwst)%start_yr = pcp(iwst)%start_yr
        slr(iwst)%end_yr = pcp(iwst)%end_yr
        slr(iwst)%start_day = pcp(iwst)%start_day
        slr(iwst)%end_day = pcp(iwst)%end_day
        slr(iwst)%yrs_start = pcp(iwst)%yrs_start

        hmd(iwst)%start_yr = pcp(iwst)%start_yr
        hmd(iwst)%end_yr = pcp(iwst)%end_yr
        hmd(iwst)%start_day = pcp(iwst)%start_day
        hmd(iwst)%end_day = pcp(iwst)%end_day
        hmd(iwst)%yrs_start = pcp(iwst)%yrs_start

        wnd(iwst)%start_yr = pcp(iwst)%start_yr
        wnd(iwst)%end_yr = pcp(iwst)%end_yr
        wnd(iwst)%start_day = pcp(iwst)%start_day
        wnd(iwst)%end_day = pcp(iwst)%end_day
        wnd(iwst)%yrs_start = pcp(iwst)%yrs_start

        ! Allocate time series arrays.
        ! Seed with -99. (SWAT+ "missing, generate instead") rather than 0., so
        ! that days the netCDF file does not cover -- a partial first or last
        ! year, or a gap in the time axis -- fall back to the weather generator
        ! instead of silently reading as zero rain / zero radiation.
        allocate (pcp(iwst)%ts(366, nbyr_nc), source = -99.)
        allocate (tmp(iwst)%ts(366, nbyr_nc), source = -99.)   ! ts for TMAX
        allocate (tmp(iwst)%ts2(366, nbyr_nc), source = -99.)  ! ts2 for TMIN
        allocate (slr(iwst)%ts(366, nbyr_nc), source = -99.)
        allocate (hmd(iwst)%ts(366, nbyr_nc), source = -99.)
        allocate (wnd(iwst)%ts(366, nbyr_nc), source = -99.)

    end subroutine setup_timeseries_arrays
    
    ! Helper subroutine to populate time series data
    !
    ! Walks the netCDF time axis and places each record on the calendar date
    ! that axis specifies. Records before the simulation start year, or beyond the last row, are skipped.

    subroutine populate_timeseries_data(iwst, target_lat_idx, target_lon_idx, ntime_total)
        integer, intent(in) :: iwst, target_lat_idx, target_lon_idx, ntime_total

        integer :: nc_yr, nc_mo, nc_dy

        do itime = 1, ntime_total
            call add_days_to_date(ref_year, ref_month, ref_day,               &
                                  int(time_vals(itime)), nc_yr, nc_mo, nc_dy)

            iyear = nc_yr - first_yr + 1
            if (iyear < 1 .or. iyear > nbyr_nc) cycle
            iday = day_of_year(nc_yr, nc_mo, nc_dy)

            ! Every assignment goes through nc_value, which routes missing data
            ! (absent variable, negative fill, NaN) to -99. and scales and
            ! clamps everything else.

            ! Precipitation - apply station scaling factor
            pcp(iwst)%ts(iday, iyear) = nc_value(                             &
                pcp_data(target_lon_idx, target_lat_idx, itime),              &
                wst(iwst)%pcp_factor, 0., NO_LIMIT)

            ! Temperature - populate both ts (TMAX) and ts2 (TMIN), neither of
            ! which is clamped
            tmp(iwst)%ts(iday, iyear) = nc_value(                             &
                tmax_data(target_lon_idx, target_lat_idx, itime),             &
                wst(iwst)%tmax_factor, -NO_LIMIT, NO_LIMIT)

            tmp(iwst)%ts2(iday, iyear) = nc_value(                            &
                tmin_data(target_lon_idx, target_lat_idx, itime),             &
                wst(iwst)%tmin_factor, -NO_LIMIT, NO_LIMIT)

            ! Solar radiation
            slr(iwst)%ts(iday, iyear) = nc_value(                             &
                slr_data(target_lon_idx, target_lat_idx, itime),              &
                wst(iwst)%slr_factor, 0., NO_LIMIT)

            ! Humidity
            hmd(iwst)%ts(iday, iyear) = nc_value(                             &
                hmd_data(target_lon_idx, target_lat_idx, itime),              &
                wst(iwst)%hmd_factor, 0., 1.)

            ! Wind speed
            wnd(iwst)%ts(iday, iyear) = nc_value(                             &
                wnd_data(target_lon_idx, target_lat_idx, itime),              &
                wst(iwst)%wnd_factor, 0., NO_LIMIT)
        end do

    end subroutine populate_timeseries_data

end subroutine cli_ncdf_meas

#else

! Stub subroutine when NetCDF support is disabled
subroutine cli_ncdf_meas
    implicit none
    
    write(*,*) "! Error: NetCDF support is not enabled in this build."
    write(*,*) "       Please rebuild SWAT+ with -DENABLE_NETCDF=ON to use NetCDF climate inputs."
    write(*,*) "       Or use traditional climate input files instead of 'netcdf.ncw'."
    write(9003,*) "! Error: NetCDF support is not enabled in this build."
    stop "NetCDF support disabled"
    
end subroutine cli_ncdf_meas

#endif
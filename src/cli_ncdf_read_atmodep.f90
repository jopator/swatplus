#ifdef USE_NETCDF
subroutine cli_ncdf_read_atmodep

    ! Reads gridded the atmospheric deposition from the netCDF defined in file.cio.
    ! Fills the same atmodep(:) arrays that cli_read_atmodep fills from atmodep.cli.

    use climate_module
    use maximum_data_module
    use time_module
    use input_file_module
    use iso_c_binding
    use nc_c_bindings_module
    use cli_ncdf_date, only: add_days_to_date, parse_time_units
    use cli_ncdf_atmodep_util, only: classify_timestep, atmodep_start_index, nearest_cell
    use ieee_arithmetic, only: ieee_is_nan

    implicit none

    ! Sentinel in an .ncw deposition column. -99 means no deposition on that station.
    real, parameter :: DEP_SKIP = -99.

    integer(c_int) :: ncid, status
    integer(c_int) :: time_dimid, lat_dimid, lon_dimid
    integer(c_int) :: time_varid, lat_varid, lon_varid
    integer(c_size_t) :: ntime_c, nlat_c, nlon_c, att_len
    integer :: ntime, nlat, nlon

    character(len=256) :: dep_file, time_units
    character(len=257, kind=c_char) :: dep_file_c, time_units_c

    real(c_float), dimension(:), allocatable :: time_vals, lat_vals, lon_vals
    real(c_float), dimension(:,:,:), allocatable :: nh4_rf_d, no3_rf_d, nh4_dry_d, no3_dry_d

    integer :: ref_year, ref_month, ref_day
    integer :: yr_init, mo_init, dy_init
    integer :: iwst, it, ilat, ilon, start_idx
    logical :: exists

    write (*,*) "reading atmospheric deposition from netcdf"
    write (9003,*) "reading atmospheric deposition from netcdf"

    dep_file = trim(adjustl(in_cli%atmo_cli))
    inquire (file=trim(dep_file), exist=exists)
    if (.not. exists) then
        write (*,*) "! error: deposition netCDF does not exist at ", trim(dep_file)
        write (9003,*) "! error: deposition netCDF does not exist at ", trim(dep_file)
        stop
    end if

    call f_to_c_string(dep_file, dep_file_c)
    status = nc_open_c(dep_file_c, NC_NOWRITE, ncid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: cannot open deposition netCDF: ", trim(dep_file), " code:", status
        write (9003,*) "! error: cannot open deposition netCDF: ", trim(dep_file), " code:", status
        stop
    end if

    call get_dim("time", time_dimid, ntime_c)
    call get_dim("lat", lat_dimid, nlat_c)
    call get_dim("lon", lon_dimid, nlon_c)
    ntime = int(ntime_c)
    nlat  = int(nlat_c)
    nlon  = int(nlon_c)

    if (ntime < 1) then
        write (*,*) "! error: deposition netCDF has an empty time dimension"
        write (9003,*) "! error: deposition netCDF has an empty time dimension"
        stop
    end if

    allocate (time_vals(ntime), lat_vals(nlat), lon_vals(nlon))
    allocate (nh4_rf_d(nlon, nlat, ntime), no3_rf_d(nlon, nlat, ntime))
    allocate (nh4_dry_d(nlon, nlat, ntime), no3_dry_d(nlon, nlat, ntime))

    call get_coord("lat", lat_varid, lat_vals)
    call get_coord("lon", lon_varid, lon_vals)

    ! Time axis: units attribute gives the epoch, values give the offsets.
    status = nc_inq_varid_c(ncid, "time" // c_null_char, time_varid)
    if (status /= NC_NOERR) then
        write (*,*) "! error: no time variable in deposition netCDF, code:", status
        write (9003,*) "! error: no time variable in deposition netCDF, code:", status
        stop
    end if

    time_units = "days since 1970-01-01 00:00:00"
    status = nc_inq_attlen_c(ncid, time_varid, "units" // c_null_char, att_len)
    if (status == NC_NOERR .and. att_len > 0 .and. att_len < 256) then
        time_units_c = repeat(c_null_char, 257)
        status = nc_get_att_text_c(ncid, time_varid, "units" // c_null_char, time_units_c)
        if (status == NC_NOERR) then
            call c_to_f_string(time_units_c, time_units)
        else
            write (*,*) "! warning: cannot read deposition time units, assuming days since 1970-01-01"
            write (9003,*) "! warning: cannot read deposition time units, assuming days since 1970-01-01"
        end if
    else
        write (*,*) "! warning: no deposition time units attribute, assuming days since 1970-01-01"
        write (9003,*) "! warning: no deposition time units attribute, assuming days since 1970-01-01"
    end if

    call parse_time_units(time_units, ref_year, ref_month, ref_day)

    status = nc_get_var_float_c(ncid, time_varid, time_vals)
    if (status /= NC_NOERR) then
        write (*,*) "! error reading deposition time values, code:", status
        write (9003,*) "! error reading deposition time values, code:", status
        stop
    end if

    call add_days_to_date(ref_year, ref_month, ref_day, int(time_vals(1)), &
                          yr_init, mo_init, dy_init)

    atmodep_cont%timestep = classify_timestep(time_vals, ntime)
    if (atmodep_cont%timestep == "??") then
        write (*,*) "! error: deposition netCDF time axis is neither monthly nor yearly."
        write (*,*) "!        Only 'aa' (one record), 'mo' and 'yr' are supported."
        write (9003,*) "! error: unsupported deposition time axis spacing"
        stop
    end if

    atmodep_cont%num_sta = db_mx%wst
    atmodep_cont%num     = ntime
    atmodep_cont%yr_init = yr_init
    atmodep_cont%mo_init = mo_init

    ! cli_atmodep_time_control only sets ts = 1 when the record starts in the
    ! very month the simulation starts, which is false for any record that
    ! begins earlier. Compute the index here and hand it a running counter.

    start_idx = atmodep_start_index(atmodep_cont%timestep, yr_init, mo_init, &
                                    time%yrc_start, time%mo_start)
    atmodep_cont%ts    = start_idx
    atmodep_cont%first = 0

    write (*,'(a,a,a,i0,a,i4.4,a,i2.2)')                                      &
        " deposition record: timestep ", trim(atmodep_cont%timestep),         &
        ", ", ntime, " steps from ", yr_init, "-", mo_init
    write (9003,'(a,a,a,i0,a,i4.4,a,i2.2)')                                   &
        " deposition record: timestep ", trim(atmodep_cont%timestep),         &
        ", ", ntime, " steps from ", yr_init, "-", mo_init

    if (start_idx < 1 .or. start_idx > ntime) then
        write (*,*) "! warning: simulation period lies outside the deposition record;"
        write (*,*) "!          atmospheric deposition will be zero"
        write (9003,*) "! warning: simulation period lies outside the deposition record"
    end if

    call get_dep_var("nh4_rf",  nh4_rf_d)
    call get_dep_var("no3_rf",  no3_rf_d)
    call get_dep_var("nh4_dry", nh4_dry_d)
    call get_dep_var("no3_dry", no3_dry_d)

    ! cli_read_atmodep already allocated these as empty; replace them.
    if (allocated(atmodep)) deallocate (atmodep)
    if (allocated(atmo_n))  deallocate (atmo_n)
    allocate (atmodep(0:db_mx%wst))
    allocate (atmo_n(db_mx%wst))
    db_mx%atmodep = db_mx%wst

    do iwst = 1, db_mx%wst
        atmodep(iwst)%name = wst(iwst)%name
        atmo_n(iwst) = wst(iwst)%name

        ! Index linkage, mirroring how cli_ncdf_meas sets wco%pgage = i.
        wst(iwst)%wco%atmodep = iwst

        call nearest_cell(lat_vals, nlat, lon_vals, nlon,                     &
                          wst(iwst)%lat, wst(iwst)%lon, ilat, ilon)

        select case (trim(atmodep_cont%timestep))

        case ("aa")
            atmodep(iwst)%nh4_rf  = dep_value(nh4_rf_d(ilon, ilat, 1),  wst(iwst)%nh4_rf_factor)
            atmodep(iwst)%no3_rf  = dep_value(no3_rf_d(ilon, ilat, 1),  wst(iwst)%no3_rf_factor)
            atmodep(iwst)%nh4_dry = dep_value(nh4_dry_d(ilon, ilat, 1), wst(iwst)%nh4_dry_factor)
            atmodep(iwst)%no3_dry = dep_value(no3_dry_d(ilon, ilat, 1), wst(iwst)%no3_dry_factor)

        case ("mo")
            allocate (atmodep(iwst)%nh4_rfmo(ntime),  source = 0.)
            allocate (atmodep(iwst)%no3_rfmo(ntime),  source = 0.)
            allocate (atmodep(iwst)%nh4_drymo(ntime), source = 0.)
            allocate (atmodep(iwst)%no3_drymo(ntime), source = 0.)
            do it = 1, ntime
                atmodep(iwst)%nh4_rfmo(it)  = dep_value(nh4_rf_d(ilon, ilat, it),  wst(iwst)%nh4_rf_factor)
                atmodep(iwst)%no3_rfmo(it)  = dep_value(no3_rf_d(ilon, ilat, it),  wst(iwst)%no3_rf_factor)
                atmodep(iwst)%nh4_drymo(it) = dep_value(nh4_dry_d(ilon, ilat, it), wst(iwst)%nh4_dry_factor)
                atmodep(iwst)%no3_drymo(it) = dep_value(no3_dry_d(ilon, ilat, it), wst(iwst)%no3_dry_factor)
            end do

        case ("yr")
            allocate (atmodep(iwst)%nh4_rfyr(ntime),  source = 0.)
            allocate (atmodep(iwst)%no3_rfyr(ntime),  source = 0.)
            allocate (atmodep(iwst)%nh4_dryyr(ntime), source = 0.)
            allocate (atmodep(iwst)%no3_dryyr(ntime), source = 0.)
            do it = 1, ntime
                atmodep(iwst)%nh4_rfyr(it)  = dep_value(nh4_rf_d(ilon, ilat, it),  wst(iwst)%nh4_rf_factor)
                atmodep(iwst)%no3_rfyr(it)  = dep_value(no3_rf_d(ilon, ilat, it),  wst(iwst)%no3_rf_factor)
                atmodep(iwst)%nh4_dryyr(it) = dep_value(nh4_dry_d(ilon, ilat, it), wst(iwst)%nh4_dry_factor)
                atmodep(iwst)%no3_dryyr(it) = dep_value(no3_dry_d(ilon, ilat, it), wst(iwst)%no3_dry_factor)
            end do

        end select
    end do

    status = nc_close_c(ncid)
    if (status /= NC_NOERR) write (*,*) "warning: error closing deposition netCDF, code:", status

    deallocate (time_vals, lat_vals, lon_vals)
    deallocate (nh4_rf_d, no3_rf_d, nh4_dry_d, no3_dry_d)

    write (*,'(a,i0,a)') " atmospheric deposition populated for ", db_mx%wst, " stations"
    write (9003,'(a,i0,a)') " atmospheric deposition populated for ", db_mx%wst, " stations"

    return

contains

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

    subroutine c_to_f_string(c_str, f_str)
        character(len=*, kind=c_char), intent(in) :: c_str
        character(len=*), intent(out) :: f_str
        integer :: i, c_len

        c_len = 0
        do i = 1, len(c_str)
            if (c_str(i:i) == c_null_char) then
                c_len = i - 1
                exit
            end if
        end do
        if (c_len == 0) c_len = len(c_str)

        do i = 1, min(c_len, len(f_str))
            f_str(i:i) = c_str(i:i)
        end do
        if (c_len < len(f_str)) f_str(c_len+1:) = ' '

    end subroutine c_to_f_string

    subroutine get_dim(dim_name, dimid, dlen)
        character(len=*), intent(in) :: dim_name
        integer(c_int), intent(out) :: dimid
        integer(c_size_t), intent(out) :: dlen

        status = nc_inq_dimid_c(ncid, trim(dim_name) // c_null_char, dimid)
        if (status /= NC_NOERR) then
            write (*,*) "! error: deposition netCDF has no '", trim(dim_name), "' dimension, code:", status
            write (9003,*) "! error: deposition netCDF has no '", trim(dim_name), "' dimension, code:", status
            stop
        end if
        status = nc_inq_dimlen_c(ncid, dimid, dlen)
        if (status /= NC_NOERR) then
            write (*,*) "! error: cannot read '", trim(dim_name), "' length, code:", status
            stop
        end if

    end subroutine get_dim

    subroutine get_coord(var_name, varid, vals)
        character(len=*), intent(in) :: var_name
        integer(c_int), intent(out) :: varid
        real(c_float), dimension(:), intent(out) :: vals

        status = nc_inq_varid_c(ncid, trim(var_name) // c_null_char, varid)
        if (status /= NC_NOERR) then
            write (*,*) "! error: deposition netCDF has no '", trim(var_name), "' variable, code:", status
            write (9003,*) "! error: deposition netCDF has no '", trim(var_name), "' variable, code:", status
            stop
        end if
        status = nc_get_var_float_c(ncid, varid, vals)
        if (status /= NC_NOERR) then
            write (*,*) "! error reading '", trim(var_name), "', code:", status
            stop
        end if

    end subroutine get_coord

    ! A deposition variable must be present not like weather data that can go to the wgn.
    subroutine get_dep_var(var_name, var_data)
        character(len=*), intent(in) :: var_name
        real(c_float), dimension(:,:,:), intent(out) :: var_data
        integer(c_int) :: varid

        status = nc_inq_varid_c(ncid, trim(var_name) // c_null_char, varid)
        if (status /= NC_NOERR) then
            write (*,*) "! error: deposition netCDF has no '", trim(var_name), "' variable, code:", status
            write (9003,*) "! error: deposition netCDF has no '", trim(var_name), "' variable, code:", status
            stop
        end if
        status = nc_get_var_float_c(ncid, varid, var_data)
        if (status /= NC_NOERR) then
            write (*,*) "! error reading '", trim(var_name), "', code:", status
            stop
        end if

    end subroutine get_dep_var

    ! One raw netCDF value into a deposition value.
    ! -99 means no deposition. No wgn for atmo dep
    real function dep_value(raw, factor)
        real, intent(in) :: raw, factor

        if (factor <= DEP_SKIP + 1.) then
            dep_value = 0.
        else if (ieee_is_nan(raw)) then
            dep_value = 0.
        else if (raw < 0.) then
            dep_value = 0.
        else
            dep_value = raw * factor
        end if

    end function dep_value

end subroutine cli_ncdf_read_atmodep
#endif

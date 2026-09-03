#ifdef USE_NETCDF
      module recall_ncdf_module

      !! Reads point source org_min data from the single orgmin.nc netCDF file that recall_db.rec
      !! names for every record, instead of one .rec file per object.

      !! The file is one flat table. Each row carries the recall_db id and name of the
      !! point source it belongs to, its date as jday/mo/yr, and the 18 hyd_output
      !! values. Rows are sorted by (id, time), so every record sits in one contiguous
      !! span that a single pass over the id column finds.

      use iso_c_binding
      use nc_c_bindings_module
      use recall_ncdf_util, only : recall_ncdf_file

      implicit none
      private

      public :: recall_ncdf_load, recall_ncdf_free

      !! hyd_output variables
      integer, parameter :: NVAL = 18
      character(len=6), parameter :: VAL_NAME(NVAL) =                             &
        [character(len=6) :: "flo", "sed", "orgn", "sedp", "no3", "solp", "chla", "nh3", "no2", "cbod", "dox", &
        "sand", "silt", "clay", "sag", "lag","gravel", "tmp"]

      logical :: loaded = .false.
      integer :: nrow = 0, name_len = 0, tstep_len = 0

      integer(c_int), allocatable   :: id_v(:), jday_v(:), mo_v(:), yr_v(:)
      real(c_float), allocatable    :: val_v(:,:)
      character(len=:), allocatable :: name_buf, tstep_buf ! These are long strings with the timestep and len of each record

      contains

      subroutine recall_ncdf_load (irec)

      !! Fill recall(irec)%hd from the rows belonging to recall object irec.
        use hydrograph_module
        use recall_module
        use time_module

        integer, intent(in) :: irec
        integer :: i, i0, i1, iyrs, step
        character(len=13) :: nm
        character(len=4) :: tstep       ! not "ts": hydrograph_module exports a public array by that name

        if (.not. loaded) call open_and_read

        !! rows are sorted by id, so this record is one contiguous span
        i0 = 0
        i1 = -1
        do i = 1, nrow
          if (id_v(i) == irec) then
            if (i0 == 0) i0 = i
            i1 = i
          end if
        end do

        if (i0 == 0) then
          call fail ("recall record " // trim(recall_db(irec)%name) // " has no rows in " // trim(recall_ncdf_file), 0)
        end if

        !! ids come from a counter that shifts when point sources are added in a
        !! different order, so the name is what confirms we have the right record
        nm = cut_null (name_buf((i0-1)*name_len + 1 : i0*name_len))
        if (trim(nm) /= trim(recall_db(irec)%name)) then
          call fail ("recall id " // trim(recall_db(irec)%name) // " points at " // trim(nm) // " in " // trim(recall_ncdf_file), 0)
        end if

        tstep = cut_null (tstep_buf((i0-1)*tstep_len + 1 : i0*tstep_len)) !Get timestep from long string blob

        if (trim(tstep) /= trim(recall_db(irec)%org_min%tstep)) then
          call fail ("tstep for " // trim(nm) // ": recall_db.rec says " // trim(recall_db(irec)%org_min%tstep) // ", the netCDF says " // trim(tstep), 0)
        end if

        select case (trim(tstep))
          case ("day")
            allocate (recall(irec)%hd(366, time%nbyr))
          case ("mo")
            allocate (recall(irec)%hd(12, time%nbyr))
          case ("yr")
            allocate (recall(irec)%hd(1, time%nbyr))
          case default
            call fail ("recall record " // trim(nm) // " has tstep " // trim(tstep) // ", this is not recorded on orgmin.nc", 0)
        end select

        recall(irec)%start_yr = yr_v(i0)
        recall(irec)%end_yr = yr_v(i0)

        do i = i0, i1
          if (yr_v(i) < time%yrc) cycle
          if (yr_v(i) > time%yrc_end) exit

          iyrs = yr_v(i) - time%yrc + 1
          if (iyrs < 1 .or. iyrs > time%nbyr) cycle

          select case (trim(tstep))
            case ("day")
              step = jday_v(i)
            case ("mo")
              step = mo_v(i)
            case default
              step = 1
          end select

          call set_hyd (recall(irec)%hd(step, iyrs), i)
          recall(irec)%end_yr = yr_v(i)
        end do

      end subroutine recall_ncdf_load


      subroutine recall_ncdf_free

        if (allocated(id_v)) deallocate (id_v, jday_v, mo_v, yr_v)
        if (allocated(val_v)) deallocate (val_v)
        if (allocated(name_buf)) deallocate (name_buf)
        if (allocated(tstep_buf)) deallocate (tstep_buf)
        loaded = .false.

      end subroutine recall_ncdf_free


      subroutine open_and_read

      !! Read the file whole on the first record that needs it. At CoSWAT sizes this is
      !! tens of MB, against reopening the file once per point source.

        integer(c_int) :: ncid, status
        character(len=257, kind=c_char) :: path_c
        logical :: exists
        integer :: k

        inquire (file=trim(recall_ncdf_file), exist=exists)
        if (.not. exists) call fail ("no recall netCDF at " // trim(recall_ncdf_file), 0)

        call f_to_c (trim(recall_ncdf_file), path_c)
        status = nc_open_c (path_c, NC_NOWRITE, ncid)
        if (status /= NC_NOERR) call fail ("cannot open " // trim(recall_ncdf_file), int(status))

        nrow      = dim_len (ncid, "row")
        name_len  = dim_len (ncid, "name_strlen")
        tstep_len = dim_len (ncid, "tstep_strlen")

        if (nrow < 1) call fail (trim(recall_ncdf_file) // " has no rows", 0)

        allocate (id_v(nrow), jday_v(nrow), mo_v(nrow), yr_v(nrow))
        allocate (val_v(nrow, NVAL))
        allocate (character(len=nrow*name_len) :: name_buf)
        allocate (character(len=nrow*tstep_len) :: tstep_buf)

        call get_int (ncid, "id", id_v)
        call get_int (ncid, "jday", jday_v)
        call get_int (ncid, "mo", mo_v)
        call get_int (ncid, "yr", yr_v)
        call get_text (ncid, "name", name_buf)
        call get_text (ncid, "tstep", tstep_buf)

        do k = 1, NVAL
          call get_float (ncid, trim(VAL_NAME(k)), val_v(:,k))
        end do

        status = nc_close_c (ncid)
        if (status /= NC_NOERR) write (*,*) "warning: error closing recall netCDF, code:", status

        loaded = .true.

        write (*,'(a,i0,a,a)') " recall: read ", nrow, " rows from ", trim(recall_ncdf_file)
        write (9003,'(a,i0,a,a)') " recall: read ", nrow, " rows from ", trim(recall_ncdf_file)

      end subroutine open_and_read


      subroutine set_hyd (h, i)

      !! One table row into one hyd_output slot.

        use hydrograph_module, only : hyd_output

        type (hyd_output), intent(out) :: h
        integer, intent(in) :: i

        h%flo  = val_v(i,1)
        h%sed  = val_v(i,2)
        h%orgn = val_v(i,3)
        h%sedp = val_v(i,4)
        h%no3  = val_v(i,5)
        h%solp = val_v(i,6)
        h%chla = val_v(i,7)
        h%nh3  = val_v(i,8)
        h%no2  = val_v(i,9)
        h%cbod = val_v(i,10)
        h%dox  = val_v(i,11)
        h%san  = val_v(i,12)
        h%sil  = val_v(i,13)
        h%cla  = val_v(i,14)
        h%sag  = val_v(i,15)
        h%lag  = val_v(i,16)
        h%grv  = val_v(i,17)
        h%temp = val_v(i,18)

      end subroutine set_hyd


      function cut_null (s) result (out)

      !! netCDF pads char variables with nulls, and trim() leaves those in place.

        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: k

        out = s
        k = index (out, char(0))
        if (k > 0) out(k:) = " "

      end function cut_null


      integer function dim_len (ncid, nm)

        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: nm
        integer(c_int) :: dimid, status
        integer(c_size_t) :: dlen

        status = nc_inq_dimid_c (ncid, nm // c_null_char, dimid)
        if (status /= NC_NOERR) call fail ("no '" // nm // "' dimension in " // trim(recall_ncdf_file), int(status))
        status = nc_inq_dimlen_c (ncid, dimid, dlen)
        if (status /= NC_NOERR) call fail ("cannot read '" // nm // "' length", int(status))

        dim_len = int (dlen)

      end function dim_len


      integer(c_int) function var_id (ncid, nm)

        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: nm
        integer(c_int) :: status

        status = nc_inq_varid_c (ncid, nm // c_null_char, var_id)
        if (status /= NC_NOERR) call fail ("no '" // nm // "' variable in " // trim(recall_ncdf_file), int(status))

      end function var_id


      subroutine get_int (ncid, nm, buf)

        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: nm
        integer(c_int), intent(out) :: buf(*)
        integer(c_int) :: status

        status = nc_get_var_int_c (ncid, var_id(ncid, nm), buf)
        if (status /= NC_NOERR) call fail ("cannot read '" // nm // "'", int(status))

      end subroutine get_int


      subroutine get_float (ncid, nm, buf)

        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: nm
        real(c_float), intent(out) :: buf(*)
        integer(c_int) :: status

        status = nc_get_var_float_c (ncid, var_id(ncid, nm), buf)
        if (status /= NC_NOERR) call fail ("cannot read '" // nm // "'", int(status))

      end subroutine get_float


      subroutine get_text (ncid, nm, buf)

        integer(c_int), intent(in) :: ncid
        character(len=*), intent(in) :: nm
        character(len=*), intent(out) :: buf
        integer(c_int) :: status

        status = nc_get_var_text_c (ncid, var_id(ncid, nm), buf)
        if (status /= NC_NOERR) call fail ("cannot read '" // nm // "'", int(status))

      end subroutine get_text


      subroutine f_to_c (f_str, c_str)

        character(len=*), intent(in) :: f_str
        character(len=*, kind=c_char), intent(out) :: c_str
        integer :: i, n

        n = len_trim (f_str)
        do i = 1, n
          c_str(i:i) = f_str(i:i)
        end do
        c_str(n+1:n+1) = c_null_char

      end subroutine f_to_c


      subroutine fail (msg, status)

      !! A point source that quietly loads nothing is worse than a run that stops, so
      !! every problem here ends the run.

        character(len=*), intent(in) :: msg
        integer, intent(in) :: status

        write (*,*) "! error: ", msg
        write (9003,*) "! error: ", msg
        if (status /= 0) then
          write (*,*) "         netCDF code:", status
          write (9003,*) "         netCDF code:", status
        end if
        stop

      end subroutine fail

      end module recall_ncdf_module

#endif

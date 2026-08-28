      module nc_c_bindings_module

      !! bind(c) interfaces to libnetcdf's C API
      !! Originally established on cli_ncdf_meas
      !! Here to be used consistently on other routines such as
      !! The netcdf reader for atmospheric deposition
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      !! Declarations only, so this compiles and links without netCDF lib.

      use iso_c_binding
      implicit none
      public

      ! NetCDF C constants
      integer(c_int), parameter :: NC_NOERR = 0
      integer(c_int), parameter :: NC_NOWRITE = 0
      integer(c_int), parameter :: NC_FLOAT = 5
      integer(c_int), parameter :: NC_CHAR = 2
      integer(c_int), parameter :: NC_MAX_NAME = 256
      integer(c_int), parameter :: NC_MAX_VAR_DIMS = 32

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

          function c_strlen(str) bind(c, name='strlen')
              import :: c_ptr, c_size_t
              type(c_ptr), value, intent(in) :: str
              integer(c_size_t) :: c_strlen
          end function
      end interface

      end module nc_c_bindings_module

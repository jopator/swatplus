      module recall_ncdf_util

      !! This module has some small helpers for the netcdf recall reader
      !! It is a separate module because it needs to read the recall_db.rec file and see if it points to a netCDF file
      !! And if it would live in the netcdf reader subroutines/modules, it would need NETCDF ENABLED in the build all the time
      !! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      implicit none
      public

      logical :: recall_ncdf_active = .false.
      character(len=256) :: recall_ncdf_file = ""

      contains

      subroutine recall_ncdf_detect (org_min_name)
        !! Just checks if the org_min_name provided in recall_db.rec is a netCDF file.
        !! If it is, switches recall_ncdf_active to .true. and also stores the name of the nc file.
        character(len=*), intent(in) :: org_min_name
        character(len=:), allocatable :: nm
        integer :: n

        nm = trim(adjustl(org_min_name))
        n = len(nm)

        if (n > 3) then
          if (nm(n-2:n) == ".nc") then
            recall_ncdf_active = .true.
            recall_ncdf_file = nm
          end if
        end if

      end subroutine recall_ncdf_detect

      end module recall_ncdf_util

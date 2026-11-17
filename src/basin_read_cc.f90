       subroutine basin_read_cc
      
       use input_file_module
       use basin_module
      
       implicit none
       
       character (len=80) :: titldum = "" !             |title of file
       character (len=80) :: header = ""  !             |header
       integer :: eof = 0               !             |end of file
       integer :: ios = 0               !             |read status with optional lake_evap column
       character (len=1000) :: line = ""!             |record of codes.bsn
       character (len=1010) :: line2 = ""!            |record of codes.bsn with lake_evap = 0 appended
       logical :: i_exist               !             |check to determine if file exists
       
       eof = 0
      
       !! read basin
       inquire (file=in_basin%codes_bas, exist=i_exist)
       if (i_exist .or. in_basin%codes_bas /= "null") then      
       do 
         open (107,file=in_basin%codes_bas)
         read (107,*,iostat=eof) titldum
         if (eof < 0) exit
         read (107,*,iostat=eof) header
         if (eof < 0) exit
         read (107,'(a)',iostat=eof) line
         if (eof < 0 .and. len_trim(line) == 0) exit
         !! idc_till and lake_evap (last two columns) are optional - missing ones get their defaults (idc_till = 3, lake_evap = 0)
         read (line,*,iostat=ios) bsn_cc
         if (ios /= 0) then
           line2 = trim(line) // " 0"
           read (line2,*,iostat=ios) bsn_cc
         end if
         if (ios /= 0) then
           line2 = trim(line) // " 3 0"
           read (line2,*,iostat=eof) bsn_cc
         end if
         exit
       enddo
       endif
       
       if (bsn_cc%pet == 3) then 
        open (140,file = 'pet.cli')
       do
        read (140,*,iostat=eof) titldum
        if (eof < 0) exit
        read (140,*,iostat=eof) header
        if (eof < 0) exit
        read (140,*,iostat = eof) titldum
        exit
       end do
       end if

       close(107)
       
       return
      end subroutine basin_read_cc
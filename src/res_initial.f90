      subroutine res_initial

      use basin_module
      use reservoir_module
      use maximum_data_module
      use reservoir_data_module
      use hydrograph_module
      use constituent_mass_module
      use pesticide_data_module
      use water_body_module
      use res_salt_module
      use res_cs_module
      
      implicit none
      
      integer :: ires = 0    !none          |counter
      real :: lnvol = 0.     !              |
      real :: resdif = 0.    !              |
      real :: stage = 0.     !m             |water level above lake bottom (GLOBathy)
      real :: dm = 0.        !m             |mean depth for air2water parameters
      logical :: a2w_exist   !              |air2water.res exists
      integer :: eof = 0     !              |end of file
      integer :: iu = 0      !              |unit of air2water.res
      character (len=80) :: titldum = ""      !  |title/header line of air2water.res
      character (len=40) :: a2w_name = ""     !  |reservoir (hydrology.res) name in air2water.res
      real, dimension(6) :: a2w_par = 0.      !  |air2water a1..a6 read from air2water.res
      integer :: i = 0       !none          |counter
      integer :: idat = 0    !none          |counter
      integer :: icon = 0    !none          |
      integer :: init = 0    !              | 
      integer :: ipest = 0   !none          |counter
      integer :: ipath = 0   !              |
      integer :: isalt = 0   !              |counter for salt ions !rtb salt
      integer :: ipest_db = 0  !none      |counter

      do ires = 1, sp_ob%res
        !! set initial volumes for res and hru types and convert units
        res_ob(ires)%ob  = sp_ob1%res + ires - 1
        res_ob(ires)%evol = res_hyd(ires)%evol * 10000.       !! ha-m => m**3
        res_ob(ires)%pvol = res_hyd(ires)%pvol * 10000.       !! ha-m => m**3
        res_ob(ires)%esa = res_hyd(ires)%esa
        res_ob(ires)%psa = res_hyd(ires)%psa
        !! set initial weir height to principal depth - m
        res_ob(ires)%weir_hgt = res_ob(ires)%pvol / (res_ob(ires)%psa * 10000.)

        !! set initial values for Hanazaki parameters   |Jose T 2025
        !!! Non-irrigation
        res_ob(ires)%I_mean               = 0.0
        res_ob(ires)%S_ini                = res_ob(ires)%pvol                                   !Hanasaki capacity = principal volume
        res_ob(ires)%N_memory             = 5                                                   !For now it will be fixed to 5
        allocate(res_ob(ires)%I_mon_past(12*res_ob(ires)%N_memory))
        res_ob(ires)%I_mon_past           = 0.0                                                 !At the beggining of the simulation this is a 1D array of N*12 zeros
        allocate(res_ob(ires)%daily_inflow_array(1))
        res_ob(ires)%daily_inflow_array   = 0.0                                                 !Same as above, 1d array with zeros
        res_ob(ires)%c_ratio              = 0.51

        !!! Irrigation
        res_ob(ires)%d_mean               = 0.0
        allocate(res_ob(ires)%d_mon_past(12*res_ob(ires)%N_memory))                             !Same logic as for inflow arrays
        res_ob(ires)%d_mon_past           = 0.0
        allocate(res_ob(ires)%daily_demand_array(1))
        res_ob(ires)%daily_demand_array   = 0.0
        res_ob(ires)%d_irrig_day          = 0.0
        res_ob(ires)%irrig_track          = 0

        !! use br1 as lag - then compute actual br1 (no option to input actual br1)
        res_ob(ires)%lag_up   = res_hyd(ires)%br1
        res_ob(ires)%lag_down = res_hyd(ires)%br2

        !! calculate shape parameters for surface area equation
        resdif = res_hyd(ires)%evol - res_hyd(ires)%pvol

        !! GLOBathy: A = hA_a*h**hA_b (km2), V = hV_a*h**hV_b (km3)  ->  area_ha = br1 * vol_m3**br2
        !! shape (br2) from GLOBathy, level anchored so that the area at the principal volume equals psa (a_scale)
        if (res_hyd(ires)%hA_a > 0. .and. res_hyd(ires)%hV_a > 0. .and. res_hyd(ires)%hV_b > 0.) then
          res_ob(ires)%br2 = res_hyd(ires)%hA_b / res_hyd(ires)%hV_b
          res_ob(ires)%br1 = 100. * res_hyd(ires)%hA_a * (1.e-9 / res_hyd(ires)%hV_a) ** res_ob(ires)%br2
          res_ob(ires)%a_scale = res_ob(ires)%psa / (res_ob(ires)%br1 * res_ob(ires)%pvol ** res_ob(ires)%br2)
          res_ob(ires)%br1 = res_ob(ires)%br1 * res_ob(ires)%a_scale
        !! For now it is normal procedure
        else if ((res_hyd(ires)%esa - res_hyd(ires)%psa) > 0. .and. resdif > 0.) then
          lnvol = Log10(res_ob(ires)%evol) - Log10(res_ob(ires)%pvol)
          if (lnvol > 1.e-4) then
            res_ob(ires)%br2 = (Log10(res_ob(ires)%esa) - Log10(res_ob(ires)%psa)) / lnvol
          else
            res_ob(ires)%br2 = (Log10(res_ob(ires)%esa) - Log10(res_ob(ires)%psa)) / 0.001
          end if
          if (res_ob(ires)%br2 > 0.9) then
            res_ob(ires)%br2 = 0.9
            res_ob(ires)%br1 = res_ob(ires)%psa / res_ob(ires)%pvol ** 0.9
          else
            res_ob(ires)%br1 = res_ob(ires)%esa / res_ob(ires)%evol ** res_ob(ires)%br2
          end if
        else
          res_ob(ires)%br2 = 0.9
          if (res_ob(ires)%pvol > 1.e-6) then
            res_ob(ires)%br1 = res_ob(ires)%psa / res_ob(ires)%pvol ** 0.9
          else
            res_ob(ires)%br1 = .1
          end if
        end if

        !! lake evaporation with Penman + heat storage and air2water Tw (bsn_cc%lake_evap = 1)    |Jose T 2026
        !! air2water a1..a6 from air2water.res (row with the hydrology.res name, not all zero); otherwise from the mean
        !! depth at the principal spillway with the regressions of Toffolon et al. (2014) (air2water pre_processing.m)
        if (bsn_cc%lake_evap == 1) then
          res_ob(ires)%dep_mean = max(res_ob(ires)%pvol / (res_ob(ires)%psa * 10000.), 0.1)
          dm = res_ob(ires)%dep_mean
          res_ob(ires)%a2w = 0.
          inquire (file="air2water.res", exist=a2w_exist)
          if (a2w_exist) then
            open (newunit=iu, file="air2water.res")
            read (iu,*,iostat=eof) titldum
            read (iu,*,iostat=eof) titldum
            do while (eof == 0)
              read (iu,*,iostat=eof) a2w_name, a2w_par
              if (eof /= 0) exit
              if (trim(a2w_name) == trim(res_hyd(ires)%name)) then
                res_ob(ires)%a2w = a2w_par
                exit
              end if
            end do
            close (iu)
            eof = 0
          end if
          if (all(res_ob(ires)%a2w == 0.)) then
            res_ob(ires)%a2w = (/ 0.487944 - 0.096334 * log(dm), 0.207095 * dm ** (-0.67172), 0.262359 * dm ** (-0.658578),   &
                                  31.331137 * dm ** (-0.329713), 0.843294 * dm ** (-0.731699), 0.627705 - 0.030032 * log(dm) /)
            write (9003,*) "air2water parameters from the depth regression for ", trim(res_hyd(ires)%name)
          end if
          res_ob(ires)%tw = 4.
          res_ob(ires)%th = 4.
          res_ob(ires)%ta_prev = -999.
          res_ob(ires)%tw_hist = 4.
          res_ob(ires)%tw_ihist = 1
        end if

      end do
      
      do ires = 1, sp_ob%res
        !! only initialize volume and constituents if operational at start of simulation
        !! if not, assume zero volume when dam is built
        if (time%yrc > res_hyd(ires)%iyres .or. (time%mo >= res_hyd(ires)%mores   &
                                   .and. time%yrc == res_hyd(ires)%iyres)) then
          idat = res_ob(ires)%props
          i = res_dat(idat)%init
        
          !! initialize org-min in reservoir
          init = res_init(i)%org_min
          res(ires) = om_init_water(init)
          call res_convert_mass (res(ires), res_ob(ires)%pvol)
        
          !! set initial reservoir org-min to reset for soft calibration
          res_om_init(ires) = res(ires)

          !! initialize pesticides in reservoir water and benthic from input data
          init = res_init(i)%pest
          do ipest = 1, cs_db%num_pests
            ipest_db = cs_db%pest_num(ipest)
            res_water(ires)%pest(ipest) = pest_water_ini(init)%water(ipest)
            res_benthic(ires)%pest(ipest) = pest_water_ini(init)%benthic(ipest)
            !! calculate mixing velocity using molecular weight and porosity
            res_ob(ires)%aq_mix(ipest) = pestdb(ipest_db)%mol_wt * (1. - res_sed(ires)%bd / 2.65)
          end do
                  
          !! initialize pathogens in reservoir water and benthic from input data
          init = res_init(i)%path
          do ipath = 1, cs_db%num_paths
            res_water(ires)%path(ipath) = path_water_ini(init)%water(ipath)
            res_benthic(ires)%path(ipath) = path_water_ini(init)%benthic(ipath)
          end do
                        
          !! calculate initial surface area - GLOBathy stage from volume, then area from stage; else br1/br2
          if (res_hyd(ires)%hA_a > 0. .and. res_hyd(ires)%hV_a > 0. .and. res_hyd(ires)%hV_b > 0.) then
            stage = (res(ires)%flo * 1.e-9 / res_hyd(ires)%hV_a) ** (1. / res_hyd(ires)%hV_b)
            res_wat_d(ires)%area_ha = res_ob(ires)%a_scale * 100. * res_hyd(ires)%hA_a * stage ** res_hyd(ires)%hA_b
          else
            res_wat_d(ires)%area_ha = res_ob(ires)%br1 * res(ires)%flo ** res_ob(ires)%br2
          end if

          !! initialize salts in reservoir water, from database file (salt.res)
          !rtb salt
          if(cs_db%num_salts > 0) then
            idat = res_ob(ires)%props
            icon = res_dat(idat)%salt !database to use (in salt_res file)
            if(icon > 0) then
              do isalt = 1, cs_db%num_salts
                res_water(ires)%saltc(isalt) = res_salt_data(icon)%c_init(isalt) !concentration (g/m3)
                res_water(ires)%salt(isalt) = res_water(ires)%saltc(isalt) * res(ires)%flo / 1000. !g/m3 * m3 / 1000. = kg
              enddo
            else
              do isalt = 1, cs_db%num_salts
                res_water(ires)%saltc(isalt) = 0. !concentration (g/m3)
                res_water(ires)%salt(isalt) = 0. !kg
              enddo 
            endif
          endif
          
          !! initialize constituents in reservoir water, from database file (cs_res)
          !rtb cs
          if(cs_db%num_cs > 0) then
            idat = res_ob(ires)%props
            icon = res_dat(idat)%cs !database to use (in cs.res file)
            if(icon > 0) then
              !seo4
              res_water(ires)%csc(1) = res_cs_data(icon)%c_seo4 !concentration (g/m3)
              res_water(ires)%cs(1) = res_water(ires)%csc(1) * res(ires)%flo / 1000. !g/m3 * m3 / 1000. = kg
              !seo3
              res_water(ires)%csc(2) = res_cs_data(icon)%c_seo3 !concentration (g/m3)
              res_water(ires)%cs(2) = res_water(ires)%csc(2) * res(ires)%flo / 1000. !g/m3 * m3 / 1000. = kg
              !boron
              res_water(ires)%csc(3) = res_cs_data(icon)%c_born !concentration (g/m3)
              res_water(ires)%cs(3) = res_water(ires)%csc(3) * res(ires)%flo / 1000. !g/m3 * m3 / 1000. = kg
            else
              !seo4
              res_water(ires)%csc(1) = 0.
              res_water(ires)%cs(1) = 0.
              !seo3
              res_water(ires)%csc(2) = 0.
              res_water(ires)%cs(2) = 0.
              !boron
              res_water(ires)%csc(3) = 0.
              res_water(ires)%cs(3) = 0.
            endif
          endif
          
        endif
        
      end do
      close(105)

      return
      end subroutine res_initial
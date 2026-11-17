      subroutine res_control (jres)
      
      use basin_module
      use reservoir_data_module 
      use time_module
      use reservoir_module
      use climate_module
      use hydrograph_module
      use conditional_module
      use water_body_module
      use constituent_mass_module  !! added nbs
      
      implicit none

      integer :: ii = 0               !none          |counter
      integer :: jres                 !none          |reservoir number
      integer :: idat = 0             !              |
      integer :: irel = 0             !              |
      integer :: iob = 0              !none          |counter
      integer :: ictbl = 0
      integer :: icon = 0              !! nbs
      real :: pvol_m3 = 0.
      real :: evol_m3 = 0.
      real :: dep = 0.
      real :: stage = 0.              !m             |water level above lake bottom (GLOBathy)
      real :: weir_hgt = 0.
      real :: alpha_up = 0.
      real :: alpha_down = 0.
      real, parameter :: pi2 = 6.2831853               !              |2 pi
      real :: ta_d = 0.                                !deg C         |mean air temperature of the day
      real :: t0 = 0.                                  !frac          |previous day as fraction of the year
      real :: t1 = 0.                                  !frac          |current day as fraction of the year
      real :: tw0 = 0.                                 !deg C         |lake surface temperature of the previous day
      real :: tw1 = 0.                                 !deg C         |lake surface temperature of the day
      real :: dd = 0.                                  !              |air2water normalized mixed layer depth (delta)
      real :: pp = 0.                                  !deg C/d       |air2water heat flux term of the previous day
      real :: tw30 = 0.                                !deg C         |lake surface temperature 30 days ago
      real :: u2 = 0.                                  !m/s           |wind speed at 2 m (SWAT+ windsp)
      real :: eas = 0.                                 !kPa           |saturation vapour pressure at air temperature
      real :: eair = 0.                                !kPa           |actual vapour pressure of the air
      real :: pres = 0.                                !kPa           |atmospheric pressure
      real :: xlam = 0.                                !MJ/kg         |latent heat of vaporisation
      real :: gma = 0.                                 !kPa/deg C     |psychrometric constant
      real :: dlt = 0.                                 !kPa/deg C     |slope of the saturation vapour pressure curve at ta
      real :: rnet = 0.                                !MJ/m2/d       |net radiation at the water surface
      real :: rso = 0.                                 !MJ/m2/d       |clear-sky solar radiation (no lwr only)
      real :: rnl = 0.                                 !MJ/m2/d       |FAO-56 net longwave radiation (no lwr only)
      real :: area_km2 = 0.                            !km2           |lake surface area
      real :: hmix = 0.                                !m             |heat storage depth: min(mean depth, wave base)
      real :: dheat = 0.                               !MJ/m2/d       |net energy gained from heat storage
      real :: tk = 0.                                  !K             |air temperature
      real :: visc = 0.                                !m2/s          |kinematic viscosity of air
      real :: ustar = 0.                               !m/s           |friction velocity
      real :: cu = 0.                                  !mm/d/(m/s mbar)|wind coefficient (Vardavas 1987)
      real :: ea_mm = 0.                               !mm/d          |evaporation component due to wind
      real :: evap_mm = 0.                             !mm/d          |lake evaporation
      integer :: it = 0                                !              |iteration counter

      integer :: dom            = 0                                           !           |Day of month
      integer :: mon            = 0                                           !           |Month of year
      integer :: end_of_mo      = 0                                           !           |End of month flag
      integer :: n_days         = 0                                           !           |Number of days in current month
      real    :: daily_inflow   = 0.                                          !m3         |Daily inflow from the past day
      real, dimension(:), allocatable :: temp_array                           !           |Temporary to store new values
      real :: daily_demand      = 0.                                          !m3         |Daily irrigation demand
      integer :: irrig_track_b  = 0                                           !none       |Tracker to update daily irrigation demand

      ht1 = ob(icmd)%hin    !! set incoming flow
      ht2 = resz            !! zero outgoing flow, sediment and nutrients
      hcs2 = hin_csz        !! zero outgoing constituents

      if (time%yrc > res_hyd(jres)%iyres .or. (time%mo >= res_hyd(jres)%mores   &
                                   .and. time%yrc == res_hyd(jres)%iyres)) then
        iob = res_ob(jres)%ob
        iwst = ob(iob)%wst
      
        !! adjust precip and temperature for elevation using lapse rates
        w = wst(iwst)%weat
        if (bsn_cc%lapse == 1) call cli_lapse
        wst(iwst)%weat = w
      
        !! set water body pointer to res
        wbody => res(jres)
        wbody_wb => res_wat_d(jres)
        wbody_prm => res_prm(jres)
      
        !! add incoming flow to reservoir
        res(jres) = res(jres) + ht1

        dom          = time%day_mo
        mon          = time%mo
        end_of_mo    = time%end_mo
        daily_inflow = ht1%flo

        if (irrig_track_b == res_ob(jres)%irrig_track) then
            daily_demand = 0

        else
            irrig_track_b = res_ob(jres)%irrig_track
            daily_demand = res_ob(jres)%d_irrig_day

        end if

        !! Store values in daily inflow array for the month and reset if month has ended
        if (dom == 1) then
            if (allocated(res_ob(jres)%daily_inflow_array)) then                                    ! Deallocate array to start over
                deallocate(res_ob(jres)%daily_inflow_array)
            end if

            if (allocated(res_ob(jres)%daily_demand_array)) then                                    ! Deallocate array to start over
                deallocate(res_ob(jres)%daily_demand_array)
            end if

            allocate(res_ob(jres)%daily_inflow_array(1))
            res_ob(jres)%daily_inflow_array(1) = daily_inflow                                       ! Store first inflow of the month

            allocate(res_ob(jres)%daily_demand_array(1))
            res_ob(jres)%daily_demand_array(1) = daily_demand                                       ! Store first irrigation demand of the month

        else
        !! Append inflow of current day
            n_days = size(res_ob(jres)%daily_inflow_array)

            allocate(temp_array(n_days+1))
            temp_array(1:n_days)   = res_ob(jres)%daily_inflow_array(1:n_days)
            temp_array(n_days +1)  = daily_inflow

            call move_alloc(temp_array, res_ob(jres)%daily_inflow_array)                            !Replace original array with move_alloc

        !! Append irrigaton demand of current day
            allocate(temp_array(n_days+1))
            temp_array(1:n_days)   = res_ob(jres)%daily_demand_array(1:n_days)
            temp_array(n_days +1)  = daily_demand

            call move_alloc(temp_array, res_ob(jres)%daily_demand_array)                            !Replace original array with move_alloc

        end if


        !! Get mean and store in reservoir's memory if end of month
        if (end_of_mo == 1) then
            ! Shift rolling window to the left
            res_ob(jres)%I_mon_past(1:12*(res_ob(jres)%N_memory)-1)     = res_ob(jres)%I_mon_past(2:12*(res_ob(jres)%N_memory))
            res_ob(jres)%I_mon_past(12*(res_ob(jres)%N_memory))         = sum(res_ob(jres)%daily_inflow_array) / real(size(res_ob(jres)%daily_inflow_array), kind=8)

            ! Do the same for irrigation
            res_ob(jres)%d_mon_past(1:12*(res_ob(jres)%N_memory)-1)     = res_ob(jres)%d_mon_past(2:12*(res_ob(jres)%N_memory))
            res_ob(jres)%d_mon_past(12*(res_ob(jres)%N_memory))         = sum(res_ob(jres)%daily_demand_array) / real(size(res_ob(jres)%daily_demand_array), kind=8)

        end if

        !! perform reservoir water/sediment balance
        idat = res_ob(jres)%props
        if(res_ob(jres)%rel_tbl == "d") then
          !! determine reservoir outflow
          irel = res_dat(idat)%release
          d_tbl => dtbl_res(irel)
          pvol_m3 = res_ob(jres)%pvol
          evol_m3 = res_ob(jres)%evol
          if (res_wat_d(jres)%area_ha > 1.e-6) then
            dep = wbody%flo / res_wat_d(jres)%area_ha / 10000.     !m = m3 / ha / 10000m2/ha
          else
            dep = 0.
          end if
          weir_hgt = res_ob(jres)%weir_hgt
          call conditions (jres, irel)

                  !! Retrospective information -> Inflow and irrigation demand memory of the reservoir

          call res_hydro (jres, irel, pvol_m3, evol_m3)
          
          !! new lag to smooth condition jumps (volume or month conditions)
          alpha_up = Exp(-res_ob(jres)%lag_up)
          alpha_down = Exp(-res_ob(jres)%lag_down)
          !! lag outflow when flows are receding
          if (res_ob(jres)%prev_flo < ht2%flo) then
            ht2%flo = ht2%flo * alpha_up + res_ob(jres)%prev_flo * (1. - alpha_up)
          else
            ht2%flo = ht2%flo * alpha_down + res_ob(jres)%prev_flo * (1. - alpha_down)
          end if
          res_ob(jres)%prev_flo = ht2%flo
            
          call res_sediment
        else
          ictbl = res_dat(idat)%release                              !! Osvaldo
          call res_rel_conds (ictbl, res(jres)%flo, ht1%flo, 0.)
          
        endif 
        
        !! calculate water balance for day
        if (bsn_cc%lake_evap == 1) then
          !! lake surface temperature tw: air2water, 6 parameters (Piccolroaz et al. 2013; Toffolon et al. 2014),
          !! daily Crank-Nicolson step as in github.com/marcotoffolon/air2water (t in fraction of year, tw >= 0)   |Jose T 2026
          ta_d = wst(iwst)%weat%tave
          if (res_ob(jres)%ta_prev < -900.) res_ob(jres)%ta_prev = ta_d
          t0 = float(time%day - 2) / float(time%day_end_yr)
          t1 = float(time%day - 1) / float(time%day_end_yr)
          tw0 = res_ob(jres)%tw
          if (tw0 >= res_ob(jres)%th) then
            dd = exp(-(tw0 - res_ob(jres)%th) / res_ob(jres)%a2w(4))
          else
            dd = 1.
          end if
          pp = res_ob(jres)%a2w(1) + res_ob(jres)%a2w(2) * res_ob(jres)%ta_prev - res_ob(jres)%a2w(3) * tw0          &
               + res_ob(jres)%a2w(5) * cos(pi2 * (t0 - res_ob(jres)%a2w(6)))
          tw1 = (2. * dd * tw0 + pp + res_ob(jres)%a2w(1) + res_ob(jres)%a2w(2) * ta_d                               &
               + res_ob(jres)%a2w(5) * cos(pi2 * (t1 - res_ob(jres)%a2w(6)))) / (2. * dd + res_ob(jres)%a2w(3))
          res_ob(jres)%tw = max(tw1, 0.)
          res_ob(jres)%th = max(4., min(res_ob(jres)%th, res_ob(jres)%tw))
          res_ob(jres)%ta_prev = ta_d

          !! tw of 30 days ago (oldest value of the ring buffer), then store today's tw in its place
          tw30 = res_ob(jres)%tw_hist(res_ob(jres)%tw_ihist)
          res_ob(jres)%tw_hist(res_ob(jres)%tw_ihist) = res_ob(jres)%tw
          res_ob(jres)%tw_ihist = mod(res_ob(jres)%tw_ihist, 30) + 1

          !! Penman with heat storage, McMahon et al. (2013) S10.7-S10.16 (Vardavas and Fountoulakis 1996), run daily:
          !!   E = D/(D+g) (Rn + dH)/lambda + g/(D+g) Ea
          !!   Rn = (1 - 0.08) Rs + 0.97 Rld - 0.97 sigma (tw + 273.15)^4
          !!   dH = -rho c h (tw - tw_30d) / 30,  h = min(mean depth, wave base 45.7 A^0.5 / (21.4 + A^0.5)) (Hakanson 2004)
          !!   Ea = Cu u2 (es(ta) - ea), Cu from roughness lengths (S10.12-S10.16)
          !! without incoming longwave (lwr): Rn = 0.92 Rs - Rnl (FAO-56 Eq. 39) - 4 sigma (ta + 273.15)^3 (tw - ta) (McMahon S11.7)
          !! no evaporation from an ice-covered lake (air2water tw = 0)
          if (res_ob(jres)%tw <= 0. .or. res_wat_d(jres)%area_ha <= 1.e-6) then
            evap_mm = 0.
          else
            eas = 0.5 * (0.6108 * exp(17.27 * wst(iwst)%weat%tmax / (wst(iwst)%weat%tmax + 237.3))                     &
                  + 0.6108 * exp(17.27 * wst(iwst)%weat%tmin / (wst(iwst)%weat%tmin + 237.3)))                       !kPa
            eair = wst(iwst)%weat%rhum * eas                                                                        !kPa, FAO-56 Eq. 19
            pres = 101.3 * ((293. - 0.0065 * ob(iob)%elev) / 293.) ** 5.26                                           !kPa, FAO-56 Eq. 7
            xlam = 2.501 - 0.002361 * ta_d
            gma = 0.00163 * pres / xlam
            dlt = 4098. * 0.6108 * exp(17.27 * ta_d / (ta_d + 237.3)) / (ta_d + 237.3) ** 2
            if (wst(iwst)%weat%lwrad > -97.) then
              rnet = 0.92 * wst(iwst)%weat%solrad + 0.97 * wst(iwst)%weat%lwrad                                        &
                     - 0.97 * 4.903e-9 * (res_ob(jres)%tw + 273.15) ** 4
            else
              rso = (0.75 + 2.e-5 * ob(iob)%elev) * wst(iwst)%weat%solradmx * 37.59 / 30.                             !FAO-56 Eq. 37
              rnl = 4.903e-9 * 0.5 * ((wst(iwst)%weat%tmax + 273.16) ** 4 + (wst(iwst)%weat%tmin + 273.16) ** 4)        &
                    * (0.34 - 0.14 * sqrt(eair)) * (1.35 * min(wst(iwst)%weat%solrad / max(rso, 1.e-6), 1.) - 0.35)
              rnet = 0.92 * wst(iwst)%weat%solrad - rnl - 4. * 4.903e-9 * (ta_d + 273.15) ** 3 * (res_ob(jres)%tw - ta_d)
            end if
            area_km2 = res_wat_d(jres)%area_ha / 100.
            hmix = min(res(jres)%flo / (res_wat_d(jres)%area_ha * 10000.), 45.7 * sqrt(area_km2) / (21.4 + sqrt(area_km2)))
            dheat = -4.18 * hmix * (res_ob(jres)%tw - tw30) / 30.                                                    !rho c = 4.18 MJ m-3 K-1
            u2 = max(wst(iwst)%weat%windsp, 0.1)
            tk = ta_d + 273.15
            visc = 2.964e-7 * tk ** 1.5 / pres
            ustar = 0.05 * u2
            do it = 1, 50
              ustar = 0.41 * u2 / log(2. * ustar / (0.135 * visc))
            end do
            cu = 3966. / (tk * log(2. / (0.624 * visc / ustar)) * log(2. / (0.135 * visc / ustar)))
            ea_mm = cu * u2 * 10. * (eas - eair)
            evap_mm = max(0., dlt / (dlt + gma) * (rnet + dheat) / xlam + gma / (dlt + gma) * ea_mm)
          end if
          res_wat_d(jres)%evap = 10. * evap_mm * res_wat_d(jres)%area_ha                                             !m3
        else
          res_wat_d(jres)%evap = 10. * res_hyd(jres)%evrsv * wst(iwst)%weat%pet * res_wat_d(jres)%area_ha
        end if
        res_wat_d(jres)%precip = 10. * wst(iwst)%weat%precip * res_wat_d(jres)%area_ha
        if(bsn_cc%gwflow == 0) then !if gwflow active, seepage calculated in gwflow_simulate (rtb gwflow)
          res_wat_d(jres)%seep = 240. * res_hyd(jres)%k * res_wat_d(jres)%area_ha
        endif

        !! add precip to reservoir storage
        res(jres)%flo = res(jres)%flo + res_wat_d(jres)%precip

        !! subtract outflow from reservoir storage
        res(jres)%flo = res(jres)%flo - ht2%flo
        if (res(jres)%flo < 0.) then
          ht2%flo = ht2%flo + res(jres)%flo
          res(jres)%flo = 0.
        end if

        !! subtract evaporation from reservoir storage
        res(jres)%flo = res(jres)%flo - res_wat_d(jres)%evap
        if (res(jres)%flo < 0.) then
          res_wat_d(jres)%evap = res_wat_d(jres)%evap + res(jres)%flo
          res(jres)%flo = 0.
        end if
      
        !! subtract seepage from reservoir storage
        res(jres)%flo = res(jres)%flo - res_wat_d(jres)%seep
        if (res(jres)%flo < 0.) then
          res_wat_d(jres)%seep = res_wat_d(jres)%seep + res(jres)%flo
          res(jres)%flo = 0.
        end if

        !! update surface area - GLOBathy stage from volume, then area from stage (anchored to psa by a_scale); else br1/br2
        if (res(jres)%flo > 0.) then
          if (res_hyd(jres)%hA_a > 0. .and. res_hyd(jres)%hV_a > 0. .and. res_hyd(jres)%hV_b > 0.) then
            stage = (res(jres)%flo * 1.e-9 / res_hyd(jres)%hV_a) ** (1. / res_hyd(jres)%hV_b)
            res_wat_d(jres)%area_ha = res_ob(jres)%a_scale * 100. * res_hyd(jres)%hA_a * stage ** res_hyd(jres)%hA_b
          else
            res_wat_d(jres)%area_ha = res_ob(jres)%br1 * res(jres)%flo ** res_ob(jres)%br2
          end if
        else
          res_wat_d(jres)%area_ha = 0.
        end if

        !! subtract sediment leaving from reservoir
        !res(jres)%sed = max (0., res(jres)%sed - ht2%sed)
        !res(jres)%sil = max (0., res(jres)%sil - ht2%sil)
        !res(jres)%cla = max (0., res(jres)%cla - ht2%cla)
          
        !! perform reservoir nutrient balance
        call res_nutrient (iob)

        !! perform reservoir pesticide transformations
        if (cs_db%num_pests > 0) then
          call res_pest (jres)
          obcs(icmd)%hd(1)%pest = hcs2%pest
        end if

        !! perform reservoir salt process (rtb salt)
        if (cs_db%num_salts > 0) then
          call res_salt(jres)
          obcs(icmd)%hd(1)%salt = hcs2%salt
        endif
        
        !! perform reservoir constituent process (rtb cs)
        if(cs_db%num_cs > 0) then
          icon = res_dat(idat)%cs
          call res_cs(jres, icon, iob)
          obcs(icmd)%hd(1)%cs = hcs2%cs
        endif
        
        !! set values for outflow variables
        ob(icmd)%hd(1) = ht2
        if (cs_db%num_tot > 0) then
          obcs(icmd)%hd(1) = hcs2
        end if

        !! total incoming to output to SWIFT
        ob(icmd)%hin_tot = ob(icmd)%hin_tot + ob(icmd)%hin     
        !! total outgoing to output to SWIFT
        ob(icmd)%hout_tot = ob(icmd)%hout_tot + ht2
        
        if (time%step > 1) then
          do ii = 1, time%step
            ob(icmd)%ts(1,ii) = ht2 / real(time%step)
          end do
        end if

        !! set inflow and outflow variables for reservoir_output
        if (time%yrs > pco%nyskip) then
          res_in_d(jres) = ht1 
          res_out_d(jres) = ht2
          !res_in_d(jres)%flo = res_in_d(jres)%flo / 10000.          !m^3 -> ha-m
          !res_out_d(jres)%flo = res_out_d(jres)%flo / 10000.        !m^3 -> ha-m
          !res_wat_d(jres)%evap = res_wat_d(jres)%evap / 10000.      !m^3 -> ha-m
          !res_wat_d(jres)%seep = res_wat_d(jres)%seep / 10000.      !m^3 -> ha-m
          !res_wat_d(jres)%precip = res_wat_d(jres)%precip / 10000.  !m^3 -> ha-m
        end if             
        
      else
        !! reservoir has not been constructed yet
        ob(icmd)%hd(1) = ob(icmd)%hin
        if (cs_db%num_tot > 0) obcs(icmd)%hd(1) = obcs(icmd)%hin(1)
      end if
      
      !! if called from water allocation, set flag to indicate reservoir has been processed for the day
      res_ob(jres)%wallo_call = 1

      return
      end subroutine res_control

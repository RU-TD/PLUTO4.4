module prizmo
  use prizmo_commons
  use prizmo_core
  use prizmo_loaders
  use prizmo_attenuate
  integer,parameter::prizmo_nspecies=nspecies
  integer,parameter::prizmo_nphoto=nphoto
  integer,parameter::prizmo_nreactions=nreactions

  !! PREPROCESS_INDEXES

  integer,parameter::prizmo_idx_C=idx_C
  integer,parameter::prizmo_idx_CH=idx_CH
  integer,parameter::prizmo_idx_CH2=idx_CH2
  integer,parameter::prizmo_idx_CH2j=idx_CH2j
  integer,parameter::prizmo_idx_CH3=idx_CH3
  integer,parameter::prizmo_idx_CH3j=idx_CH3j
  integer,parameter::prizmo_idx_CH4=idx_CH4
  integer,parameter::prizmo_idx_CH4j=idx_CH4j
  integer,parameter::prizmo_idx_CH5j=idx_CH5j
  integer,parameter::prizmo_idx_CHj=idx_CHj
  integer,parameter::prizmo_idx_CO=idx_CO
  integer,parameter::prizmo_idx_CO_DUST=idx_CO_DUST
  integer,parameter::prizmo_idx_COj=idx_COj
  integer,parameter::prizmo_idx_Cj=idx_Cj
  integer,parameter::prizmo_idx_Cjj=idx_Cjj
  integer,parameter::prizmo_idx_Cjjj=idx_Cjjj
  integer,parameter::prizmo_idx_Cjjjj=idx_Cjjjj
  integer,parameter::prizmo_idx_E=idx_E
  integer,parameter::prizmo_idx_H=idx_H
  integer,parameter::prizmo_idx_H2=idx_H2
  integer,parameter::prizmo_idx_H2O=idx_H2O
  integer,parameter::prizmo_idx_H2O_DUST=idx_H2O_DUST
  integer,parameter::prizmo_idx_H2Oj=idx_H2Oj
  integer,parameter::prizmo_idx_H2j=idx_H2j
  integer,parameter::prizmo_idx_H3Oj=idx_H3Oj
  integer,parameter::prizmo_idx_H3j=idx_H3j
  integer,parameter::prizmo_idx_HCOj=idx_HCOj
  integer,parameter::prizmo_idx_He=idx_He
  integer,parameter::prizmo_idx_Hej=idx_Hej
  integer,parameter::prizmo_idx_Hejj=idx_Hejj
  integer,parameter::prizmo_idx_Hj=idx_Hj
  integer,parameter::prizmo_idx_O=idx_O
  integer,parameter::prizmo_idx_O2=idx_O2
  integer,parameter::prizmo_idx_O2j=idx_O2j
  integer,parameter::prizmo_idx_OH=idx_OH
  integer,parameter::prizmo_idx_OHj=idx_OHj
  integer,parameter::prizmo_idx_Oj=idx_Oj
  integer,parameter::prizmo_idx_Ojj=idx_Ojj
  integer,parameter::prizmo_idx_Ojjj=idx_Ojjj
  integer,parameter::prizmo_idx_Ojjjj=idx_Ojjjj

  !! PREPROCESS_END

contains

  ! ************************
  subroutine prizmo_init()
    implicit none

    runtime_data_folder = "runtime_data/"

    call load_energy()
    call load_all_photo_xsecs()
    call load_energy_thresholds()
    call load_all_atomic_cooling_tables()
    call load_dust_cooling_table()
    !call load_photoelectric_tables() ! FIXME
    call load_dust_kappa_opacity()
    call load_H2_cooling_tabs()
    call load_shielding_H2_table()
    call load_shielding_CO_table()
    call load_CO_cooling()
    call load_verbatim_reactions()

    print *, "everything loaded!"

    gamma_ad = 7./5.
    d2g = 1d-2
    ortho_to_para = 3. / 1.
    user_Av = 0d0
    radial_Ncol_H2 = 0d0  ! 1/cm2
    radial_Ncol_CO = 0d0  ! 1/cm2
    vertical_Ncol_CO = 0d0  ! 1/cm2
    user_cr = 5d-17  ! 1/s
    solve_thermo = .true.
    solve_chemistry = .true.

    ode_atol = 1d-20  ! default absolute tolerance
    ode_rtol = 1d-8  ! default relative tolerance

  end subroutine prizmo_init

  ! ************************
  subroutine prizmo_evolve(x, Tgas, jflux, dt)
    implicit none
    real*8,intent(inout)::x(nspecies), Tgas, jflux(nphoto)
    real*8,intent(in)::dt

    call init(x, Tgas, jflux)
    call evolve(x, Tgas, jflux, dt)

  end subroutine prizmo_evolve

  ! ************************
  subroutine prizmo_evolve_rho(xx, rho, Tgas, jflux, dt)
    use prizmo_commons
    implicit none
    real*8,intent(inout)::xx(nspecies), Tgas, jflux(nphoto)
    real*8,intent(in)::dt, rho
    real*8::x(nspecies)

    x = rho * xx / masses

    call prizmo_evolve(x, Tgas, jflux, dt)

    xx = x * masses / rho

  end subroutine prizmo_evolve_rho

  ! **********************
  subroutine prizmo_rt(x, Tgas, jflux, ds)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, ds
    real*8,intent(inout)::jflux(nphoto)

    call attenuate(x, Tgas, jflux, ds)

  end subroutine prizmo_rt

  ! **********************
  subroutine prizmo_rt_rho(xx, rho, Tgas, jflux, ds)
    implicit none
    real*8,intent(in)::xx(nspecies), rho, Tgas, ds
    real*8,intent(inout)::jflux(nphoto)
    real*8::x(nspecies)

    x = rho * xx / masses

    call attenuate(x, Tgas, jflux, ds)

  end subroutine prizmo_rt_rho

  ! **************
  function prizmo_get_rho(x) result(rho)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::rho

    rho = sum(masses * x)

  end function prizmo_get_rho

  ! **************
  subroutine prizmo_set_atol_all(val)
    implicit none
    real*8,intent(in)::val

    ode_atol = val

  end subroutine prizmo_set_atol_all

  ! **************
  subroutine prizmo_set_atol_idx(val, idx)
    implicit none
    real*8,intent(in)::val
    integer,intent(in)::idx

    ode_atol(idx) = val

  end subroutine prizmo_set_atol_idx

  ! **************
  function prizmo_load_radiation_field(filename) result(field)
    implicit none
    character(len=*),intent(in)::filename
    real*8::field(nphoto)
    integer::i, unit

    open(newunit=unit, file=trim(filename), status="old")
    do i=1,nphoto
      read(unit, *) field(i)
    end do
    close(unit)

  end function prizmo_load_radiation_field

  ! ! ********************
  ! function prizmo_get_atomic_cooling_array(x, Tgas) result(coola)
  !   use prizmo_cooling_atomic
  !   implicit none
  !   real*8,intent(in)::x(nspecies), Tgas
  !   real*8::coola(atomic_cooling_nvec)
  !
  !   coola = cooling_atomic_array(x, log10(Tgas))
  !
  ! end function prizmo_get_atomic_cooling_array

  ! ******************
  function prizmo_get_cooling_array(x, Tgas, Tdust, jflux) result(cools)
    use prizmo_cooling
    use prizmo_flux
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::fluxes(nreactions), cools(ncooling)

    call init(x, Tgas, jflux)

    fluxes(:) = get_flux(x, Tgas, Tdust)
    cools(:) = cooling_array(x, Tgas, Tdust, jflux, fluxes)

  end function prizmo_get_cooling_array

  ! ******************
  function prizmo_get_heating_array(x, Tgas, Tdust, jflux) result(heats)
    use prizmo_heating
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::heats(nheating)

    call init(x, Tgas, jflux)
    heats(:) = heating_array(x, Tgas, Tdust, jflux)

  end function prizmo_get_heating_array

  ! ***************************
  subroutine prizmo_save_cooling_function(x, jflux, fname)
    use prizmo_cooling
    use prizmo_flux
    implicit none
    integer,parameter::ntemp=100
    real*8,intent(in)::x(nspecies), jflux(nphoto)
    character(len=*),intent(in)::fname
    real*8::fluxes(nreactions), cools(ncooling), tmin, tmax, Tgas, Tdust
    integer::unit, i

    tmin = log10(1d0)
    tmax = log10(1d5)

    open(newunit=unit, file=trim(fname), status="replace")
      do i=1, ntemp
        tgas = 1d1**((i - 1) * (tmax - tmin) / (ntemp - 1) + tmin)
        tdust = tgas
        fluxes(:) = get_flux(x, Tgas, Tdust)
        cools(:) = cooling_array(x, Tgas, Tdust, jflux, fluxes)
        write(unit, '(99e17.8e3)') tgas, cools
      end do
    close(unit)

  end subroutine prizmo_save_cooling_function

  ! ***************************
  subroutine prizmo_save_heating_function(x, jflux, fname)
    use prizmo_heating
    implicit none
    integer,parameter::ntemp=100
    real*8,intent(in)::x(nspecies), jflux(nphoto)
    character(len=*),intent(in)::fname
    real*8::heats(nheating), tmin, tmax, Tgas, Tdust
    integer::unit, i

    tmin = log10(1d0)
    tmax = log10(1d5)

    open(newunit=unit, file=trim(fname), status="replace")
      do i=1, ntemp
        tgas = 1d1**((i - 1) * (tmax - tmin) / (ntemp - 1) + tmin)
        tdust = tgas
        heats(:) = heating_array(x, Tgas, Tdust, jflux)
        write(unit, '(99e17.8e3)') tgas, heats
      end do
    close(unit)

  end subroutine prizmo_save_heating_function

  ! ****************************
  function prizmo_get_energy()
    implicit none
    real*8::prizmo_get_energy(nphoto)

    prizmo_get_energy = energy

  end function prizmo_get_energy

  ! ****************************
  function prizmo_get_energy_ev()
    implicit none
    real*8::prizmo_get_energy_ev(nphoto)

    prizmo_get_energy_ev = energy / ev2erg

  end function prizmo_get_energy_ev

  ! ****************************
  subroutine prizmo_set_user_Av(Av)
    implicit none
    real*8,intent(in)::Av

    user_Av = Av

  end subroutine prizmo_set_user_Av

  ! ****************************
  subroutine prizmo_set_d2g(val)
    implicit none
    real*8,intent(in)::val

    d2g = val

  end subroutine prizmo_set_d2g

  ! ****************************
  subroutine prizmo_set_crate(val)
    implicit none
    real*8,intent(in)::val

    user_cr = val

  end subroutine prizmo_set_crate

  ! ****************************
  subroutine prizmo_set_radial_Ncol_H2(val)
    implicit none
    real*8,intent(in)::val

    radial_Ncol_H2 = val

  end subroutine prizmo_set_radial_Ncol_H2

  ! ****************************
  subroutine prizmo_set_radial_Ncol_CO(val)
    implicit none
    real*8,intent(in)::val

    radial_Ncol_CO = val

  end subroutine prizmo_set_radial_Ncol_CO

  ! ****************************
  subroutine prizmo_set_vertical_Ncol_CO(val)
    implicit none
    real*8,intent(in)::val

    vertical_Ncol_CO = val

  end subroutine prizmo_set_vertical_Ncol_CO

  ! ********************
  function prizmo_get_rates(x, Tgas, Tdust, jflux) result(k)
    use prizmo_rates
    use prizmo_rates_photo
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::k(nreactions)

    call compute_rates(x, Tgas, Tdust)
    call compute_photorates(x, Tgas, jflux)

    k = kall

  end function prizmo_get_rates

  ! ***************
  function prizmo_get_tdust(x, Tgas, jflux) result(tdust)
    implicit none
    real*8,intent(in)::jflux(nphoto), Tgas, x(nspecies)
    real*8::tdust, log_Ea, Eabsorption, f(nphoto), log_ngas, log_tgas

    log_Tgas = log10(Tgas)
    log_ngas = log10(sum(x))

    call compute_Eabsorption(jflux)

    log_Tgas = min(log_Tgas, tdust_table_data%ymax*0.999999)  ! FIXME

    tdust = 1d1**interp_3dfit(log_Eabsorption, log_Tgas, log_ngas, tdust_table_data, &
      dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3)

  end function prizmo_get_tdust

  ! ******************
  subroutine prizmo_set_solve_thermo(val)
    implicit none
    logical,intent(in)::val

    solve_thermo = val

  end subroutine prizmo_set_solve_thermo

  ! ******************
  subroutine prizmo_set_solve_chemistry(val)
    implicit none
    logical,intent(in)::val

    solve_chemistry = val

  end subroutine prizmo_set_solve_chemistry

  ! ***********************
  subroutine print_ranked_fluxes(x, Tgas, Tdust, jflux)
    use prizmo_rates
    use prizmo_rates_photo
    use prizmo_flux
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto)
    real*8::fluxes(nreactions)

    call compute_rates(x, Tgas, Tdust)
    call compute_photorates(x, Tgas, jflux)

    fluxes = get_flux(x, Tgas, Tdust)

    call ranker(fluxes, 60, reactions_verbatim)

  end subroutine print_ranked_fluxes

  ! ***************
  function prizmo_get_Xnuclei(x, atom) result(nX)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies)
    character(len=*),intent(in)::atom
    real*8::nX

    if(trim(atom) == "H") then
      nX = get_Hnuclei(x)
    else if(trim(atom) == "C") then
      nX = get_Cnuclei(x)
    else if(trim(atom) == "O") then
      nX = get_Onuclei(x)
    else if(trim(atom) == "He") then
      nX = get_HEnuclei(x)
    else
      print *, "ERROR: unknown atom "//trim(atom)//" in prizmo_get_Xnuclei"
      stop
    end if

  end function prizmo_get_Xnuclei

  ! ************************
  function prizmo_get_chi_FUV() result(xfuv)
    implicit none
    real*8::xfuv

    xfuv = chi_FUV

  end function prizmo_get_chi_FUV

end module prizmo

module prizmo_c
  use prizmo
  use iso_c_binding
contains

  ! *********************
  ! initialize prizmo, to be called once
  subroutine prizmo_init_c() bind(C)

    call prizmo_init()

  end subroutine prizmo_init_c

  ! ************************
  ! evolve chemistry and temperature
  ! x(prizmo_nspecies): species abundances, 1/cm3
  ! Tgas: gas temperature, K
  ! jflux(prizmo_nphoto): radiation per bin, erg/cm2/s/Hz
  ! dt: time interval, s
  ! returns: updated x(prizmo_nspecies) and Tgas
  subroutine prizmo_evolve_c(x, Tgas, jflux, dt) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::x(nspecies), Tgas, jflux(nphoto)
    real(C_DOUBLE),intent(in)::dt

    call prizmo_evolve(x, Tgas, jflux, dt)

  end subroutine prizmo_evolve_c

  ! ************************
  ! same as prizmo_evolve_c, but with mass fractions
  ! xx(prizmo_nspecies): species mass fractions, no dimensions
  ! rho: gas mass density, g/cm3
  ! Tgas: gas temperature, K
  ! jflux(prizmo_nphoto): radiation per bin, erg/cm2/s/Hz
  ! dt: time interval, s
  ! returns: updated xx(prizmo_nspecies) and Tgas
  subroutine prizmo_evolve_rho_c(xx, rho, Tgas, jflux, dt) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::xx(nspecies), Tgas, jflux(nphoto)
    real(C_DOUBLE),intent(in)::dt, rho

    call prizmo_evolve_rho(xx, rho, Tgas, jflux, dt)

  end subroutine prizmo_evolve_rho_c

  ! ************************
  ! alias of prizmo_evolve_rho_c
  subroutine prizmo_rho_c(xx, rho, Tgas, jflux, dt) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::xx(nspecies), Tgas, jflux(nphoto)
    real(C_DOUBLE),intent(in)::dt, rho

    call prizmo_evolve_rho_c(xx, rho, Tgas, jflux, dt)

  end subroutine prizmo_rho_c

  ! **********************
  ! attenuate radiation accordingly to chemistry and dust content
  ! x(prizmo_nspecies): species abundances, 1/cm3
  ! Tgas: gas temperature, K
  ! jflux(prizmo_nphoto): radiation per bin, erg/cm2/s/Hz
  ! ds: cell size, cm
  ! returns: updated jflux(prizmo_nphoto)
  subroutine prizmo_rt_c(x, Tgas, jflux, ds) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::x(nspecies), Tgas, ds
    real(C_DOUBLE),intent(inout)::jflux(nphoto)

    call prizmo_rt(x, Tgas, jflux, ds)

  end subroutine prizmo_rt_c

  ! **********************
  ! equivalent to prizmo_rt_c, but using mass fractions
  ! xx(prizmo_nspecies): species mass fractions, no dimensions
  ! rho: gas mass density, g/cm3
  ! Tgas: gas temperature, K
  ! jflux(prizmo_nphoto): radiation per bin, erg/cm2/s/Hz
  ! ds: cell size, cm
  ! returns: updated jflux(prizmo_nphoto)
  subroutine prizmo_rt_rho_c(xx, rho, Tgas, jflux, ds) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::xx(nspecies), rho, Tgas, ds
    real(C_DOUBLE),intent(inout)::jflux(nphoto)

    call prizmo_rt_rho(xx, rho, Tgas, jflux, ds)

  end subroutine prizmo_rt_rho_c

  ! **********************
  ! compute mass fraction from abundances
  ! x(prizmo_nspecies): species abundances, 1/cm3
  ! returns: rho, gas mass density, g/cm3
  subroutine prizmo_get_rho_c(x, rho) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(in)::x(nspecies)
    real(C_DOUBLE),intent(inout)::rho

    rho = prizmo_get_rho(x)

  end subroutine prizmo_get_rho_c

  ! **********************
  ! convert number densities to fractions
  ! x(prizmo_nspecies): species abundances, 1/cm3
  ! returns: xx(prizmo_nspecies), species mass fractions, no dimensions
  subroutine prizmo_n2frac_c(x, xx) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(in)::x(nspecies)
    real(C_DOUBLE),intent(inout)::xx(nspecies)

    xx = x * masses / sum(x * masses)

  end subroutine prizmo_n2frac_c

  ! **********************
  ! convert fractions to number densities
  ! xx(prizmo_nspecies): mass fractions, no dimension
  ! rho: gas mass density, g/cm3
  ! returns: x(prizmo_nspecies), species abundances, 1/cm3
  subroutine prizmo_frac2n_c(xx, rho, x) bind(C)
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(in)::xx(nspecies), rho
    real(C_DOUBLE),intent(inout)::x(nspecies)

    x = xx * rho / masses

  end subroutine prizmo_frac2n_c

  ! ******************
  ! set if cooling and heating will be solved
  ! 0 = NO cooling nor heating solved
  ! 1 = cooling and heating are solved
  subroutine prizmo_set_solve_thermo_c(val) bind(C)
    implicit none
    integer(C_INT),intent(in)::val

    if(val == 0) then
      solve_thermo = .false.
    elseif(val == 1) then
      solve_thermo = .true.
    else
      print *, "ERROR: solve thermo can be either 0 or 1!"
      stop
    end if

  end subroutine prizmo_set_solve_thermo_c

  ! ******************
  ! set if chemistry will be solved
  ! 0 = NO chemistry solved
  ! 1 = chemistry solved
  subroutine prizmo_set_solve_chemistry_c(val) bind(C)
    implicit none
    integer(C_INT),intent(in)::val

    if(val == 0) then
      solve_chemistry = .false.
    elseif(val == 1) then
      solve_chemistry = .true.
    else
      print *, "ERROR: solve chemistry can be either 0 or 1!"
      stop
    end if

  end subroutine prizmo_set_solve_chemistry_c

  ! ****************************
  ! set the radial H2 column density to compute shielding
  ! val in 1/cm2
  subroutine prizmo_set_radial_ncol_h2_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_radial_Ncol_H2(val)

  end subroutine prizmo_set_radial_ncol_h2_c

  ! ****************************
  ! set the radial CO column density to compute shielding
  ! val in 1/cm2
  subroutine prizmo_set_radial_ncol_co_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_radial_Ncol_CO(val)

  end subroutine prizmo_set_radial_ncol_co_c

  ! ****************************
  ! set the radial CO column density to compute CO cooling escape
  ! val in 1/cm2
  subroutine prizmo_set_vertical_ncol_co_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_vertical_Ncol_CO(val)

  end subroutine prizmo_set_vertical_ncol_co_c

  ! ****************************
  ! set dust/gas mass ratio, default: 1e-2
  subroutine prizmo_set_d2g_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_d2g(val)

  end subroutine prizmo_set_d2g_c

  ! ****************************
  ! set cosmic rays ionization rate, default 5e-17 1/s
  subroutine prizmo_set_crate_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_crate(val)

  end subroutine prizmo_set_crate_c

  ! ****************************
  ! set absolute tolerances for all the species
  subroutine prizmo_set_atol_all_c(val) bind(C)
    implicit none
    real(C_DOUBLE),intent(in)::val

    call prizmo_set_atol_all(val)

  end subroutine prizmo_set_atol_all_c

end module prizmo_c

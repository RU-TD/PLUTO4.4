module prizmo_main_c
  use prizmo_main
  use prizmo_user_photo

contains

  ! *******************
  ! initialize LKRM, C interface
  subroutine prizmo_init_c() bind(C)
    use iso_c_binding

    call prizmo_init()

  end subroutine prizmo_init_c

  ! *******************
  ! do chemistry, C interface
  subroutine prizmo_c(n, Tgas, jflux, dt) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::n(nmols), Tgas
    real(C_DOUBLE),intent(in)::dt, jflux(nphoto)

    ! call chemistry
    call prizmo(n(:), Tgas, jflux(:), dt)

  end subroutine prizmo_c

  ! *********************
  ! call chemistry using mass fractions (x),
  ! mass density (rho, g/cm3), temperature (Tgas, K)
  ! and time-step (dt, s)
  subroutine prizmo_rho_c(x, rho, tgas, jflux, dt) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::x(nmols), tgas
    real(C_DOUBLE),intent(in)::dt, rho, jflux(nphoto)
    real*8::n(nmols)

    ! convert mass fractions to number density
    n(:) = x(:) * rho / mass(:)

    ! call chemistry
    call prizmo(n(:), Tgas, jflux(:), dt)

    ! convert number density to mass fractions
    x(:) = n(:) * mass(:) / rho

  end subroutine prizmo_rho_c

  ! *********************
  subroutine prizmo_set_cooling_mode_c(mode) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    integer(C_INT),intent(in)::mode

    call prizmo_set_cooling_mode(mode)

  end subroutine prizmo_set_cooling_mode_c

  ! *********************
  subroutine prizmo_set_heating_mode_c(mode) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    integer(C_INT),intent(in)::mode

    call prizmo_set_heating_mode(mode)

  end subroutine prizmo_set_heating_mode_c

  ! *********************
  subroutine prizmo_set_tdust_mode_c(mode) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    integer(C_INT),intent(in)::mode

    call prizmo_set_tdust_mode(mode)

  end subroutine prizmo_set_tdust_mode_c

  ! ***********************
  subroutine prizmo_get_species_mass_c(masses) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(out)::masses(nmols)

    masses(:) = prizmo_get_species_mass()

  end subroutine prizmo_get_species_mass_c

  ! ***********************
  subroutine prizmo_get_cooling_c(cooling, n, Tgas, jflux)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(in)::n(nmols), Tgas, jflux(nphoto)
    real(C_DOUBLE),intent(out)::cooling

    cooling = prizmo_get_cooling(n(:), Tgas, jflux(:))

  end subroutine prizmo_get_cooling_c

  ! ********************************
  subroutine prizmo_get_draine_flux_c(jflux) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(out)::jflux(nphoto)

    jflux(:) = prizmo_get_draine_flux()

  end subroutine prizmo_get_draine_flux_c

  ! ************************
  subroutine prizmo_attenuate_rho_c(jflux, x, Tgas, length, rho) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::jflux(nphoto)
    real(C_DOUBLE),intent(in)::x(nmols), Tgas, length, rho
    real*8::n(nmols)

    n(:) = x(:) * rho / mass(:)
    call prizmo_attenuate(jflux(:), n(:), Tgas, length)

  end subroutine prizmo_attenuate_rho_c

  ! ************************
  subroutine prizmo_x2n_c(x, n, rho) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::x(nmols)
    real(C_DOUBLE),intent(in)::rho
    real(C_DOUBLE),intent(out)::n(nmols)

    n(:) = x(:) * rho / mass(:)

  end subroutine prizmo_x2n_c

  ! ************************
  subroutine prizmo_n2x_c(n, x, rho) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::n(nmols)
    real(C_DOUBLE),intent(out)::x(nmols), rho

    rho = sum(n(:) * mass(:))
    x(:) = n(:) * mass(:) / rho

  end subroutine prizmo_n2x_c

  ! ************************
  subroutine prizmo_x2ngas_c(x, rho, ngas) bind(C)
    use iso_c_binding
    use prizmo_commons
    implicit none
    real(C_DOUBLE),intent(inout)::x(nmols)
    real(C_DOUBLE),intent(in)::rho
    real(C_DOUBLE),intent(out)::ngas

    ngas = sum(x(:) * rho / mass(:))

  end subroutine prizmo_x2ngas_c


  !!BEGIN_USER_VARIABLES_FUNCTIONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! *************************
    ! user variable VARIABLE_AV initialization function (with C binding)
    subroutine prizmo_set_variable_av_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_Av(arg)

    end subroutine prizmo_set_variable_av_c

    ! *************************
    ! user variable VARIABLE_G0 initialization function (with C binding)
    subroutine prizmo_set_variable_g0_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_G0(arg)

    end subroutine prizmo_set_variable_g0_c

    ! *************************
    ! user variable VARIABLE_CRFLUX initialization function (with C binding)
    subroutine prizmo_set_variable_crflux_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_crflux(arg)

    end subroutine prizmo_set_variable_crflux_c

    ! *************************
    ! user variable VARIABLE_NCO_INCOMING initialization function (with C binding)
    subroutine prizmo_set_variable_nco_incoming_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_NCO_incoming(arg)

    end subroutine prizmo_set_variable_nco_incoming_c

    ! *************************
    ! user variable VARIABLE_NCO_ESCAPING initialization function (with C binding)
    subroutine prizmo_set_variable_nco_escaping_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_NCO_escaping(arg)

    end subroutine prizmo_set_variable_nco_escaping_c

    ! *************************
    ! user variable VARIABLE_NH2_INCOMING initialization function (with C binding)
    subroutine prizmo_set_variable_nh2_incoming_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_NH2_incoming(arg)

    end subroutine prizmo_set_variable_nh2_incoming_c

    ! *************************
    ! user variable VARIABLE_NH2_ESCAPING initialization function (with C binding)
    subroutine prizmo_set_variable_nh2_escaping_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_NH2_escaping(arg)

    end subroutine prizmo_set_variable_nh2_escaping_c

    ! *************************
    ! user variable VARIABLE_NH2O_ESCAPING initialization function (with C binding)
    subroutine prizmo_set_variable_nh2o_escaping_c(arg) bind(C)
        use iso_c_binding
        implicit none
        real(C_DOUBLE),intent(in)::arg

        call prizmo_set_variable_NH2O_escaping(arg)

    end subroutine prizmo_set_variable_nh2o_escaping_c

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_USER_VARIABLES_FUNCTIONS

end module prizmo_main_c
module prizmo_H2_dust
integer,parameter::ngrid=1000
real*8::H2_table_tgas(ngrid), xmin, xfact, invdx
real*8::H2_table_tdust(ngrid), ymin, yfact, invdy
real*8::H2_table_eff(ngrid, ngrid)

contains

  ! ******************************
  ! compute H+H->H2 rate
  function get_H2_dust(Tgas_in, Tdust_in, pre_scale, nH) result(k)
    use prizmo_commons
    use prizmo_fit
    use prizmo_dust_evaporation
    implicit none
    real*8,intent(in)::Tgas_in, Tdust_in, pre_scale, nH
    real*8::k, eff, Tgas, Tdust

    Tgas = min(Tgas_in, H2_table_tgas(ngrid))
    Tdust = min(Tdust_in, H2_table_tdust(ngrid))

    ! this is vgas * epsilon * stick
    eff = 1d1**interp_2d(log10(Tgas), log10(Tdust), H2_table_eff(:, :), &
      xmin, xfact, invdx, &
      ymin, yfact, invdy)

    ! rate cm3/s (divided by the number density of H, because RHS is k*nH*nH)
    k = 0.5 * pre_scale * eff / (nH + 1d-40) * get_dust_fevap(Tdust)

  end function get_H2_dust

  ! ******************************
  function get_H2_dust_n(Tgas, Tdust, n) result(k)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::Tgas, Tdust, n(nmols)
    real*8::ntot, mu, rho_dust, k, pre_scale

    ! get total number density
    ntot = get_ntot(n(:))

    ! mean molecular weight
    mu = get_mu(n(:))

    ! dust mass density, g/cm3
    rho_dust = ntot * mu * dust_d2g * proton_mass

    ! dust content scaling factor
    pre_scale = rho_dust * dust_fact3 / dust_fact4 / pi43 / dust_rho_bulk

    k = get_H2_dust(Tgas, Tdust, pre_scale, n(idx_H))

  end function get_H2_dust_n

  ! ******************************
  subroutine load_H2_dust_table()
    use prizmo_commons
    use prizmo_fit
    implicit none

    print *, "loading H2 dust table..."

    call load_data_2d(runtime_folder//"H2_dust_factor.dat", H2_table_eff(:, :), &
       H2_table_tgas(:), ngrid, xmin, xfact, invdx, &
       H2_table_tdust(:), ngrid, ymin, yfact, invdy)

  end subroutine load_H2_dust_table

end module prizmo_H2_dust

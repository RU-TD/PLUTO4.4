module prizmo_heating_photoelectric
  use prizmo_commons
  real*8::phe_fact(nphoto, dust_zmin:dust_zmax), phe_fact_sum(dust_zmin:dust_zmax)
  real*8::kpe_fact(nphoto, dust_zmin:dust_zmax), kpe_fact_sum(dust_zmin:dust_zmax)
  real*8::kpd_fact(nphoto, dust_zmin:dust_zmax), kpd_fact_sum(dust_zmin:dust_zmax)
  real*8::latest_fZ(dust_zmin:dust_zmax)
  integer,parameter::imin=1

  !!BEGIN_PHOTOHEATING_FIT_VARIABLES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    integer,parameter::fit_jei_nsteps=100
    integer,parameter::particle_zmax=2

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_PHOTOHEATING_FIT_VARIABLES

  real*8::fit_jei_xmin, fit_jei_xfact, fit_jei_invdx
  real*8::fit_jei_fdata(fit_jei_nsteps, dust_zmin:dust_zmax, -1:particle_zmax)

contains

  ! *********************
  ! compute photoelectric heating using G0 and Av, erg/s/cm3
  function heating_photoelectric(n, Tgas, G0, Av) result(heat)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, G0, Av
    real*8::heat, eps, psi, nH

    if(dust_d2g < 1d-6) then
      heat = 0d0
      return
    end if

    nH = get_Hnuclei(n(:))
    psi = G0 * sqrt(Tgas) / (n(idx_e) + 1d-40)
    eps = 4.9d-2 / (1d0 + (psi / 1.925d3)**0.73) &
         + 3.7d-2 * (Tgas / 1d4)**0.7 / (1d0 + psi / 5d3)
    heat = 1.3d-24 * eps * G0 * nH * exp(-2.5 * Av)

  end function heating_photoelectric

  ! *********************
  ! compute multifreqency photoelectric heating, erg/s/cm3
  function heating_photoelectric_multi(n, Tgas, Tdust, jflux) result(heat)
    use prizmo_commons
    use prizmo_utils
    use prizmo_dust_evaporation
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::heat, f(dust_ncharge), dust_cnorm, ndust, fevap
    integer::i, iz

    fevap = get_dust_fevap(Tdust)

    ! skip calulation for small amounts of dust
    if(dust_d2g * fevap < 1d-6) then
      heat = 0d0
      return
    end if

    ! precompute heating rate and Jpe rate
    call precompute_fact_sum(n(:), jflux(:))

    ! compute fraction of dust in each charge state
    ! note that the array is 1-based (instead of zmin-based)
    f(:) = compute_fractions(n(:), Tgas)

    heat = 0d0

    do i=dust_zmin, dust_zmax
      iz = i + 1 - dust_zmin
      !print *, i, phe_fact_sum(i), f(iz)
      heat = heat + phe_fact_sum(i) * f(iz)
    end do

    ! normalize to the actual amount of dust
    dust_cnorm = get_rho(n(:)) * dust_d2g / pi43 / dust_rho_bulk / dust_fact4
    ndust = dust_cnorm * dust_fact1  ! 1/cm3
    heat = max(heat * ndust * fevap, 0d0) ! FIXME

  end function heating_photoelectric_multi

  ! *********************
  ! note that arrays here are 1-based (instead of zmin-based)
  function compute_fractions(n, Tgas_in) result(f)
    use prizmo_commons
    use prizmo_linear_solver
    use prizmo_fit
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas_in
    real*8::a(dust_ncharge, dust_ncharge), f(dust_ncharge), Je_rate, Ji_rate
    real*8::logTgas, Tgas
    integer::i, iz, j

    ! check Tgas range
    Tgas = max(min(Tgas_in, 3.9d3), Tgas_min)

    logTgas = log10(Tgas)

    !print *, "**********"
    a(:, :) = 0d0
    a(dust_ncharge, :) = 1d0
    do i=1, dust_ncharge - 1
      iz = i - 1 + dust_zmin
      ! electron-dust recombination rate
      Je_rate = n(idx_e) / sqrt(mass(idx_e)) &
        * 1d1**interp_1d(logTgas, fit_jei_fdata(:, iz+1, -1), fit_jei_xmin, fit_jei_xfact, fit_jei_invdx)

      ! cation-dust recombination rate
      Ji_rate = 0d0
      Ji_rate = Ji_rate + n(idx_Hj) / sqrt(mass(idx_Hj)) * 1d1**interp_1d(logTgas, fit_jei_fdata(:, iz, 1), fit_jei_xmin, fit_jei_xfact, fit_jei_invdx)
      Ji_rate = Ji_rate + n(idx_Hej) / sqrt(mass(idx_Hej)) * 1d1**interp_1d(logTgas, fit_jei_fdata(:, iz, 1), fit_jei_xmin, fit_jei_xfact, fit_jei_invdx)

      !print '(i4,99e17.8e3)', iz, Tgas, Je_rate, Ji_rate, kpe_fact_sum(iz)

      ! prepare matrix for the bidiagonal solver
      a(i, i) = kpd_fact_sum(iz) + kpe_fact_sum(iz) + Ji_rate
      a(i, i + 1) = Je_rate
    end do

    !a(:, :) = a(:, :) + 1d-40

    !do j=1,size(a, 2)
    !  do i=1,size(a, 1)
    !    print *, i, j, a(i, j)
    !  end do
    !end do

    f(:) = bidiagonal_solver(a(:, :), 1d0, dust_ncharge)

  end function compute_fractions

  ! *********************
  ! precompute the integrals of Jpe and heating
  subroutine precompute_fact_sum(n, jflux)
    use prizmo_commons
    implicit none
    real*8,intent(in)::n(nmols), jflux(nphoto)
    real*8::rho
    integer::i

    !FIXME: starting index of sum on energy bins based on Emin

    ! loop on charges
    do i=dust_zmin, dust_zmax
      phe_fact_sum(i) = sum(jflux(imin:nphoto) * phe_fact(imin:nphoto, i))  ! heating
      kpe_fact_sum(i) = sum(jflux(imin:nphoto) * kpe_fact(imin:nphoto, i))  ! rate pe
      kpd_fact_sum(i) = sum(jflux(imin:nphoto) * kpd_fact(imin:nphoto, i))  ! rate pd
      !print *, "pre", i, kpe_fact_sum(i), phe_fact_sum(i)
    end do

  end subroutine precompute_fact_sum

  ! ********************
  ! load factors for integrals
  subroutine load_photoelectric_data()
    use prizmo_commons
    use prizmo_fit
    implicit none
    integer::unit, unit_kpe, unit_kpd, zcharge, i, j, k, zparticle
    real*8::energy, xdata(fit_jei_nsteps)

    print *, "load photoelectric heating data..."

    ! LOAD Jpe rate and Heating data
    open(newunit=unit, file=runtime_folder//"photoelectric_heating_factors.dat", status="old")
    open(newunit=unit_kpe, file=runtime_folder//"photoelectric_rate_pe_factors.dat", status="old")
    open(newunit=unit_kpd, file=runtime_folder//"photoelectric_rate_pd_factors.dat", status="old")
    do i=dust_zmin, dust_zmax
      do j=1,nphoto
        read(unit, *) zcharge, energy, phe_fact(j, i)
        read(unit_kpe, *) zcharge, energy, kpe_fact(j, i)
        read(unit_kpd, *) zcharge, energy, kpd_fact(j, i)
      end do
    end do
    close(unit)
    close(unit_kpe)

    ! LOAD Je rate data
    open(newunit=unit, file=runtime_folder//"photoelectric_jei_rate.dat", status="old")
    do k=-1, particle_zmax
      do i=dust_zmin, dust_zmax
        do j=1, fit_jei_nsteps
          read(unit, *) zparticle, zcharge, xdata(j), fit_jei_fdata(j, i, k)
        end do
      end do
    end do
    close(unit)

    ! fit is in log space
    xdata(:) = log10(xdata(:))

    ! prepare fit parameters (same for all Z charge)
    call fit_prepare(xdata(:), fit_jei_nsteps, fit_jei_xmin, fit_jei_xfact, fit_jei_invdx)

    ! fit is in log space
    fit_jei_fdata(:, :, :) = log10(fit_jei_fdata(:, :, :) + 1d-80)

  end subroutine load_photoelectric_data

end module prizmo_heating_photoelectric

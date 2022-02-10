! test code to call prizmo
program test
  use prizmo_main
  use prizmo_user_photo
  use prizmo_fit
  implicit none
  integer,parameter::imax=200
  real*8,parameter::spy=365.*24.*3600.
  real*8,parameter::pc2cm=3.08567758128e18, pi=acos(-1d0)
  real*8,parameter::N2Av=6.289d-22
  real*8::x(prizmo_nmols), tend, t, dt, xall(prizmo_nmols, imax)
  real*8::Tgas(imax), jflux(prizmo_nphoto), tau(prizmo_nphoto)
  real*8::crflux, dz, NCO, NH2, NH2O, Ntot, ngas, z(imax)
  real*8::Gnot, k(prizmo_nrea), xold(prizmo_nmols), energy(prizmo_nphoto)
  real*8::zmin_log, zmax_log, zold, Av(imax)
  real*8::Tdust, rout(12+prizmo_nmols), fluxes(prizmo_nrea)
  integer::istep, unit, i, unit_final, j, unit_cooling, unit_cooling_at, unit_emission, ios

  ! init prizmo
  call prizmo_init()

  ! 1/s
  crflux = 5d-17

  ! 1/cm3
  ngas = 1d3

  Gnot = 1d1 / 1.15 ! / 2d0 !/ 1.16

  jflux(:) = 0d0
  energy(:) = prizmo_get_energy()

  x(:) = 0d0
  x(prizmo_idx_H2) = 0.5 * ngas
  x(prizmo_idx_O) = 3d-4 * ngas
  x(prizmo_idx_C) = 1d-4 * ngas
  x(prizmo_idx_Hej) = 1d-1 * ngas
  x(prizmo_idx_E) = prizmo_get_electrons(x(:))

  do i=1,imax
    xall(:, i) = x(:)
  end do

  do i=1,imax
    zmin_log = log10(1e-6 / N2Av / ngas)
    zmax_log = log10(3e1 / N2Av / ngas)
    z(i) = 1e1**((i - 1) * (zmax_log - zmin_log) / (imax - 1) + zmin_log)
    Av(i) = z(i) * ngas * N2Av
  end do

  ! K
  Tgas(:) = 5d1

  open(newunit=unit, file="slab_final.out", status="replace")
  open(newunit=unit_cooling, file="cool.out", status="replace")

  t = 0d0
  dt = spy
  tend = spy * 1d7
  do
    jflux(:) = prizmo_get_draine_flux() * Gnot * 2d0 * pi

    xold(:) = xall(:, 1)
    zold = 0d0
    NH2 = 0d0
    NCO = 0d0
    NH2O = 0d0
    Ntot = 0d0
    dt = dt * 1.1
    print *, t / spy

    do i=1,imax

      !print *, i, imax, Av(i), prizmo_get_Gnot(jflux(:))

      dz = z(i) - zold

      NH2 = NH2 + xold(prizmo_idx_H2) * dz
      NCO = NCO + xold(prizmo_idx_CO) * dz
      NH2O = NH2O + xold(prizmo_idx_H2O) * dz
      Ntot = Ntot + (2d0 * xold(prizmo_idx_H2) + xold(prizmo_idx_H)) * dz

      call prizmo_attenuate(jflux(:), xold(:), Tgas(i), dz)

      call prizmo_set_variable_crflux(crflux)
      call prizmo_set_variable_av(Av(i))
      call prizmo_set_variable_g0(Gnot)

      call prizmo_set_variable_NH2_incoming(NH2)
      call prizmo_set_variable_NCO_incoming(NCO)

      call prizmo_set_variable_NH2_escaping(NH2)
      call prizmo_set_variable_NCO_escaping(NCO)
      call prizmo_set_variable_NH2O_escaping(NH2O)

      ! do chemistry
      call prizmo(xall(:, i), Tgas(i), jflux(:), dt)

      ! store some data for ouput
      k(:) = prizmo_get_rates(xall(:, i), Tgas(i), jflux(:))
      Tdust = prizmo_get_Tdust(xall(:, i), Tgas(i), jflux(:))
      fluxes(:) = prizmo_get_fluxes(xall(:, i), Tgas(i), jflux(:))

      ! save data to disk
      write(unit, '(9999E17.8e3)') t/spy, z(i)/pc2cm, Tgas(i), Tdust, Ntot, NH2, NCO, Av(i), crflux, &
      k(prizmo_idx_H2_photodissociation), &
      k(prizmo_idx_CO_photodissociation), &
      k(prizmo_idx_C_photoionization), &
      xall(:, i), fluxes(:)


      write(unit_cooling, '(999E17.8e3)') t/spy, z(i)/pc2cm, Av(i), &
        prizmo_get_cooling_array(xall(:, i), Tgas(i), Tdust, jflux(:)), &
        prizmo_get_heating_array(xall(:, i), Tgas(i), Tdust, jflux(:))

      flush(unit)
      flush(unit_cooling)

      zold = z(i)
      xold(:) = xall(:, i)
    end do
    t = t + dt
    if(t > tend) exit
  end do
  close(unit)

  ! say goodbye
  print *, "EVERYTHING DONE!"

end program test

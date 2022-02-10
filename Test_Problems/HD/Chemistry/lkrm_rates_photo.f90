module lkrm_rates_photo
contains

  ! **************************
  subroutine compute_rates_photo(n, Tgas, Tdust, jflux)
    use lkrm_commons, only: kall, nmols, nphoto
    implicit none
    real*8,intent(in)::n(nmols)
    real*8,intent(in)::Tgas, Tdust
    real*8,intent(in)::jflux(nphoto)

    !!BEGIN_RATES_PHOTO
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_RATES_PHOTO

  end subroutine compute_rates_photo

  ! *************************
  ! compute rate reaction, 1/s
  function integrate_rate_photo(jflux, idx) result(k)
    use lkrm_commons, only: nphoto, xsecs, &
         energy_grid, inv_energy_grid, hplanck_eV
    implicit none
    real*8,intent(in)::jflux(nphoto)
    integer,intent(in)::idx
    integer::i
    real*8::k, f0, f1, x0, x1

    k = 0d0  ! eV
    f0 = xsecs(1, idx) * jflux(1) * inv_energy_grid(1)
    x0 = energy_grid(1)  ! eV
    do i=2,nphoto
       f1 = xsecs(i, idx) * jflux(i) * inv_energy_grid(i)
       x1 = energy_grid(i)
       k = k + 0.5 * (f0 + f1) * (x1 - x0)
       f0 = f1
       x0 = x1
    end do

    ! convert to 1/s
    k = k / hplanck_eV

  end function integrate_rate_photo


end module lkrm_rates_photo

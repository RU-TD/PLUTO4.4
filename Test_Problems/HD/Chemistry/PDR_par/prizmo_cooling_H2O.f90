module prizmo_cooling_H2O
  !!BEGIN_COOLING_H2O_VARS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    integer,parameter::nx = 40
    integer,parameter::ny = 40
    integer,parameter::nz = 40

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_COOLING_H2O_VARS

  ! H2O cooling variables for interpolation
  real*8::cooling_H2O_logT(nx)
  real*8::cooling_H2O_logHnuclei(ny)
  real*8::cooling_H2O_logNH2O(nz)
  real*8::cooling_H2O_logCool(nx, ny, nz)
  real*8::cooling_H2O_x0min, cooling_H2O_x1min, cooling_H2O_x2min
  real*8::cooling_H2O_x0pre, cooling_H2O_x1pre, cooling_H2O_x2pre
  real*8::cooling_H2O_ixd0, cooling_H2O_ixd1, cooling_H2O_ixd2
  real*8::xmin, ymin, zmin
  real*8::xfact, yfact, zfact
  real*8::invdx, invdy, invdz
contains

  ! ***********************
  ! get H2O cooling, erg/cm3/s
  function cooling_H2O(n, Tgas, NH2O) result(cool)
    use prizmo_commons
    use prizmo_utils
    use prizmo_fit
    implicit none
    real*8,intent(in)::n(nmols), Tgas, NH2O
    real*8::cool, invT, Htot
    real*8::v0, v1, v2

    cool = 0d0  ! default cooling, erg/cm3/s

    Htot = get_Hnuclei(n(:))
    Htot = max(Htot, 1d-40)

    ! log of fit variables
    v0 = log10(Tgas)
    v1 = log10(Htot)
    v2 = log10(NH2O + 1d-40)

    ! check limits on Tgas (small temperature -> inefficient cooling)
    if(v0 .le. xmin) then
       return
    end if
    v0 = min(v0, cooling_H2O_logT(nx) - 1d-5)
    v0 = max(v0, xmin)

    ! check limits on Hnuclei (large gas density -> inefficient cooling)
    if(v1 .ge. cooling_H2O_logHnuclei(ny)) then
       return
    end if
    v1 = max(v1, cooling_H2O_logHnuclei(1))

    ! check limits on NCO (large H2O column density -> inefficient cooling)
    if(v2 .ge. cooling_H2O_logNH2O(nz)) then
       return
    end if
    v2 = max(v2, cooling_H2O_logNH2O(1))

    cool = 1d1**interp_3d(v0, v1, v2, cooling_H2O_logCool(:, :, :), &
       xmin, xfact, invdx, &
       ymin, yfact, invdy, &
       zmin, zfact, invdz) * Htot * max(n(idx_H2O), 0d0)

    ! avoid constant cooling at lower temperatures
    !if(Tgas < 1e1) then
    !   cool = cool * (tanh(Tgas - 3e0) + 1e0) / 2e0
    !end if

  end function cooling_H2O

  ! ***********************
  ! load H2O cooling from data:
  ! log10(Tgas), log10(nHnuclei/cm3), log10(NH2O/cm2), log10(cooling/erg/cm3/s)
  subroutine load_cooling_H2O()
    use prizmo_commons
    use prizmo_fit
    implicit none
    integer::i, j, k, unit, ios
    real*8::rout(4), dx, dx_old

    call load_data_3d(runtime_folder//"cooling_H2O.dat", cooling_H2O_logCool(:,:,:), &
       cooling_H2O_logT(:), nx, xmin, xfact, invdx, &
       cooling_H2O_logHnuclei(:), ny, ymin, yfact, invdy, &
       cooling_H2O_logNH2O(:), nz, zmin, zfact, invdz)

  end subroutine load_cooling_H2O

end module prizmo_cooling_H2O

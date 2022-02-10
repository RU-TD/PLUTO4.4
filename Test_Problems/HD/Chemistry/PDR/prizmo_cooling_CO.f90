module prizmo_cooling_CO
  !!BEGIN_COOLING_CO_VARS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    integer,parameter::nx = 40
    integer,parameter::ny = 40
    integer,parameter::nz = 40

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_COOLING_CO_VARS

  ! CO cooling variables for interpolation
  real*8::cooling_CO_logT(nx)
  real*8::cooling_CO_logHnuclei(ny)
  real*8::cooling_CO_logNCO(nz)
  real*8::cooling_CO_logCool(nx, ny, nz)
  real*8::cooling_CO_x0min, cooling_CO_x1min, cooling_CO_x2min
  real*8::cooling_CO_x0pre, cooling_CO_x1pre, cooling_CO_x2pre
  real*8::cooling_CO_ixd0, cooling_CO_ixd1, cooling_CO_ixd2

  real*8::xmin, ymin, zmin
  real*8::xfact, yfact, zfact
  real*8::invdx, invdy, invdz
contains

  ! ***********************
  ! get CO cooling, erg/cm3/s
  function cooling_CO(n, Tgas, NCO) result(cool)
    use prizmo_commons
    use prizmo_utils
    use prizmo_fit
    implicit none
    real*8,intent(in)::n(nmols), Tgas, NCO
    real*8::cool
    real*8::v0, v1, v2, Htot

    cool = 0d0  ! default cooling, erg/cm3/s

    Htot = get_Hnuclei(n(:))
    Htot = max(Htot, 1d-40)

    ! log of fit variables
    v0 = log10(Tgas)
    v1 = log10(Htot)
    v2 = log10(NCO + 1d-40)

    ! check limits on Tgas (small temperature -> inefficient cooling)
    if(v0 .le. xmin) then
       return
    end if
    v0 = min(v0, cooling_CO_logT(nx) - 1d-5)
    v0 = max(v0, xmin)

    ! check limits on Hnuclei (large gas density -> inefficient cooling)
    if(v1 .ge. cooling_CO_logHnuclei(ny)) then
       return
    end if
    v1 = max(v1, cooling_CO_logHnuclei(1))

    ! check limits on NCO (large CO column density -> inefficient cooling)
    if(v2 .ge. cooling_CO_logNCO(nz)) then
       return
    end if
    v2 = max(v2, cooling_CO_logNCO(1))

    if(n(idx_CO) < -1d-20) then
      print *, "WARNING: negative CO in CO cooling", n(idx_CO)
    end if

    cool = 1d1**interp_3d(v0, v1, v2, cooling_CO_logCool, &
       xmin, xfact, invdx, &
       ymin, yfact, invdy, &
       zmin, zfact, invdz)  * Htot * max(n(idx_CO), 0d0)

  end function cooling_CO

  ! ***********************
  ! load CO cooling from data:
  ! log10(Tgas), log10(nHnuclei/cm3), log10(NCO/cm2), log10(cooling/erg/cm3/s)
  subroutine load_cooling_CO()
    use prizmo_commons
    use prizmo_fit
    implicit none
    character(len=30)::fname

    fname = "cooling_CO.dat"

    call load_data_3d(runtime_folder//trim(fname), cooling_CO_logCool(:,:,:), &
       cooling_CO_logT(:), nx, xmin, xfact, invdx, &
       cooling_CO_logHnuclei(:), ny, ymin, yfact, invdy, &
       cooling_CO_logNCO(:), nz, zmin, zfact, invdz)

    print *, "CO cooling loaded from file "//trim(fname)//" Trange:", &
      1d1**cooling_CO_logT(1),1d1**cooling_CO_logT(nx)

  end subroutine load_cooling_CO

end module prizmo_cooling_CO
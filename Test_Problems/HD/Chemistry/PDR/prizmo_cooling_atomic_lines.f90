module prizmo_cooling_atomic_lines

contains

  ! **********************
  ! get atomic cooling, erg/cm3/s
  ! NOTE: this is no longer used, see emission subroutine instead
  function cooling_atomic_lines(n, Tgas) result(cool)
    use prizmo_emission
    use prizmo_commons
    real*8,intent(in)::n(nmols), Tgas
    real*8::cool

    cool = get_atomic_cooling(n(:), Tgas)

  end function cooling_atomic_lines

  ! **********************
  ! get atomic cooling, erg/cm3/s
  function cooling_atomic_lines_old(n, inTgas) result(cool)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    real*8,intent(in)::n(nmols), inTgas
    real*8::cool, invT, lnT, Tgas
    integer::i

    !!BEGIN_COOLING_ATOMIC_LINES_DEFINITIONS

    !!END_COOLING_ATOMIC_LINES_DEFINITIONS

    ! define local Tgas value and upper limit
    Tgas = min(inTgas, Tgas)

    invT = 1d0 / Tgas
    lnT = log(Tgas)

    ! default cooling, erg/cm3/s
    cool = 0d0

    !!BEGIN_COOLING_ATOMIC_LINES

    !!END_COOLING_ATOMIC_LINES

    if(cool < 0d0) then
       print *, "ERROR: atomic cooling negative", cool
       stop
    end if

  end function cooling_atomic_lines_old

end module prizmo_cooling_atomic_lines

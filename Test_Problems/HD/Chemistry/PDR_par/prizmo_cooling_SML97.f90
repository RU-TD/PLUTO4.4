module prizmo_cooling_SML97
contains

  ! ***********************
  ! Smith+McLow 1997 cooling, erg/cm3/s
  ! https://ui.adsabs.harvard.edu/abs/1997A&A...326..801S/abstract
  function cooling_SML97(n, Tgas) result(cool)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::cool, ntot

    ntot = get_ntot(n(:))

    cool = 4.2d-31 * ntot * Tgas**3.3
    cool = max(cool, 0d0)

  end function cooling_SML97

end module prizmo_cooling_SML97

module prizmo_heating
  implicit none
  integer,parameter::heating_number=5
contains

  ! **********************
  function heating(n, Tgas, Tdust, jflux) result(heat)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::heat, heat_array(heating_number), ntot
    character(len=max_character_len)::names(nmols)
    integer::i

    heat_array(:) = get_heating_array(n(:), Tgas, Tdust, jflux(:))

    !ntot = get_ntot(n)
    !if(maxval(abs(heat_array) / ntot / kboltzmann) > 1d0) then
    !   print *, heat_array(:) / ntot / kboltzmann
    !   print *, "ERROR: large heating!"
    !   stop
    !end if

    heat = sum(heat_array)

    ! print '(a10,99E17.8e3)', "heat", heat_array(:), heat
    ! print '(a10,99E17.8e3)', "Tgas, heat", Tgas, heat

    if(minval(heat_array) < 0d0) then
      print *, "ERROR: negative heating!"
      print *, "Tgas, Tdust", Tgas, Tdust
      print '(a10,99E17.8e3)', "CR", heat_array(1)
      print '(a10,99E17.8e3)', "phe", heat_array(2)
      print '(a10,99E17.8e3)', "photo", heat_array(3)
      print '(a10,99E17.8e3)', "PAH", heat_array(4)
      print '(a10,99E17.8e3)', "chem", heat_array(5)
      names(:) = get_species_names()
      do i=1, nmols
        print *, i, names(i), n(i)
      end do
      stop
    end	if

  end function heating

  ! **********************
  function get_heating_array(n, Tgas, Tdust, jflux) result(heat_array)
    use prizmo_commons
    use prizmo_heating_CR
    use prizmo_heating_photoelectric
    use prizmo_heating_photo
    use prizmo_heating_PAH
    use prizmo_heating_chemical
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::heat_array(heating_number)

    if(heating_mode == 0) then
      heat_array(1) = heating_CR(n(:), Tgas)
      !heat_array(2) = heating_photoelectric(n(:), Tgas, variable_G0, variable_Av)
      heat_array(2) = heating_photoelectric_multi(n(:), Tgas, Tdust, jflux(:))
      heat_array(3) = heating_photo(n(:), Tgas, jflux(:))
      heat_array(4) = heating_PAH(n(:), Tgas, variable_G0, variable_Av)
      heat_array(5) = heating_chemical(n(:), Tgas)
    else if(heating_mode == 1) then
      heat_array(:) = 0d0
      heat_array(1) = heating_CR(n(:), Tgas)
    else if(heating_mode == -1) then
      heat_array(:) = 0d0
    else
      print *, "ERROR: unknown heating mode", heating_mode
      stop
    end if

  end function get_heating_array

end module prizmo_heating

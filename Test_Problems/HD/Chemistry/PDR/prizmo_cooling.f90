module prizmo_cooling
  implicit none
  integer,parameter::cooling_number=7

contains

  ! *************************
  ! compute cooling in erg/cm3/s
  function cooling(n, Tgas, Tdust, jflux) result(cool)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::cool, cool_array(cooling_number), ntot
    character(len=max_character_len)::names(nmols)
    integer::i

    cool_array(:) = get_cooling_array(n(:), Tgas, Tdust, jflux(:))

    ! ntot = get_ntot(n)
    ! if(maxval(abs(cool_array) / ntot / kboltzmann) > 1d0) then
    !    print *, cool_array(:) / ntot / kboltzmann
    !    print *, "ERROR: large cooling!"
    !    stop
    ! end if

    cool = sum(cool_array)  ! * (tanh(Tgas - 1d1) + 1d0) / 2d0

    ! check if any of the cooling function is negative
    ! (exclude dust cooling that can be negative)
    if(minval(cool_array(1:4)) < 0d0 .or. minval(cool_array(6:)) < 0d0) then
      print *, "ERROR: negative cooling!"
      print *, "Tgas, Tdust", Tgas, Tdust
      print '(a10,99E17.8e3)', "H2", cool_array(1)
      print '(a10,99E17.8e3)', "CO", cool_array(2)
      print '(a10,99E17.8e3)', "atomic", cool_array(3)
      print '(a10,99E17.8e3)', "chem", cool_array(4)
      print '(a10,99E17.8e3)', "dust", cool_array(5)
      print '(a10,99E17.8e3)', "H2O", cool_array(6)
      print '(a10,99E17.8e3)', "BS", cool_array(7)
      names(:) = get_species_names()
      do i=1, nmols
        print *, i, names(i), n(i)
      end do
      stop
    end	if


    ! print '(a10,99E17.8e3)', "cool", cool_array(:), cool
    ! print '(a10,99E17.8e3)', "Tgas, cool", Tgas, cool

  end function cooling

  ! ***********************
  function get_cooling_array(n, Tgas, Tdust, jflux) result(cool_array)
    use prizmo_commons
    use prizmo_cooling_H2
    use prizmo_cooling_CO
    use prizmo_cooling_H2O
    use prizmo_cooling_atomic_lines
    use prizmo_cooling_chemical
    use prizmo_cooling_dust
    use prizmo_cooling_bremsstrahlung
    use prizmo_cooling_SML97
    implicit none
    real*8,intent(in)::n(nmols), Tgas, Tdust, jflux(nphoto)
    real*8::cool_array(cooling_number)

    if(cooling_mode == 0) then
      cool_array(1) = cooling_H2(n(:), Tgas)
      cool_array(2) = cooling_CO(n(:), Tgas, variable_NCO_escaping)
      cool_array(3) = cooling_atomic_lines(n(:), Tgas)
      cool_array(4) = cooling_chemical(n(:), Tgas)
      cool_array(5) = cooling_dust(n(:), Tgas, Tdust, jflux(:))
      cool_array(6) = cooling_H2O(n(:), Tgas, variable_NH2O_escaping)
      cool_array(7) = cooling_bremsstrahlung(n(:), Tgas)
    else if(cooling_mode == 1) then
      cool_array(:) = 0d0
      cool_array(1) = cooling_SML97(n(:), Tgas)
    else if(cooling_mode == -1) then
      cool_array(:) = 0d0
    else
      print *, "ERROR: unknown cooling mode", heating_mode
      stop
    end if

  end function get_cooling_array

end module prizmo_cooling

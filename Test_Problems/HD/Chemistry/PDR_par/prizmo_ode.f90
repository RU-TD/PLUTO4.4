module prizmo_ode
  use prizmo_commons
  real*8::jflux_arg(nphoto)
contains

  !************************
  !evolve chemistry for a time-step dt (s)
  ! n(:) are species number densities
  subroutine dochem(n, Tgas, jflux, dt)
    use prizmo_commons
    use prizmo_rates
    use prizmo_utils
    use prizmo_emission
    use prizmo_rates_evaluate_once
    use prizmo_cooling
    use prizmo_heating
    use prizmo_tdust
    implicit none
    integer,parameter::neq=nmols+1
    real*8,intent(inout)::n(nmols), Tgas
    real*8,intent(in)::dt, jflux(nphoto)
    real*8::y(neq), emission(nphoto), y_old(neq), Tdust_old
    real*8::cooling_array(cooling_number), heating_array(heating_number)
    integer::istep

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer,parameter::lrw=20+3*neq**2+20*neq
    integer,parameter::liw=30!+nmols+nmols**2

    !tolerances
    real*8::rtol(neq)
    real*8::atol(neq)

    integer::neqa(1),itol,itask,istate,iopt,mf
    integer::iwork(liw), i
    real*8::rwork(lrw), tloc, Tdust

    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are arrays

    !number of equations is a single sized array
    neqa(:) = neq

    !  FROM DLSODES manual
    !  Name    Location   Meaning and default value
    !  ------  ---------  -----------------------------------------------
    !  H0      RWORK(5)   Step size to be attempted on the first step.
    !                     The default value is determined by the solver.
    !  HMAX    RWORK(6)   Maximum absolute step size allowed.  The
    !                     default value is infinite.
    !  HMIN    RWORK(7)   Minimum absolute step size allowed.  The
    !                     default value is 0.  (This lower bound is not
    !                     enforced on the final step before reaching
    !                     TCRIT when ITASK = 4 or 5.)
    !  MAXORD  IWORK(5)   Maximum order to be allowed.  The default value
    !                     is 12 if METH = 1, and 5 if METH = 2. (See the
    !                     MF description above for METH.)  If MAXORD
    !                     exceeds the default value, it will be reduced
    !                     to the default value.  If MAXORD is changed
    !                     during the problem, it may cause the current
    !                     order to be reduced.
    !  MXSTEP  IWORK(6)   Maximum number of (internally defined) steps
    !                     allowed during one call to the solver.  The
    !                     default value is 500.
    !  MXHNIL  IWORK(7)   Maximum number of messages printed (per
    !                     problem) warning that T + H = T on a step
    !                     (H = step size).  This must be positive to
    !                     result in a nondefault value.  The default
    !                     value is 10.


    itask = 1
    iopt = 1
    iwork(6) = int(1e6) !maximum number of iteration before warning
    iwork(5) = 2 !maximum integration order, zero means default
    ! 222: JAC and IA/JA both generated
    ! 121: IA/JA generated, JAC provided
    !  21: JAC and IA/JA both provided
    MF = 222

    ! tolerances
    atol(1:nmols) = 1d-12 !1d-12
    atol(idx_Tgas) = 1d-12
    rtol(1:nmols) = 1d-6
    rtol(idx_Tgas) = 1d-6

    ! copy variables to dlosdes varibale array
    y(1:nmols) = n(:)
    y(idx_Tgas) = Tgas

    ! evaluate rates that do not depend on Tgas or ntot
    call init_evaluate_once()

    ! use common variable to pass to fex
    jflux_arg(:) = jflux(:)

    emission(:) = 0d0

    !!BEGIN_SPARSITY

    !!END_SPARSITY

    ! store input variables to compute cooling/heating in case of problems
    y_old(:) = y(:)
    Tdust_old = get_Tdust(n(:), Tgas, jflux_arg(:))

    istate = 1
    tloc = 0d0
    istep = 0
    do
       istep = istep + 1
       !call solver, see DLSODES documentation
       CALL DLSODES(fex, neqa(:), y(:), tloc, dt, &
            ITOL, RTOL, ATOL, ITASK, ISTATE, IOPT, RWORK, LRW, IWORK, &
            LIW, JES, MF)
       !check solver output state
       if(istate==-1) then
          !maximum number of iteration reached, continue
          istate = 1
          cycle
       elseif(istate==2) then
          !success integration
          exit
       elseif(istate==-5 .or. istate==-4) then
          !problem with sparsity, need to recompute
          istate = 3
          cycle
       else
          !unknonw problem stop program
          print *, "unknown istate ", istate
          print *, "stopping..."
          stop
       end if
    end do

    ! copy dlsodes variables back to physical variables
    n(:) = y(1:nmols)
    Tgas = y(idx_Tgas)

    ! limit ouput temperature with user limits
    ! default range is [-1e99, 1e99] K, i.e. no limits
    Tgas = min(max(Tgas, Tgas_min_forced), Tgas_max_forced)

    ! check negative temperature, if so writes negative species
    ! and stop program execution
    if(Tgas < 0d0) then
       print *, "ERROR: negative temperature (prizmo_ode.f90)!", Tgas
       print *, "input Tgas, Tdust:", y_old(idx_Tgas), Tdust_old
       print *, "calls to DLSODES:", istep
       do i=1,nmols
          if(n(i) < 0d0) then
             print *, "negative species:", i, n(i)
          end if
       end do
       cooling_array(:) = get_cooling_array(y_old(1:nmols), y_old(idx_Tgas), Tdust_old, jflux(:))
       heating_array(:) = get_heating_array(y_old(1:nmols), y_old(idx_Tgas), Tdust_old, jflux(:))
       do i=1, cooling_number
         print *, "cooling", i, cooling_array(i)
       end do
       print *, "cooling tot", sum(cooling_array)
       do i=1, heating_number
         print *, "heating", i, heating_array(i)
       end do
       print *, "heating tot", sum(heating_array)
       stop
    end if

    ! check negative species
    do i=1,nmols
       if(n(i) < 0d0) then
          if(abs(n(i)) > atol(i) / 1d1) then
             print *, "WARNING: negative species (prizmo_ode.f90)!"
             print *, i, n(i), get_species_name(i)
             ! stop
          end if
       end if
       n(i) = max(n(i), 0d0)
    end do

    ! compute electrons to ensure neutrality
    n(idx_E) = get_electrons(n(1:nmols))

    ! add emissions from current chemical composition
    ! to the global emission array
    !call add_emission(emission(:), n(:), Tgas)

  end subroutine dochem

  ! ********************
  function get_timescale(n, Tgas) result(dt)
    use prizmo_commons
    use prizmo_rates
    use prizmo_rates_evaluate_once
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::y(nmols+1), dy(nmols+1), tt, dt
    integer,parameter::neq=nmols+1
    integer::i

    tt = 0d0
    y(1:nmols) = n(:)
    y(idx_Tgas) = Tgas

    call init_evaluate_once()

    call fex(neq, tt, y(:), dy(:))

    dt = 1d99
    do i=1,neq
       if(abs(y(i) * dy(i)) < 1d-40) cycle
       dt = min(dt, y(i) / abs(dy(i)))
    end do

  end function get_timescale

  !*********************
  !differential equations, returns dn(:)
  ! see DLSODES documentation
  subroutine fex(neq, tt, y, dy)
    use prizmo_commons
    use prizmo_cooling
    use prizmo_heating
    use prizmo_utils
    use prizmo_tdust
    use prizmo_self_shielding
    use prizmo_rates
    implicit none
    integer::neq, i
    real*8::y(neq), tt, dy(neq), fluxes(nrea)
    real*8::ntot, gamma_ad, Tgas, Tdust, n(nmols), dn(nmols)

    !if(y(idx_Tgas) < 0d0) then
    !   print *, "T<0", y(idx_Tgas)
    !end if

    ! get temperature as local variable, K, check ranges
    ! this check is to avoid strange unphysical behaviour when the solver
    ! computes jacobian internally
    Tgas = min(max(y(idx_Tgas), Tgas_min), Tgas_max)

    ! y are the solver variables, that includes temperature
    ! while n are just the chemical species (used for the sake of clarity)
    n(:) = y(1:nmols)

    ! get adiabatic index from utility function
    gamma_ad = get_gamma_ad(n(:), Tgas)

    ! get total density, cm-3
    ntot = get_ntot(n(:))

    ! use the same temperature for dust and gas
    Tdust = get_Tdust(n(:), Tgas, jflux_arg(:))

    ! compute non-photochemical rates
    call compute_rates(n(:), Tgas, Tdust, jflux_arg(:))

    ! temperature ODE
    dy(idx_Tgas) = (gamma_ad - 1d0) * (heating(n(:), Tgas, Tdust, jflux_arg(:)) - &
        cooling(n(:), Tgas, Tdust, jflux_arg(:))) / ntot / kboltzmann

    !!BEGIN_FLUXES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    fluxes(1) = n(idx_H)*n(idx_CH)
    fluxes(2) = n(idx_H)*n(idx_CH3)
    fluxes(3) = n(idx_H)*n(idx_CH4)
    fluxes(4) = n(idx_H)*n(idx_OH)
    fluxes(5) = n(idx_H)*n(idx_H2O)
    fluxes(6) = n(idx_H)*n(idx_CO)
    fluxes(7) = n(idx_H)*n(idx_O2)
    fluxes(8) = n(idx_H2)*n(idx_C)
    fluxes(9) = n(idx_H2)*n(idx_CH)
    fluxes(10) = n(idx_H2)*n(idx_CH2)
    fluxes(11) = n(idx_H2)*n(idx_CH3)
    fluxes(12) = n(idx_H2)*n(idx_O)
    fluxes(13) = n(idx_H2)*n(idx_OH)
    fluxes(14) = n(idx_H2)*n(idx_O2)
    fluxes(15) = n(idx_C)*n(idx_CH2)
    fluxes(16) = n(idx_C)*n(idx_OH)
    fluxes(17) = n(idx_C)*n(idx_OH)
    fluxes(18) = n(idx_CH)*n(idx_O)
    fluxes(19) = n(idx_CH)*n(idx_CH4)
    fluxes(20) = n(idx_CH2)*n(idx_CH2)
    fluxes(21) = n(idx_CH2)*n(idx_O)
    fluxes(22) = n(idx_CH2)*n(idx_CH4)
    fluxes(23) = n(idx_CH2)*n(idx_OH)
    fluxes(24) = n(idx_CH2)*n(idx_OH)
    fluxes(25) = n(idx_CH2)*n(idx_O2)
    fluxes(26) = n(idx_CH3)*n(idx_CH3)
    fluxes(27) = n(idx_CH3)*n(idx_OH)
    fluxes(28) = n(idx_CH3)*n(idx_OH)
    fluxes(29) = n(idx_CH3)*n(idx_H2O)
    fluxes(30) = n(idx_O)*n(idx_CH4)
    fluxes(31) = n(idx_O)*n(idx_OH)
    fluxes(32) = n(idx_O)*n(idx_H2O)
    fluxes(33) = n(idx_CH4)*n(idx_OH)
    fluxes(34) = n(idx_OH)*n(idx_OH)
    fluxes(35) = n(idx_H)*n(idx_CH2j)
    fluxes(36) = n(idx_H)*n(idx_CH3j)
    fluxes(37) = n(idx_H2)*n(idx_Hej)
    fluxes(38) = n(idx_H2)*n(idx_Cj)
    fluxes(39) = n(idx_Hej)*n(idx_CO)
    fluxes(40) = n(idx_Hj)*n(idx_H2)
    fluxes(41) = n(idx_H)*n(idx_Hej)
    fluxes(42) = n(idx_H)*n(idx_Oj)
    fluxes(43) = n(idx_Hj)*n(idx_O)
    fluxes(44) = n(idx_Hej)*n(idx_C)
    fluxes(45) = n(idx_Oj)*n(idx_CO)
    fluxes(46) = n(idx_H2j)*n(idx_E)
    fluxes(47) = n(idx_H3j)*n(idx_E)
    fluxes(48) = n(idx_H3j)*n(idx_E)
    fluxes(49) = n(idx_CHj)*n(idx_E)
    fluxes(50) = n(idx_CH2j)*n(idx_E)
    fluxes(51) = n(idx_CH2j)*n(idx_E)
    fluxes(52) = n(idx_CH2j)*n(idx_E)
    fluxes(53) = n(idx_CH3j)*n(idx_E)
    fluxes(54) = n(idx_CH3j)*n(idx_E)
    fluxes(55) = n(idx_CH3j)*n(idx_E)
    fluxes(56) = n(idx_CH4j)*n(idx_E)
    fluxes(57) = n(idx_CH4j)*n(idx_E)
    fluxes(58) = n(idx_OHj)*n(idx_E)
    fluxes(59) = n(idx_CH5j)*n(idx_E)
    fluxes(60) = n(idx_CH5j)*n(idx_E)
    fluxes(61) = n(idx_H2Oj)*n(idx_E)
    fluxes(62) = n(idx_H2Oj)*n(idx_E)
    fluxes(63) = n(idx_H2Oj)*n(idx_E)
    fluxes(64) = n(idx_H3Oj)*n(idx_E)
    fluxes(65) = n(idx_H3Oj)*n(idx_E)
    fluxes(66) = n(idx_H3Oj)*n(idx_E)
    fluxes(67) = n(idx_H3Oj)*n(idx_E)
    fluxes(68) = n(idx_COj)*n(idx_E)
    fluxes(69) = n(idx_HCOj)*n(idx_E)
    fluxes(70) = n(idx_O2j)*n(idx_E)
    fluxes(71) = n(idx_Hj)*n(idx_E)
    fluxes(72) = n(idx_Hej)*n(idx_E)
    fluxes(73) = n(idx_Cj)*n(idx_E)
    fluxes(74) = n(idx_CH3j)*n(idx_E)
    fluxes(75) = n(idx_Oj)*n(idx_E)
    fluxes(76) = n(idx_Hj)*n(idx_H)
    fluxes(77) = n(idx_H)*n(idx_O)
    fluxes(78) = n(idx_H)*n(idx_OH)
    fluxes(79) = n(idx_H2)*n(idx_Cj)
    fluxes(80) = n(idx_H2)*n(idx_CH)
    fluxes(81) = n(idx_H2)*n(idx_CH3j)
    fluxes(82) = n(idx_O)*n(idx_O)
    fluxes(83) = n(idx_H)*n(idx_H2)
    fluxes(84) = n(idx_H)*n(idx_CH)
    fluxes(85) = n(idx_H)*n(idx_OH)
    fluxes(86) = n(idx_H)*n(idx_H2O)
    fluxes(87) = n(idx_H)*n(idx_O2)
    fluxes(88) = n(idx_H2)*n(idx_E)
    fluxes(89) = n(idx_H2)*n(idx_H2)
    fluxes(90) = n(idx_H2)*n(idx_CH)
    fluxes(91) = n(idx_H2)*n(idx_OH)
    fluxes(92) = n(idx_H2)*n(idx_H2O)
    fluxes(93) = n(idx_H2)*n(idx_O2)
    fluxes(94) = n(idx_CH)*n(idx_O)
    fluxes(95) = n(idx_H)*n(idx_CH2)
    fluxes(96) = n(idx_C)*n(idx_O2)
    fluxes(97) = n(idx_CH)*n(idx_O)
    fluxes(98) = n(idx_CH)*n(idx_O2)
    fluxes(99) = n(idx_CH2)*n(idx_O)
    fluxes(100) = n(idx_CH2)*n(idx_O)
    fluxes(101) = n(idx_H)*n(idx_CHj)
    fluxes(102) = n(idx_H2)*n(idx_CHj)
    fluxes(103) = n(idx_Cj)*n(idx_OH)
    fluxes(104) = n(idx_Hj)*n(idx_H2O)
    fluxes(105) = n(idx_Hj)*n(idx_O2)
    fluxes(106) = n(idx_H2)*n(idx_Hej)
    fluxes(107) = n(idx_H)*n(idx_C)
    fluxes(108) = n(idx_H)*n(idx_Cj)
    fluxes(109) = n(idx_H2)*n(idx_C)
    fluxes(110) = n(idx_C)*n(idx_O)
    fluxes(111) = n(idx_Cj)*n(idx_O)
    fluxes(112) = n(idx_Hj)*n(idx_CH2)
    fluxes(113) = n(idx_Hj)*n(idx_CH4)
    fluxes(114) = n(idx_H)*n(idx_CH4j)
    fluxes(115) = n(idx_H)*n(idx_CH5j)
    fluxes(116) = n(idx_H2j)*n(idx_H2)
    fluxes(117) = n(idx_H2j)*n(idx_C)
    fluxes(118) = n(idx_H2j)*n(idx_CH)
    fluxes(119) = n(idx_H2j)*n(idx_CH2)
    fluxes(120) = n(idx_H2)*n(idx_CH2j)
    fluxes(121) = n(idx_H2j)*n(idx_O)
    fluxes(122) = n(idx_H2)*n(idx_Oj)
    fluxes(123) = n(idx_H2j)*n(idx_CH4)
    fluxes(124) = n(idx_H2)*n(idx_CH4j)
    fluxes(125) = n(idx_H2j)*n(idx_CH4)
    fluxes(126) = n(idx_H2j)*n(idx_OH)
    fluxes(127) = n(idx_H2)*n(idx_OHj)
    fluxes(128) = n(idx_H2j)*n(idx_H2O)
    fluxes(129) = n(idx_H2)*n(idx_H2Oj)
    fluxes(130) = n(idx_H2j)*n(idx_CO)
    fluxes(131) = n(idx_H2)*n(idx_COj)
    fluxes(132) = n(idx_H3j)*n(idx_C)
    fluxes(133) = n(idx_H3j)*n(idx_CH)
    fluxes(134) = n(idx_H3j)*n(idx_CH2)
    fluxes(135) = n(idx_H3j)*n(idx_CH3)
    fluxes(136) = n(idx_H3j)*n(idx_O)
    fluxes(137) = n(idx_H3j)*n(idx_CH4)
    fluxes(138) = n(idx_H3j)*n(idx_OH)
    fluxes(139) = n(idx_H3j)*n(idx_H2O)
    fluxes(140) = n(idx_H3j)*n(idx_CO)
    fluxes(141) = n(idx_Hej)*n(idx_CH)
    fluxes(142) = n(idx_Hej)*n(idx_CH2)
    fluxes(143) = n(idx_Hej)*n(idx_CH2)
    fluxes(144) = n(idx_Hej)*n(idx_CH3)
    fluxes(145) = n(idx_Hej)*n(idx_CH4)
    fluxes(146) = n(idx_Hej)*n(idx_CH4)
    fluxes(147) = n(idx_Hej)*n(idx_CH4)
    fluxes(148) = n(idx_Hej)*n(idx_CH4)
    fluxes(149) = n(idx_Hej)*n(idx_OH)
    fluxes(150) = n(idx_Hej)*n(idx_H2O)
    fluxes(151) = n(idx_Hej)*n(idx_H2O)
    fluxes(152) = n(idx_Hej)*n(idx_CO)
    fluxes(153) = n(idx_Hej)*n(idx_O2)
    fluxes(154) = n(idx_C)*n(idx_OHj)
    fluxes(155) = n(idx_Cj)*n(idx_OH)
    fluxes(156) = n(idx_C)*n(idx_CH5j)
    fluxes(157) = n(idx_C)*n(idx_H2Oj)
    fluxes(158) = n(idx_Cj)*n(idx_H2O)
    fluxes(159) = n(idx_C)*n(idx_H3Oj)
    fluxes(160) = n(idx_C)*n(idx_HCOj)
    fluxes(161) = n(idx_Cj)*n(idx_O2)
    fluxes(162) = n(idx_Cj)*n(idx_O2)
    fluxes(163) = n(idx_C)*n(idx_O2j)
    fluxes(164) = n(idx_CHj)*n(idx_O)
    fluxes(165) = n(idx_CH)*n(idx_Oj)
    fluxes(166) = n(idx_CH)*n(idx_OHj)
    fluxes(167) = n(idx_CHj)*n(idx_OH)
    fluxes(168) = n(idx_CH)*n(idx_CH5j)
    fluxes(169) = n(idx_CHj)*n(idx_H2O)
    fluxes(170) = n(idx_CH)*n(idx_H2Oj)
    fluxes(171) = n(idx_CHj)*n(idx_H2O)
    fluxes(172) = n(idx_CH)*n(idx_H3Oj)
    fluxes(173) = n(idx_CH)*n(idx_COj)
    fluxes(174) = n(idx_CH)*n(idx_HCOj)
    fluxes(175) = n(idx_CHj)*n(idx_O2)
    fluxes(176) = n(idx_CHj)*n(idx_O2)
    fluxes(177) = n(idx_CH)*n(idx_O2j)
    fluxes(178) = n(idx_CH2j)*n(idx_O)
    fluxes(179) = n(idx_CH2)*n(idx_OHj)
    fluxes(180) = n(idx_CH2)*n(idx_CH5j)
    fluxes(181) = n(idx_CH2)*n(idx_H2Oj)
    fluxes(182) = n(idx_CH2)*n(idx_H3Oj)
    fluxes(183) = n(idx_CH2)*n(idx_COj)
    fluxes(184) = n(idx_CH2)*n(idx_HCOj)
    fluxes(185) = n(idx_CH2j)*n(idx_O2)
    fluxes(186) = n(idx_CH3j)*n(idx_O)
    fluxes(187) = n(idx_Oj)*n(idx_CH4)
    fluxes(188) = n(idx_O)*n(idx_CH4j)
    fluxes(189) = n(idx_Oj)*n(idx_OH)
    fluxes(190) = n(idx_O)*n(idx_OHj)
    fluxes(191) = n(idx_O)*n(idx_CH5j)
    fluxes(192) = n(idx_O)*n(idx_H2Oj)
    fluxes(193) = n(idx_CH4j)*n(idx_CH4)
    fluxes(194) = n(idx_CH4)*n(idx_OHj)
    fluxes(195) = n(idx_CH4)*n(idx_OHj)
    fluxes(196) = n(idx_CH4j)*n(idx_H2O)
    fluxes(197) = n(idx_CH4)*n(idx_H2Oj)
    fluxes(198) = n(idx_CH4j)*n(idx_CO)
    fluxes(199) = n(idx_CH4)*n(idx_COj)
    fluxes(200) = n(idx_OHj)*n(idx_OH)
    fluxes(201) = n(idx_OH)*n(idx_CH5j)
    fluxes(202) = n(idx_OHj)*n(idx_H2O)
    fluxes(203) = n(idx_OH)*n(idx_H2Oj)
    fluxes(204) = n(idx_OHj)*n(idx_CO)
    fluxes(205) = n(idx_OH)*n(idx_COj)
    fluxes(206) = n(idx_OH)*n(idx_HCOj)
    fluxes(207) = n(idx_CH5j)*n(idx_H2O)
    fluxes(208) = n(idx_CH5j)*n(idx_CO)
    fluxes(209) = n(idx_H2Oj)*n(idx_H2O)
    fluxes(210) = n(idx_H2Oj)*n(idx_CO)
    fluxes(211) = n(idx_H2O)*n(idx_COj)
    fluxes(212) = n(idx_H2O)*n(idx_HCOj)
    fluxes(213) = n(idx_H)*n(idx_H2j)
    fluxes(214) = n(idx_Hj)*n(idx_CH)
    fluxes(215) = n(idx_Hj)*n(idx_CH2)
    fluxes(216) = n(idx_Hj)*n(idx_CH3)
    fluxes(217) = n(idx_Hj)*n(idx_CH4)
    fluxes(218) = n(idx_Hj)*n(idx_OH)
    fluxes(219) = n(idx_H)*n(idx_COj)
    fluxes(220) = n(idx_H2j)*n(idx_CH)
    fluxes(221) = n(idx_H2j)*n(idx_CH2)
    fluxes(222) = n(idx_H2j)*n(idx_CH4)
    fluxes(223) = n(idx_H2j)*n(idx_OH)
    fluxes(224) = n(idx_H2j)*n(idx_H2O)
    fluxes(225) = n(idx_H2j)*n(idx_CO)
    fluxes(226) = n(idx_H2j)*n(idx_O2)
    fluxes(227) = n(idx_Hej)*n(idx_CH)
    fluxes(228) = n(idx_Hej)*n(idx_CH4)
    fluxes(229) = n(idx_Hej)*n(idx_H2O)
    fluxes(230) = n(idx_Hej)*n(idx_O2)
    fluxes(231) = n(idx_Cj)*n(idx_CH)
    fluxes(232) = n(idx_Cj)*n(idx_CH2)
    fluxes(233) = n(idx_C)*n(idx_COj)
    fluxes(234) = n(idx_C)*n(idx_O2j)
    fluxes(235) = n(idx_CH)*n(idx_Oj)
    fluxes(236) = n(idx_CH)*n(idx_OHj)
    fluxes(237) = n(idx_CH)*n(idx_H2Oj)
    fluxes(238) = n(idx_CH)*n(idx_COj)
    fluxes(239) = n(idx_CH)*n(idx_O2j)
    fluxes(240) = n(idx_CH2)*n(idx_Oj)
    fluxes(241) = n(idx_CH2)*n(idx_OHj)
    fluxes(242) = n(idx_CH2)*n(idx_H2Oj)
    fluxes(243) = n(idx_CH2)*n(idx_COj)
    fluxes(244) = n(idx_CH2)*n(idx_O2j)
    fluxes(245) = n(idx_Oj)*n(idx_CH4)
    fluxes(246) = n(idx_Oj)*n(idx_OH)
    fluxes(247) = n(idx_Oj)*n(idx_H2O)
    fluxes(248) = n(idx_O)*n(idx_COj)
    fluxes(249) = n(idx_Oj)*n(idx_O2)
    fluxes(250) = n(idx_CH4)*n(idx_COj)
    fluxes(251) = n(idx_CH4j)*n(idx_O2)
    fluxes(252) = n(idx_OHj)*n(idx_H2O)
    fluxes(253) = n(idx_OH)*n(idx_COj)
    fluxes(254) = n(idx_OHj)*n(idx_O2)
    fluxes(255) = n(idx_H2O)*n(idx_COj)
    fluxes(256) = n(idx_H2Oj)*n(idx_O2)
    fluxes(257) = n(idx_COj)*n(idx_O2)
    fluxes(258) = n(idx_H2)
    fluxes(259) = n(idx_CO)
    fluxes(260) = n(idx_H2j)
    fluxes(261) = n(idx_H3j)
    fluxes(262) = n(idx_H3j)
    fluxes(263) = n(idx_C)
    fluxes(264) = n(idx_CH)
    fluxes(265) = n(idx_CH)
    fluxes(266) = n(idx_CHj)
    fluxes(267) = n(idx_CH2)
    fluxes(268) = n(idx_CH2)
    fluxes(269) = n(idx_CH2j)
    fluxes(270) = n(idx_CH3)
    fluxes(271) = n(idx_CH3)
    fluxes(272) = n(idx_CH3)
    fluxes(273) = n(idx_CH4)
    fluxes(274) = n(idx_CH4)
    fluxes(275) = n(idx_CH4)
    fluxes(276) = n(idx_OH)
    fluxes(277) = n(idx_OH)
    fluxes(278) = n(idx_OHj)
    fluxes(279) = n(idx_H2O)
    fluxes(280) = n(idx_H2O)
    fluxes(281) = n(idx_O2)
    fluxes(282) = n(idx_O2)
    fluxes(283) = n(idx_H2)
    fluxes(284) = n(idx_H2)
    fluxes(285) = n(idx_H2)
    fluxes(286) = n(idx_H)
    fluxes(287) = n(idx_He)
    fluxes(288) = n(idx_CH3j)*n(idx_E)
    fluxes(289) = n(idx_H)*n(idx_H)

    ! multiply rate coefficient (vectorized)
    fluxes(:) = fluxes(:) * kall(:)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_FLUXES

    !!BEGIN_ODE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    dn(idx_H) = -fluxes(1) &
        -fluxes(2) &
        -fluxes(3) &
        -fluxes(4) &
        -fluxes(5) &
        -fluxes(6) &
        -fluxes(7) &
        +fluxes(8) &
        +fluxes(9) &
        +fluxes(10) &
        +fluxes(11) &
        +fluxes(12) &
        +fluxes(13) &
        +fluxes(17) &
        +fluxes(31) &
        -fluxes(35) &
        -fluxes(36) &
        +fluxes(37) &
        +fluxes(38) &
        +fluxes(40) &
        -fluxes(41) &
        -fluxes(42) &
        +fluxes(43) &
        +fluxes(46) &
        +fluxes(46) &
        +fluxes(47) &
        +fluxes(47) &
        +fluxes(47) &
        +fluxes(48) &
        +fluxes(49) &
        +fluxes(50) &
        +fluxes(51) &
        +fluxes(51) &
        +fluxes(53) &
        +fluxes(55) &
        +fluxes(55) &
        +fluxes(56) &
        +fluxes(57) &
        +fluxes(57) &
        +fluxes(58) &
        +fluxes(60) &
        +fluxes(61) &
        +fluxes(61) &
        +fluxes(63) &
        +fluxes(64) &
        +fluxes(65) &
        +fluxes(65) &
        +fluxes(67) &
        +fluxes(69) &
        +fluxes(71) &
        -fluxes(76) &
        -fluxes(77) &
        -fluxes(78) &
        -fluxes(83) &
        +fluxes(83) &
        +fluxes(83) &
        +fluxes(83) &
        -fluxes(84) &
        +fluxes(84) &
        +fluxes(84) &
        -fluxes(85) &
        +fluxes(85) &
        +fluxes(85) &
        -fluxes(86) &
        +fluxes(86) &
        +fluxes(86) &
        -fluxes(87) &
        +fluxes(87) &
        +fluxes(88) &
        +fluxes(88) &
        +fluxes(89) &
        +fluxes(89) &
        +fluxes(90) &
        +fluxes(91) &
        +fluxes(92) &
        -fluxes(95) &
        +fluxes(97) &
        +fluxes(99) &
        +fluxes(99) &
        -fluxes(101) &
        +fluxes(102) &
        +fluxes(104) &
        +fluxes(105) &
        -fluxes(107) &
        -fluxes(108) &
        -fluxes(114) &
        -fluxes(115) &
        +fluxes(116) &
        +fluxes(117) &
        +fluxes(118) &
        +fluxes(119) &
        +fluxes(120) &
        +fluxes(121) &
        +fluxes(122) &
        +fluxes(123) &
        +fluxes(124) &
        +fluxes(125) &
        +fluxes(126) &
        +fluxes(127) &
        +fluxes(128) &
        +fluxes(129) &
        +fluxes(130) &
        +fluxes(131) &
        +fluxes(141) &
        +fluxes(143) &
        +fluxes(145) &
        +fluxes(148) &
        +fluxes(149) &
        +fluxes(151) &
        +fluxes(155) &
        +fluxes(158) &
        +fluxes(164) &
        +fluxes(165) &
        +fluxes(178) &
        +fluxes(189) &
        +fluxes(190) &
        -fluxes(213) &
        +fluxes(214) &
        +fluxes(215) &
        +fluxes(216) &
        +fluxes(217) &
        +fluxes(218) &
        -fluxes(219) &
        +fluxes(258) &
        +fluxes(258) &
        +fluxes(260) &
        +fluxes(261) &
        +fluxes(265) &
        +fluxes(266) &
        +fluxes(268) &
        +fluxes(269) &
        +fluxes(271) &
        +fluxes(273) &
        +fluxes(275) &
        +fluxes(277) &
        +fluxes(279) &
        +fluxes(283) &
        +fluxes(284) &
        +fluxes(284) &
        -fluxes(286) &
        +fluxes(288) &
        -fluxes(289) &
        -fluxes(289)

    dn(idx_CH) = -fluxes(1) &
        +fluxes(8) &
        -fluxes(9) &
        +fluxes(15) &
        +fluxes(15) &
        +fluxes(16) &
        -fluxes(18) &
        -fluxes(19) &
        +fluxes(20) &
        +fluxes(21) &
        +fluxes(24) &
        +fluxes(50) &
        +fluxes(54) &
        +fluxes(55) &
        -fluxes(80) &
        -fluxes(84) &
        -fluxes(90) &
        -fluxes(94) &
        +fluxes(95) &
        -fluxes(97) &
        -fluxes(98) &
        +fluxes(107) &
        -fluxes(118) &
        -fluxes(133) &
        -fluxes(141) &
        -fluxes(165) &
        -fluxes(166) &
        -fluxes(168) &
        -fluxes(170) &
        -fluxes(172) &
        -fluxes(173) &
        -fluxes(174) &
        -fluxes(177) &
        +fluxes(183) &
        -fluxes(214) &
        -fluxes(220) &
        -fluxes(227) &
        -fluxes(231) &
        -fluxes(235) &
        -fluxes(236) &
        -fluxes(237) &
        -fluxes(238) &
        -fluxes(239) &
        -fluxes(264) &
        -fluxes(265) &
        +fluxes(268) &
        +fluxes(272) &
        +fluxes(275)

    dn(idx_C) = +fluxes(1) &
        +fluxes(6) &
        -fluxes(8) &
        -fluxes(15) &
        -fluxes(16) &
        -fluxes(17) &
        +fluxes(18) &
        +fluxes(39) &
        -fluxes(44) &
        +fluxes(49) &
        +fluxes(51) &
        +fluxes(52) &
        +fluxes(68) &
        +fluxes(73) &
        +fluxes(84) &
        +fluxes(90) &
        -fluxes(96) &
        -fluxes(107) &
        -fluxes(109) &
        -fluxes(110) &
        -fluxes(117) &
        -fluxes(132) &
        -fluxes(154) &
        -fluxes(156) &
        -fluxes(157) &
        -fluxes(159) &
        -fluxes(160) &
        -fluxes(163) &
        +fluxes(169) &
        +fluxes(173) &
        +fluxes(231) &
        +fluxes(232) &
        -fluxes(233) &
        -fluxes(234) &
        +fluxes(259) &
        -fluxes(263) &
        +fluxes(265) &
        +fluxes(288)

    dn(idx_H2) = +fluxes(1) &
        +fluxes(2) &
        +fluxes(3) &
        +fluxes(4) &
        +fluxes(5) &
        -fluxes(8) &
        -fluxes(9) &
        -fluxes(10) &
        -fluxes(11) &
        -fluxes(12) &
        -fluxes(13) &
        -fluxes(14) &
        +fluxes(35) &
        +fluxes(36) &
        -fluxes(37) &
        -fluxes(38) &
        -fluxes(40) &
        +fluxes(48) &
        +fluxes(52) &
        +fluxes(54) &
        +fluxes(59) &
        +fluxes(62) &
        +fluxes(64) &
        +fluxes(66) &
        -fluxes(79) &
        -fluxes(80) &
        -fluxes(81) &
        -fluxes(83) &
        -fluxes(88) &
        -fluxes(89) &
        -fluxes(89) &
        +fluxes(89) &
        -fluxes(90) &
        +fluxes(90) &
        -fluxes(91) &
        +fluxes(91) &
        -fluxes(92) &
        +fluxes(92) &
        -fluxes(93) &
        +fluxes(93) &
        +fluxes(95) &
        +fluxes(100) &
        +fluxes(101) &
        -fluxes(102) &
        -fluxes(106) &
        -fluxes(109) &
        +fluxes(112) &
        +fluxes(113) &
        +fluxes(114) &
        +fluxes(115) &
        -fluxes(116) &
        -fluxes(120) &
        -fluxes(122) &
        -fluxes(124) &
        +fluxes(125) &
        -fluxes(127) &
        -fluxes(129) &
        -fluxes(131) &
        +fluxes(132) &
        +fluxes(133) &
        +fluxes(134) &
        +fluxes(135) &
        +fluxes(136) &
        +fluxes(137) &
        +fluxes(138) &
        +fluxes(139) &
        +fluxes(140) &
        +fluxes(142) &
        +fluxes(144) &
        +fluxes(145) &
        +fluxes(146) &
        +fluxes(159) &
        +fluxes(167) &
        +fluxes(171) &
        +fluxes(186) &
        +fluxes(192) &
        +fluxes(213) &
        +fluxes(220) &
        +fluxes(221) &
        +fluxes(222) &
        +fluxes(223) &
        +fluxes(224) &
        +fluxes(225) &
        +fluxes(226) &
        -fluxes(258) &
        +fluxes(262) &
        +fluxes(272) &
        +fluxes(274) &
        +fluxes(275) &
        -fluxes(283) &
        -fluxes(284) &
        -fluxes(285) &
        +fluxes(288) &
        +fluxes(289)

    dn(idx_CH3) = -fluxes(2) &
        +fluxes(3) &
        +fluxes(10) &
        -fluxes(11) &
        +fluxes(19) &
        +fluxes(20) &
        +fluxes(22) &
        +fluxes(22) &
        +fluxes(23) &
        -fluxes(26) &
        -fluxes(26) &
        -fluxes(27) &
        -fluxes(28) &
        -fluxes(29) &
        +fluxes(30) &
        +fluxes(33) &
        +fluxes(56) &
        +fluxes(59) &
        +fluxes(74) &
        +fluxes(80) &
        -fluxes(135) &
        -fluxes(144) &
        +fluxes(147) &
        +fluxes(193) &
        +fluxes(196) &
        +fluxes(197) &
        +fluxes(198) &
        +fluxes(199) &
        -fluxes(216) &
        -fluxes(270) &
        -fluxes(271) &
        -fluxes(272) &
        +fluxes(273)

    dn(idx_CH2) = +fluxes(2) &
        +fluxes(9) &
        -fluxes(10) &
        -fluxes(15) &
        +fluxes(19) &
        -fluxes(20) &
        -fluxes(20) &
        -fluxes(21) &
        -fluxes(22) &
        -fluxes(23) &
        -fluxes(24) &
        -fluxes(25) &
        +fluxes(26) &
        +fluxes(28) &
        +fluxes(53) &
        +fluxes(57) &
        -fluxes(95) &
        -fluxes(99) &
        -fluxes(100) &
        +fluxes(109) &
        -fluxes(112) &
        -fluxes(119) &
        -fluxes(134) &
        -fluxes(142) &
        -fluxes(143) &
        -fluxes(179) &
        -fluxes(180) &
        -fluxes(181) &
        -fluxes(182) &
        -fluxes(183) &
        -fluxes(184) &
        +fluxes(191) &
        +fluxes(195) &
        -fluxes(215) &
        -fluxes(221) &
        -fluxes(232) &
        -fluxes(240) &
        -fluxes(241) &
        -fluxes(242) &
        -fluxes(243) &
        -fluxes(244) &
        -fluxes(267) &
        -fluxes(268) &
        +fluxes(271) &
        +fluxes(274)

    dn(idx_CH4) = -fluxes(3) &
        +fluxes(11) &
        -fluxes(19) &
        -fluxes(22) &
        +fluxes(26) &
        +fluxes(27) &
        +fluxes(29) &
        -fluxes(30) &
        -fluxes(33) &
        +fluxes(60) &
        -fluxes(113) &
        -fluxes(123) &
        -fluxes(125) &
        -fluxes(137) &
        -fluxes(145) &
        -fluxes(146) &
        -fluxes(147) &
        -fluxes(148) &
        +fluxes(156) &
        +fluxes(168) &
        +fluxes(180) &
        -fluxes(187) &
        -fluxes(193) &
        -fluxes(194) &
        -fluxes(195) &
        -fluxes(197) &
        -fluxes(199) &
        +fluxes(201) &
        +fluxes(207) &
        +fluxes(208) &
        -fluxes(217) &
        -fluxes(222) &
        -fluxes(228) &
        -fluxes(245) &
        -fluxes(250) &
        +fluxes(251) &
        -fluxes(273) &
        -fluxes(274) &
        -fluxes(275)

    dn(idx_OH) = -fluxes(4) &
        +fluxes(5) &
        +fluxes(6) &
        +fluxes(7) &
        +fluxes(12) &
        -fluxes(13) &
        +fluxes(14) &
        +fluxes(14) &
        -fluxes(16) &
        -fluxes(17) &
        +fluxes(18) &
        +fluxes(21) &
        -fluxes(23) &
        -fluxes(24) &
        -fluxes(27) &
        -fluxes(28) &
        +fluxes(29) &
        +fluxes(30) &
        -fluxes(31) &
        +fluxes(32) &
        +fluxes(32) &
        -fluxes(33) &
        -fluxes(34) &
        -fluxes(34) &
        +fluxes(63) &
        +fluxes(65) &
        +fluxes(66) &
        +fluxes(77) &
        -fluxes(78) &
        -fluxes(85) &
        +fluxes(86) &
        -fluxes(91) &
        +fluxes(92) &
        +fluxes(98) &
        -fluxes(103) &
        -fluxes(126) &
        -fluxes(138) &
        -fluxes(149) &
        +fluxes(150) &
        -fluxes(155) &
        +fluxes(157) &
        -fluxes(167) &
        +fluxes(170) &
        +fluxes(175) &
        +fluxes(181) &
        +fluxes(185) &
        +fluxes(187) &
        +fluxes(188) &
        -fluxes(189) &
        -fluxes(200) &
        -fluxes(201) &
        -fluxes(203) &
        -fluxes(205) &
        -fluxes(206) &
        +fluxes(209) &
        +fluxes(210) &
        +fluxes(211) &
        -fluxes(218) &
        -fluxes(223) &
        +fluxes(236) &
        +fluxes(241) &
        -fluxes(246) &
        +fluxes(252) &
        -fluxes(253) &
        +fluxes(254) &
        -fluxes(276) &
        -fluxes(277) &
        +fluxes(279)

    dn(idx_O) = +fluxes(4) &
        +fluxes(7) &
        -fluxes(12) &
        +fluxes(16) &
        -fluxes(18) &
        -fluxes(21) &
        +fluxes(23) &
        +fluxes(27) &
        -fluxes(30) &
        -fluxes(31) &
        -fluxes(32) &
        +fluxes(34) &
        +fluxes(42) &
        -fluxes(43) &
        +fluxes(45) &
        +fluxes(58) &
        +fluxes(61) &
        +fluxes(62) &
        +fluxes(64) &
        +fluxes(68) &
        +fluxes(70) &
        +fluxes(70) &
        +fluxes(75) &
        -fluxes(77) &
        -fluxes(82) &
        -fluxes(82) &
        +fluxes(85) &
        +fluxes(87) &
        +fluxes(87) &
        +fluxes(91) &
        +fluxes(93) &
        +fluxes(93) &
        -fluxes(94) &
        +fluxes(96) &
        -fluxes(97) &
        -fluxes(99) &
        -fluxes(100) &
        -fluxes(110) &
        -fluxes(111) &
        -fluxes(121) &
        -fluxes(136) &
        +fluxes(152) &
        +fluxes(153) &
        +fluxes(154) &
        +fluxes(161) &
        +fluxes(163) &
        -fluxes(164) &
        +fluxes(166) &
        +fluxes(176) &
        +fluxes(177) &
        -fluxes(178) &
        +fluxes(179) &
        -fluxes(186) &
        -fluxes(188) &
        -fluxes(190) &
        -fluxes(191) &
        -fluxes(192) &
        +fluxes(194) &
        +fluxes(200) &
        +fluxes(202) &
        +fluxes(203) &
        +fluxes(204) &
        +fluxes(205) &
        +fluxes(235) &
        +fluxes(240) &
        +fluxes(245) &
        +fluxes(246) &
        +fluxes(247) &
        -fluxes(248) &
        +fluxes(249) &
        +fluxes(259) &
        +fluxes(277) &
        +fluxes(278) &
        +fluxes(281) &
        +fluxes(281)

    dn(idx_H2O) = -fluxes(5) &
        +fluxes(13) &
        +fluxes(24) &
        +fluxes(25) &
        +fluxes(28) &
        -fluxes(29) &
        -fluxes(32) &
        +fluxes(33) &
        +fluxes(34) &
        +fluxes(67) &
        +fluxes(78) &
        -fluxes(86) &
        -fluxes(92) &
        -fluxes(104) &
        -fluxes(128) &
        -fluxes(139) &
        -fluxes(150) &
        -fluxes(151) &
        -fluxes(158) &
        -fluxes(169) &
        -fluxes(171) &
        +fluxes(172) &
        +fluxes(182) &
        -fluxes(196) &
        -fluxes(202) &
        -fluxes(207) &
        -fluxes(209) &
        -fluxes(211) &
        -fluxes(212) &
        -fluxes(224) &
        -fluxes(229) &
        +fluxes(237) &
        +fluxes(242) &
        -fluxes(247) &
        -fluxes(252) &
        -fluxes(255) &
        +fluxes(256) &
        -fluxes(279) &
        -fluxes(280)

    dn(idx_CO) = -fluxes(6) &
        +fluxes(17) &
        +fluxes(25) &
        -fluxes(39) &
        -fluxes(45) &
        +fluxes(69) &
        +fluxes(96) &
        +fluxes(97) &
        +fluxes(98) &
        +fluxes(99) &
        +fluxes(100) &
        +fluxes(103) &
        +fluxes(110) &
        -fluxes(130) &
        -fluxes(140) &
        -fluxes(152) &
        +fluxes(160) &
        +fluxes(162) &
        +fluxes(174) &
        +fluxes(184) &
        -fluxes(198) &
        -fluxes(204) &
        +fluxes(206) &
        -fluxes(208) &
        -fluxes(210) &
        +fluxes(212) &
        +fluxes(219) &
        -fluxes(225) &
        +fluxes(233) &
        +fluxes(238) &
        +fluxes(243) &
        +fluxes(248) &
        +fluxes(250) &
        +fluxes(253) &
        +fluxes(255) &
        +fluxes(257) &
        -fluxes(259)

    dn(idx_O2) = -fluxes(7) &
        -fluxes(14) &
        -fluxes(25) &
        +fluxes(31) &
        +fluxes(82) &
        -fluxes(87) &
        -fluxes(93) &
        -fluxes(96) &
        -fluxes(98) &
        -fluxes(105) &
        -fluxes(153) &
        -fluxes(161) &
        -fluxes(162) &
        -fluxes(175) &
        -fluxes(176) &
        -fluxes(185) &
        -fluxes(226) &
        -fluxes(230) &
        +fluxes(234) &
        +fluxes(239) &
        +fluxes(244) &
        -fluxes(249) &
        -fluxes(251) &
        -fluxes(254) &
        -fluxes(256) &
        -fluxes(257) &
        -fluxes(281) &
        -fluxes(282)

    dn(idx_CH2j) = -fluxes(35) &
        +fluxes(36) &
        -fluxes(50) &
        -fluxes(51) &
        -fluxes(52) &
        +fluxes(79) &
        +fluxes(102) &
        +fluxes(118) &
        -fluxes(120) &
        +fluxes(133) &
        +fluxes(146) &
        +fluxes(166) &
        +fluxes(168) &
        +fluxes(170) &
        +fluxes(172) &
        +fluxes(174) &
        -fluxes(178) &
        -fluxes(185) &
        +fluxes(215) &
        +fluxes(221) &
        +fluxes(232) &
        +fluxes(240) &
        +fluxes(241) &
        +fluxes(242) &
        +fluxes(243) &
        +fluxes(244) &
        +fluxes(267) &
        -fluxes(269)

    dn(idx_CHj) = +fluxes(35) &
        +fluxes(38) &
        -fluxes(49) &
        -fluxes(101) &
        -fluxes(102) &
        +fluxes(108) &
        +fluxes(112) &
        +fluxes(117) &
        +fluxes(132) &
        +fluxes(143) &
        +fluxes(144) &
        +fluxes(145) &
        +fluxes(154) &
        +fluxes(156) &
        +fluxes(157) &
        +fluxes(160) &
        -fluxes(164) &
        -fluxes(167) &
        -fluxes(169) &
        -fluxes(171) &
        -fluxes(175) &
        -fluxes(176) &
        +fluxes(214) &
        +fluxes(220) &
        +fluxes(227) &
        +fluxes(231) &
        +fluxes(235) &
        +fluxes(236) &
        +fluxes(237) &
        +fluxes(238) &
        +fluxes(239) &
        +fluxes(264) &
        -fluxes(266) &
        +fluxes(269)

    dn(idx_CH3j) = -fluxes(36) &
        -fluxes(53) &
        -fluxes(54) &
        -fluxes(55) &
        -fluxes(74) &
        -fluxes(81) &
        +fluxes(113) &
        +fluxes(114) &
        +fluxes(119) &
        +fluxes(120) &
        +fluxes(125) &
        +fluxes(134) &
        +fluxes(148) &
        +fluxes(179) &
        +fluxes(180) &
        +fluxes(181) &
        +fluxes(182) &
        +fluxes(184) &
        -fluxes(186) &
        +fluxes(187) &
        +fluxes(188) &
        +fluxes(216) &
        +fluxes(270) &
        -fluxes(288)

    dn(idx_Hej) = -fluxes(37) &
        -fluxes(39) &
        -fluxes(41) &
        -fluxes(44) &
        -fluxes(72) &
        -fluxes(106) &
        -fluxes(141) &
        -fluxes(142) &
        -fluxes(143) &
        -fluxes(144) &
        -fluxes(145) &
        -fluxes(146) &
        -fluxes(147) &
        -fluxes(148) &
        -fluxes(149) &
        -fluxes(150) &
        -fluxes(151) &
        -fluxes(152) &
        -fluxes(153) &
        -fluxes(227) &
        -fluxes(228) &
        -fluxes(229) &
        -fluxes(230) &
        +fluxes(287)

    dn(idx_He) = +fluxes(37) &
        +fluxes(39) &
        +fluxes(41) &
        +fluxes(44) &
        +fluxes(72) &
        +fluxes(106) &
        +fluxes(141) &
        +fluxes(142) &
        +fluxes(143) &
        +fluxes(144) &
        +fluxes(145) &
        +fluxes(146) &
        +fluxes(147) &
        +fluxes(148) &
        +fluxes(149) &
        +fluxes(150) &
        +fluxes(151) &
        +fluxes(152) &
        +fluxes(153) &
        +fluxes(227) &
        +fluxes(228) &
        +fluxes(229) &
        +fluxes(230) &
        -fluxes(287)

    dn(idx_Hj) = +fluxes(37) &
        -fluxes(40) &
        +fluxes(41) &
        +fluxes(42) &
        -fluxes(43) &
        -fluxes(71) &
        -fluxes(76) &
        +fluxes(103) &
        -fluxes(104) &
        -fluxes(105) &
        -fluxes(112) &
        -fluxes(113) &
        +fluxes(147) &
        +fluxes(150) &
        +fluxes(213) &
        -fluxes(214) &
        -fluxes(215) &
        -fluxes(216) &
        -fluxes(217) &
        -fluxes(218) &
        +fluxes(219) &
        +fluxes(260) &
        +fluxes(262) &
        +fluxes(278) &
        +fluxes(283) &
        +fluxes(286)

    dn(idx_Cj) = -fluxes(38) &
        +fluxes(44) &
        -fluxes(73) &
        -fluxes(79) &
        +fluxes(101) &
        -fluxes(103) &
        -fluxes(108) &
        -fluxes(111) &
        +fluxes(141) &
        +fluxes(142) &
        +fluxes(152) &
        -fluxes(155) &
        -fluxes(158) &
        -fluxes(161) &
        -fluxes(162) &
        -fluxes(231) &
        -fluxes(232) &
        +fluxes(233) &
        +fluxes(234) &
        +fluxes(263) &
        +fluxes(266)

    dn(idx_Oj) = +fluxes(39) &
        -fluxes(42) &
        +fluxes(43) &
        -fluxes(45) &
        -fluxes(75) &
        -fluxes(122) &
        +fluxes(149) &
        +fluxes(153) &
        +fluxes(162) &
        -fluxes(165) &
        -fluxes(187) &
        -fluxes(189) &
        -fluxes(235) &
        -fluxes(240) &
        -fluxes(245) &
        -fluxes(246) &
        -fluxes(247) &
        +fluxes(248) &
        -fluxes(249)

    dn(idx_H2j) = +fluxes(40) &
        -fluxes(46) &
        +fluxes(76) &
        +fluxes(106) &
        -fluxes(116) &
        -fluxes(117) &
        -fluxes(118) &
        -fluxes(119) &
        -fluxes(121) &
        -fluxes(123) &
        -fluxes(125) &
        -fluxes(126) &
        -fluxes(128) &
        -fluxes(130) &
        -fluxes(213) &
        -fluxes(220) &
        -fluxes(221) &
        -fluxes(222) &
        -fluxes(223) &
        -fluxes(224) &
        -fluxes(225) &
        -fluxes(226) &
        -fluxes(260) &
        +fluxes(261) &
        +fluxes(285)

    dn(idx_COj) = +fluxes(45) &
        -fluxes(68) &
        +fluxes(111) &
        -fluxes(131) &
        +fluxes(155) &
        +fluxes(161) &
        +fluxes(163) &
        +fluxes(164) &
        +fluxes(165) &
        +fluxes(167) &
        -fluxes(173) &
        +fluxes(175) &
        -fluxes(183) &
        -fluxes(199) &
        -fluxes(205) &
        -fluxes(211) &
        -fluxes(219) &
        +fluxes(225) &
        -fluxes(233) &
        -fluxes(238) &
        -fluxes(243) &
        -fluxes(248) &
        -fluxes(250) &
        -fluxes(253) &
        -fluxes(255) &
        -fluxes(257)

    dn(idx_E) = -fluxes(46) &
        -fluxes(47) &
        -fluxes(48) &
        -fluxes(49) &
        -fluxes(50) &
        -fluxes(51) &
        -fluxes(52) &
        -fluxes(53) &
        -fluxes(54) &
        -fluxes(55) &
        -fluxes(56) &
        -fluxes(57) &
        -fluxes(58) &
        -fluxes(59) &
        -fluxes(60) &
        -fluxes(61) &
        -fluxes(62) &
        -fluxes(63) &
        -fluxes(64) &
        -fluxes(65) &
        -fluxes(66) &
        -fluxes(67) &
        -fluxes(68) &
        -fluxes(69) &
        -fluxes(70) &
        -fluxes(71) &
        -fluxes(72) &
        -fluxes(73) &
        -fluxes(74) &
        -fluxes(75) &
        -fluxes(88) &
        +fluxes(88) &
        +fluxes(94) &
        +fluxes(263) &
        +fluxes(264) &
        +fluxes(267) &
        +fluxes(270) &
        +fluxes(276) &
        +fluxes(280) &
        +fluxes(282) &
        +fluxes(283) &
        +fluxes(285) &
        +fluxes(286) &
        +fluxes(287) &
        -fluxes(288)

    dn(idx_H3j) = -fluxes(47) &
        -fluxes(48) &
        +fluxes(116) &
        -fluxes(132) &
        -fluxes(133) &
        -fluxes(134) &
        -fluxes(135) &
        -fluxes(136) &
        -fluxes(137) &
        -fluxes(138) &
        -fluxes(139) &
        -fluxes(140) &
        -fluxes(261) &
        -fluxes(262)

    dn(idx_CH4j) = -fluxes(56) &
        -fluxes(57) &
        -fluxes(114) &
        +fluxes(115) &
        -fluxes(124) &
        +fluxes(135) &
        -fluxes(188) &
        -fluxes(193) &
        -fluxes(196) &
        -fluxes(198) &
        +fluxes(217) &
        +fluxes(222) &
        +fluxes(228) &
        +fluxes(245) &
        +fluxes(250) &
        -fluxes(251)

    dn(idx_OHj) = -fluxes(58) &
        +fluxes(121) &
        +fluxes(122) &
        -fluxes(127) &
        +fluxes(136) &
        +fluxes(151) &
        -fluxes(154) &
        -fluxes(166) &
        -fluxes(179) &
        -fluxes(190) &
        -fluxes(194) &
        -fluxes(195) &
        -fluxes(200) &
        -fluxes(202) &
        -fluxes(204) &
        +fluxes(218) &
        +fluxes(223) &
        -fluxes(236) &
        -fluxes(241) &
        +fluxes(246) &
        -fluxes(252) &
        +fluxes(253) &
        -fluxes(254) &
        +fluxes(276) &
        -fluxes(278)

    dn(idx_CH5j) = -fluxes(59) &
        -fluxes(60) &
        +fluxes(81) &
        -fluxes(115) &
        +fluxes(123) &
        +fluxes(124) &
        +fluxes(137) &
        -fluxes(156) &
        -fluxes(168) &
        -fluxes(180) &
        -fluxes(191) &
        +fluxes(193) &
        +fluxes(194) &
        -fluxes(201) &
        -fluxes(207) &
        -fluxes(208)

    dn(idx_H2Oj) = -fluxes(61) &
        -fluxes(62) &
        -fluxes(63) &
        +fluxes(104) &
        +fluxes(126) &
        +fluxes(127) &
        -fluxes(129) &
        +fluxes(138) &
        -fluxes(157) &
        -fluxes(170) &
        -fluxes(181) &
        -fluxes(192) &
        -fluxes(197) &
        +fluxes(200) &
        +fluxes(201) &
        -fluxes(203) &
        +fluxes(206) &
        -fluxes(209) &
        -fluxes(210) &
        +fluxes(224) &
        +fluxes(229) &
        -fluxes(237) &
        -fluxes(242) &
        +fluxes(247) &
        +fluxes(252) &
        +fluxes(255) &
        -fluxes(256) &
        +fluxes(280)

    dn(idx_H3Oj) = -fluxes(64) &
        -fluxes(65) &
        -fluxes(66) &
        -fluxes(67) &
        +fluxes(128) &
        +fluxes(129) &
        +fluxes(139) &
        -fluxes(159) &
        +fluxes(169) &
        -fluxes(172) &
        -fluxes(182) &
        +fluxes(191) &
        +fluxes(195) &
        +fluxes(196) &
        +fluxes(197) &
        +fluxes(202) &
        +fluxes(203) &
        +fluxes(207) &
        +fluxes(209) &
        +fluxes(212)

    dn(idx_HCOj) = -fluxes(69) &
        +fluxes(94) &
        +fluxes(130) &
        +fluxes(131) &
        +fluxes(140) &
        +fluxes(158) &
        +fluxes(159) &
        -fluxes(160) &
        +fluxes(171) &
        +fluxes(173) &
        -fluxes(174) &
        +fluxes(176) &
        +fluxes(177) &
        +fluxes(178) &
        +fluxes(183) &
        -fluxes(184) &
        +fluxes(185) &
        +fluxes(186) &
        +fluxes(198) &
        +fluxes(199) &
        +fluxes(204) &
        +fluxes(205) &
        -fluxes(206) &
        +fluxes(208) &
        +fluxes(210) &
        +fluxes(211) &
        -fluxes(212)

    dn(idx_O2j) = -fluxes(70) &
        +fluxes(105) &
        -fluxes(163) &
        -fluxes(177) &
        +fluxes(189) &
        +fluxes(190) &
        +fluxes(192) &
        +fluxes(226) &
        +fluxes(230) &
        -fluxes(234) &
        -fluxes(239) &
        -fluxes(244) &
        +fluxes(249) &
        +fluxes(251) &
        +fluxes(254) &
        +fluxes(256) &
        +fluxes(257) &
        +fluxes(282)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ODE

    ! chemical species differntials are copied back to the solver variables
    dy(1:nmols) = dn(:)

  end subroutine fex

  ! ***************************
  ! Jacobian, pd(i,j)=df(i)/dx(j), see DLSODES documentation
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use prizmo_commons
    implicit none
    integer::neq, j, ian, jan
    real*8::tt, n(neq), pdj(neq)

    !!BEGIN_JACOBIAN

    !!END_JACOBIAN

  end subroutine jes

end module prizmo_ode

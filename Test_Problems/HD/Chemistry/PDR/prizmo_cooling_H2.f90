module prizmo_cooling_H2
contains

  function cooling_H2(n, Tgas)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas
    real*8::T3, logt3, logt, cool, w14, w24, cooling_H2
    real*8::logt32, logt33, logt34, logt35, logt36, logt37, logt38
    real*8::HDLR, HDLV, HDL, LDL, dump14, temp
    real*8::fH2H, fH2Hp, fH2H2, fH2e, fH2He

    temp = Tgas
    cooling_H2 = 0d0

    T3 = temp * 1d-3
    logt3 = log10(T3)
    logt = log10(temp)
    cool = 0d0

    logt32 = logt3 * logt3
    logt33 = logt32 * logt3
    logt34 = logt33 * logt3
    logt35 = logt34 * logt3
    logt36 = logt35 * logt3
    logt37 = logt36 * logt3
    logt38 = logt37 * logt3

    w14 = wCool(logt, 1d0, 4d0)
    w24 = wCool(logt, 2d0, 4d0)

    !//H2-H
    if(temp <= 1d2) then
       fH2H = 1d1**(-16.818342D0 +3.7383713D1*logt3 &
            +5.8145166D1*logt32 +4.8656103D1*logt33 &
            +2.0159831D1*logt34 +3.8479610D0*logt35) * n(idx_H)
    elseif(temp > 1d2 .and. temp <= 1d3) then
       fH2H = 1d1**(-2.4311209D1 +3.5692468D0*logt3 &
            -1.1332860D1*logt32 -2.7850082D1*logt33 &
            -2.1328264D1*logt34 -4.2519023D0*logt35) * n(idx_H)
    elseif(temp > 1d3 .and. temp <= 6d3) then
       fH2H = 1d1**(-2.4311209D1 +4.6450521D0*logt3 &
            -3.7209846D0*logt32 +5.9369081D0*logt33 &
            -5.5108049D0*logt34 +1.5538288D0*logt35) * n(idx_H)
    else
       fH2H = 1.862314467912518E-022 * wCool(logt, 1d0, log10(6d3)) * n(idx_H)
    end if
    cool = cool + fH2H

    !//H2-Hp
    if(temp > 1d1 .and. temp <= 1d4) then
       fH2Hp = 1d1**(-2.2089523d1 +1.5714711d0*logt3 &
            +0.015391166d0*logt32 -0.23619985d0*logt33 &
            -0.51002221d0*logt34 +0.32168730d0*logt35) * n(idx_Hj)
    else
       fH2Hp = 1.182509139382060E-021 * n(idx_Hj) * w14
    endif
    cool = cool + fH2Hp

    !//H2-H2
    fH2H2 = w24 * 1d1**(-2.3962112D1 +2.09433740D0*logt3 &
         -.77151436D0*logt32 +.43693353D0*logt33 &
         -.14913216D0*logt34 -.033638326D0*logt35) * n(idx_H2)
    cool = cool + fH2H2

    !//H2-e
    fH2e = 0d0
    if(temp <= 5d2) then
       fH2e = 1d1**(min(-2.1928796d1 + 1.6815730d1*logt3 &
            +9.6743155d1*logt32 +3.4319180d2*logt33 &
            +7.3471651d2*logt34 +9.8367576d2*logt35 &
            +8.0181247d2*logt36 +3.6414446d2*logt37 &
            +7.0609154d1*logt38, 3d1)) * n(idx_e)
    elseif(temp > 5d2) then
       fH2e = 1d1**(-2.2921189D1 +1.6802758D0*logt3 &
            +.93310622D0*logt32 +4.0406627d0*logt33 &
            -4.7274036d0*logt34 -8.8077017d0*logt35 &
            +8.9167183*logt36 + 6.4380698*logt37 &
            -6.3701156*logt38) * n(idx_e)
    end if
    cool = cool + fH2e * w24

    !//H2-He
    if(temp > 1d1 .and. temp <= 1d4) then
       fH2He = 1d1**(-2.3689237d1 +2.1892372d0*logt3 &
            -.81520438d0*logt32 +.29036281d0*logt33 -.16596184d0*logt34 &
            +.19191375d0*logt35) * n(idx_He)
    else
       fH2He = 1.002560385050777E-022 * n(idx_He) * w14
    endif
    cool = cool + fH2He


    !this to avoid negative, overflow and useless calculations below
    if(cool <= 0d0) then
       cooling_H2 = 0d0
       return
    end if

    !high density limit from HM79, GP98 below Tgas = 2d3
    !UPDATED USING GLOVER 2015 for high temperature corrections, MNRAS
    !IN THE HIGH DENSITY REGIME LAMBDA_H2 = LAMBDA_H2(LTE) = HDL
    !the following mix of functions ensures the right behaviour
    ! at low (T<10 K) and high temperatures (T>2000 K) by
    ! using both the original Hollenbach and the new Glover data
    ! merged in a smooth way.
    if(temp < 2d3) then
       HDLR = ((9.5e-22 * t3**3.76) / (1. + 0.12 * t3**2.1)*exp(-(0.13 / t3)**3) &
            + 3.e-24 * exp(-0.51 / t3)) !erg/s
       HDLV = (6.7e-19 * exp(-5.86 / t3) + 1.6e-18 * exp(-11.7 / t3)) !erg/s
       HDL  = HDLR + HDLV !erg/s
    elseif(temp >= 2d3 .and. temp <= 1d4) then
       HDL = 1d1**(-2.0584225d1 + 5.0194035*logt3 &
            -1.5738805*logt32 -4.7155769*logt33 &
            +2.4714161*logt34 +5.4710750*logt35 &
            -3.9467356*logt36 -2.2148338*logt37 &
            +1.8161874*logt38)
    else
       dump14 = 1d0 / (1d0 + exp(min((temp-3d4) * 2d-4, 3d2)))
       HDL = 5.531333679406485e-19 * dump14
    endif

    LDL = cool !erg/s
    if (HDL == 0d0) then
       cooling_H2 = 0d0
    else
       cooling_H2 = max(n(idx_H2), 0d0) / (1d0 / HDL + 1d0 / LDL) !erg/cm3/s
    endif

  end function cooling_H2

  !*****************
  !sigmoid function with x0 shift and s steepness
  function sigmoid(x, x0, s)
    implicit none
    real*8,intent(in)::x, x0, s
    real*8::sigmoid

    sigmoid = 1d1 / (1d1 + exp(-s * (x - x0)))

  end function sigmoid

  !*******************
  !window function for H2 cooling to smooth limits
  function wCool(logTgas, logTmin, logTmax)
    implicit none
    real*8,intent(in)::logTgas, logTmin, logTmax
    real*8::wCool, x

    x = (logTgas - logTmin) / (logTmax - logTmin)
    wCool = 1d1**(2d2 * (sigmoid(x, -2d-1, 5d1) * sigmoid(-x, -1.2d0, 5d1) - 1d0))
    if(wCool < 1d-199) wCool = 0d0
    if(wCool > 1d0) then
       print *, "ERROR: H2 cooling window is wCool > 1"
       stop
    end if

  end function wCool


end module prizmo_cooling_H2

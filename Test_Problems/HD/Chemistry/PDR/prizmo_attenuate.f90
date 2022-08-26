module prizmo_attenuate
  use prizmo_commons
contains

  ! *******************
  subroutine attenuate(x, Tgas, jflux, ds)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, ds
    real*8,intent(inout)::jflux(nphoto)
    real*8::tau(nphoto), rhogas

    tau = 0d0

    !! PREPROCESS_ATTENUATE

    tau = tau + x(idx_H) * photo_xsecs(:, 283)
    tau = tau + x(idx_He) * photo_xsecs(:, 284)
    tau = tau + x(idx_Hej) * photo_xsecs(:, 285)
    tau = tau + x(idx_O) * photo_xsecs(:, 286)
    tau = tau + x(idx_Oj) * photo_xsecs(:, 287)
    tau = tau + x(idx_Ojj) * photo_xsecs(:, 288)
    tau = tau + x(idx_Ojjj) * photo_xsecs(:, 289)
    tau = tau + x(idx_C) * photo_xsecs(:, 290)
    tau = tau + x(idx_Cj) * photo_xsecs(:, 291)
    tau = tau + x(idx_Cjj) * photo_xsecs(:, 292)
    tau = tau + x(idx_Cjjj) * photo_xsecs(:, 293)
    tau = tau + x(idx_H2j) * photo_xsecs(:, 296)
    tau = tau + x(idx_CH) * photo_xsecs(:, 297)
    tau = tau + x(idx_CH) * photo_xsecs(:, 298)
    tau = tau + x(idx_CHj) * photo_xsecs(:, 299)
    tau = tau + x(idx_CH2) * photo_xsecs(:, 300)
    tau = tau + x(idx_CH2j) * photo_xsecs(:, 301)
    tau = tau + x(idx_CH3) * photo_xsecs(:, 302)
    tau = tau + x(idx_CH3) * photo_xsecs(:, 303)
    tau = tau + x(idx_CH4) * photo_xsecs(:, 304)
    tau = tau + x(idx_OH) * photo_xsecs(:, 305)
    tau = tau + x(idx_OHj) * photo_xsecs(:, 306)
    tau = tau + x(idx_H2O) * photo_xsecs(:, 307)
    tau = tau + x(idx_H2O) * photo_xsecs(:, 308)
    tau = tau + x(idx_O2) * photo_xsecs(:, 309)
    tau = tau + x(idx_O2) * photo_xsecs(:, 310)

    !! PREPROCESS_END

    rhogas = get_rho(x)
    tau = tau + dust_kappa_opacity * rhogas * d2g

    jflux = jflux * exp(-tau * ds)

  end subroutine

end module prizmo_attenuate

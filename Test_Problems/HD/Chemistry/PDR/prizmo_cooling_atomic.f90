module prizmo_cooling_atomic
  use prizmo_commons
  use prizmo_fit
contains

  ! ***************
  function cooling_atomic(x, log_Tgas_in, Tgas) result(cool)
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas_in, Tgas
    real*8::cool !, coola2d(atomic_cooling_nvec2d), coola3d(atomic_cooling_nvec3d), coola4d(atomic_cooling_nvec4d)
    real*8::log_xnr, log_xpr, log_xer, log_tgas
    real*8::nmin, nmax, tmin, tmax

    nmin = -6d0
    nmax = 2d1 * 0.99999
    tmin = 0d0
    tmax = 6d0 * 0.99999

    log_xnr = min(max(log10(x(idx_H) + 1d-40), nmin), nmax)
    log_xpr = min(max(log10(x(idx_Hj) + 1d-40), nmin), nmax)
    log_xer = min(max(log10(x(idx_E) + 1d-40), nmin), nmax)
    log_tgas = min(max(log_Tgas_in, tmin), tmax)

    !log_xnr = min(max(log10(x(idx_H) + 1d-40), atomic_cooling_table%xmin), atomic_cooling_table%xmax*0.999999)
    !log_xpr = min(max(log10(x(idx_Hj) + 1d-40), atomic_cooling_table%ymin), atomic_cooling_table%ymax*0.999999)
    !log_xer = min(max(log10(x(idx_E) + 1d-40), atomic_cooling_table%zmin), atomic_cooling_table%zmax*0.999999)
    !log_tgas = min(log_Tgas_in, atomic_cooling_table%umax*0.99999)

    ! coola2d = cooling_atomic_array2d(log_xer, log_tgas)
    ! coola3d = cooling_atomic_array3d(log_xnr, log_xer, log_tgas)
    ! coola4d = cooling_atomic_array4d(log_xpr, log_xnr, log_xer, log_tgas)

    cool = 0d0

    ! Ly-alpha H cooling
    cool = cool + 7.39e-19 * x(idx_H) * x(idx_E) * exp(-118400. / tgas)

    !! PREPROCESS_ATOMIC_COOLING

    cool = cool + atomic_cooling_C(x, log_Tgas)
    cool = cool + atomic_cooling_O(x, log_Tgas)
    cool = cool + atomic_cooling_Cj(x, log_Tgas)
    cool = cool + atomic_cooling_Oj(x, log_Tgas)

    !! PREPROCESS_END

  end function cooling_atomic

  !! PREPROCESS_ATOMIC_COOLING_FUNCTIONS

  ! *****************
  function atomic_cooling_C(x, log_Tgas) result(cool)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas
    real*8::cool, b(3), A(3, 3), n(3), H2or, H2pa
    real*8::kfit_H(atomic_cooling_3lev_nvec), kfit_H2or(atomic_cooling_3lev_nvec), kfit_H2pa(atomic_cooling_3lev_nvec), kfit_Hj(atomic_cooling_3lev_nvec), kfit_e(atomic_cooling_3lev_nvec)
  
    H2or = x(idx_H2) * ortho_to_para / (ortho_to_para + 1d0)
    H2pa = x(idx_H2) / (ortho_to_para + 1d0)
  
    kfit_Hj = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_C_Hj, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_H = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_C_H, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_e = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_C_e, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_H2or = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_C_H2or, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_H2pa = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_C_H2pa, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    A(2, 1) = -kfit_Hj(1) * x(idx_Hj) - kfit_H(1) * x(idx_H) - kfit_e(1) * x(idx_e) - kfit_H2or(1) * H2or - kfit_H2pa(1) * H2pa
    A(2, 2) = 7.900000e-08 + kfit_Hj(2) * x(idx_Hj) + kfit_H(2) * x(idx_H) + kfit_e(2) * x(idx_e) + kfit_H2or(2) * H2or + kfit_H2pa(2) * H2pa
    A(2, 3) = 2.100000e-14 + kfit_Hj(3) * x(idx_Hj) + kfit_H(3) * x(idx_H) + kfit_e(3) * x(idx_e) + kfit_H2or(3) * H2or + kfit_H2pa(3) * H2pa
    A(3, 1) = kfit_Hj(4) * x(idx_Hj) + kfit_H(4) * x(idx_H) + kfit_e(4) * x(idx_e) + kfit_H2or(4) * H2or + kfit_H2pa(4) * H2pa
    A(3, 2) = kfit_Hj(5) * x(idx_Hj) + kfit_H(5) * x(idx_H) + kfit_e(5) * x(idx_e) + kfit_H2or(5) * H2or + kfit_H2pa(5) * H2pa
    A(3, 3) = -2.700000e-07 - 2.100000e-14 - kfit_Hj(6) * x(idx_Hj) - kfit_H(6) * x(idx_H) - kfit_e(6) * x(idx_e) - kfit_H2or(6) * H2or - kfit_H2pa(6) * H2pa
  
    b = (/1d0, 0d0, 0d0/)
    n = linear_solver_n3(A, b)
  
    cool = n(3) * (2.100000e-14 * 8.623590e-15 + 2.700000e-07 * 5.362476e-15) + n(2) * 7.900000e-08 * 3.261114e-15
    cool = cool * x(idx_C)
  
  end function atomic_cooling_C
  
  ! *****************
  function atomic_cooling_O(x, log_Tgas) result(cool)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas
    real*8::cool, b(3), A(3, 3), n(3), H2or, H2pa
    real*8::kfit_H(atomic_cooling_3lev_nvec), kfit_Hj(atomic_cooling_3lev_nvec), kfit_e(atomic_cooling_3lev_nvec)
  
    kfit_H = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_O_H, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_Hj = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_O_Hj, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    kfit_e = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_O_e, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    A(2, 1) = -kfit_H(1) * x(idx_H) - kfit_Hj(1) * x(idx_Hj) - kfit_e(1) * x(idx_e)
    A(2, 2) = 8.900000e-05 + kfit_H(2) * x(idx_H) + kfit_Hj(2) * x(idx_Hj) + kfit_e(2) * x(idx_e)
    A(2, 3) = 1.300000e-10 + kfit_H(3) * x(idx_H) + kfit_Hj(3) * x(idx_Hj) + kfit_e(3) * x(idx_e)
    A(3, 1) = kfit_H(4) * x(idx_H) + kfit_Hj(4) * x(idx_Hj) + kfit_e(4) * x(idx_e)
    A(3, 2) = kfit_H(5) * x(idx_H) + kfit_Hj(5) * x(idx_Hj) + kfit_e(5) * x(idx_e)
    A(3, 3) = -1.800000e-05 - 1.300000e-10 - kfit_H(6) * x(idx_H) - kfit_Hj(6) * x(idx_Hj) - kfit_e(6) * x(idx_e)
  
    b = (/1d0, 0d0, 0d0/)
    n = linear_solver_n3(A, b)
  
    cool = n(3) * (1.300000e-10 * 4.508815e-14 + 1.800000e-05 * 1.364932e-14) + n(2) * 8.900000e-05 * 3.143883e-14
    cool = cool * x(idx_O)
  
  end function atomic_cooling_O
  
  ! *****************
  function atomic_cooling_Cj(x, log_Tgas) result(cool)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas
    real*8::cool, b(2), A(2, 2), n(2)
    real*8::kfit_H(atomic_cooling_2lev_nvec), kfit_e(atomic_cooling_2lev_nvec)
  
    kfit_e = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_Cj_e, & 
      atomic_cooling_2lev_nvec, atomic_cooling_n1)
  
    kfit_H = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_Cj_H, & 
      atomic_cooling_2lev_nvec, atomic_cooling_n1)
  
    A(2, 1) = -kfit_e(1) * x(idx_e) - kfit_H(1) * x(idx_H)
    A(2, 2) = 2.400000e-06 + kfit_e(2) * x(idx_e) + kfit_H(2) * x(idx_H)
  
    b = (/1d0, 0d0/)
    n = linear_solver_n2(A, b)
  
    cool = n(2) * 2.400000e-06 * 1.259850e-14
    cool = cool * x(idx_Cj)
  
  end function atomic_cooling_Cj
  
  ! *****************
  function atomic_cooling_Oj(x, log_Tgas) result(cool)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    real*8,intent(in)::x(nspecies), log_Tgas
    real*8::cool, b(3), A(3, 3), n(3), H2or, H2pa
    real*8::kfit_e(atomic_cooling_3lev_nvec)
  
    kfit_e = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_Oj_e, & 
      atomic_cooling_3lev_nvec, atomic_cooling_n1)
  
    A(2, 1) = -kfit_e(1) * x(idx_e)
    A(2, 2) = 5.100000e-05 + kfit_e(2) * x(idx_e)
    A(2, 3) = 1.700000e-04 + kfit_e(3) * x(idx_e)
    A(3, 1) = kfit_e(4) * x(idx_e)
    A(3, 2) = kfit_e(5) * x(idx_e)
    A(3, 3) = -1.300000e-07 - 1.700000e-04 - kfit_e(6) * x(idx_e)
  
    b = (/1d0, 0d0, 0d0/)
    n = linear_solver_n3(A, b)
  
    cool = n(3) * (1.700000e-04 * 5.329782e-12 + 1.300000e-07 * 3.976295e-15) + n(2) * 5.100000e-05 * 5.325805e-12
    cool = cool * x(idx_Oj)
  
  end function atomic_cooling_Oj

  !! PREPROCESS_END

  ! ! **********************
  ! function cooling_atomic_array2d(log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec2d)
  !
  !   coola(:) = 1d1**interp_2dfit_vec(log_xer, log_Tgas, atomic_cooling_table_2d, &
  !         atomic_cooling_nvec2d, atomic_cooling2d_n1, atomic_cooling2d_n2)
  !
  ! end function cooling_atomic_array2d
  !
  ! ! **********************
  ! function cooling_atomic_array3d(log_xnr, log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xnr, log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec3d)
  !
  !   coola(:) = 1d1**interp_3dfit_vec(log_xnr, log_xer, log_Tgas, atomic_cooling_table_3d, &
  !         atomic_cooling_nvec3d, atomic_cooling3d_n1, atomic_cooling3d_n2, atomic_cooling3d_n3)
  !
  ! end function cooling_atomic_array3d
  !
  ! ! **********************
  ! function cooling_atomic_array4d(log_xpr, log_xnr, log_xer, log_tgas) result(coola)
  !   implicit none
  !   real*8,intent(in)::log_xpr, log_xnr, log_xer, log_tgas
  !   real*8::coola(atomic_cooling_nvec4d)
  !
  !   coola(:) = 1d1**interp_4dfit_vec(log_xnr, log_xpr, log_xer, log_Tgas, atomic_cooling_table_4d, &
  !         atomic_cooling_nvec4d, atomic_cooling4d_n1, atomic_cooling4d_n2, atomic_cooling4d_n3, atomic_cooling4d_n4)
  !
  ! end function cooling_atomic_array4d

end module prizmo_cooling_atomic

module prizmo_equilibrium
  real*8::ntot_H, ntot_C, ntot_O, ntot_He
contains

! ***************************
subroutine NR_equilibrium(x)
  use prizmo_commons
  use prizmo_ode
  use prizmo_utils
  implicit none
  integer,parameter::neq=nmols+1
  real*8,intent(inout)::x(neq)
  real*8::tt, df(nmols+1)
  integer::step, i

  ntot_H = get_Hnuclei(x(1:nmols))
  ntot_He = get_Henuclei(x(1:nmols))
  ntot_C = get_Cnuclei(x(1:nmols))
  ntot_O = get_Onuclei(x(1:nmols))

  tt = 0d0
  step = 0
  do
    call NR_step(x(:), 1d-10)

    do i=1,neq
      x(i) = max(x(i), 0d0)
    end do

    df = f(x(:))

    print *, step, maxval(abs(df(:))), minval(x), maxval(x)

    if(maxval(abs(df(:))) < 1d-15) then
      exit
    end if
    step = step + 1
  end do

end subroutine NR_equilibrium

! *************************
subroutine NR_step(x, h2)
  use prizmo_commons
  use prizmo_ode
  use prizmo_linear_solver
  implicit none
  integer,parameter::neq=nmols+1
  real*8,intent(inout)::x(neq)
  real*8,intent(in)::h2
  real*8::jac(neq, neq), df(neq)
  integer::i,j

  df(:) = -f(x(:))
  jac(:, :) = jex(x(:), h2)

  !do i=1,nmols
  !  do j=1,nmols
  !    write(22, *) i, j, jac(i, j)
  !  end do
  !  write(22, *)
  !end do

  !x(1:nmols) = matmul(jac, x(1:nmols))
  !
  do i=1,nmols+1
    print *, i, -df(i)
  end do

  df(:) = linear_solver_nn(jac(:, :), df(:), neq)

  x(:) = df(:) + x(:)

end subroutine NR_step

! *******************
function f(x)
  use prizmo_commons
  use prizmo_ode
  use prizmo_utils
  implicit none
  real*8,intent(in)::x(nmols+1)
  real*8::df(nmols+1), f(nmols+1), tt

  tt = 0d0

  call fex(nmols+1, tt, x(:), df(:))
  f(:) = df(:)

  f(idx_H3Oj) = get_Hnuclei(x(1:nmols)) - ntot_H
  f(idx_Hejj) = get_Henuclei(x(1:nmols)) - ntot_He
  f(idx_CH4j) = get_Cnuclei(x(1:nmols)) - ntot_C
  f(idx_Ojjj) = get_Onuclei(x(1:nmols)) - ntot_O
  f(idx_Cjjj) = get_electrons(x(1:nmols)) - x(idx_E)

end function

! ********************
function jex(x, h2in)
  use prizmo_commons
  use prizmo_ode
  implicit none
  real*8,intent(in)::x(nmols+1), h2in
  real*8::jex(nmols+1, nmols+1), tt, yorg, idh(nmols+1), h2(nmols+1)
  real*8::dxu(nmols+1), dxl(nmols+1), y(nmols+1), lb(nmols+1), ub(nmols+1)
  integer::i

  h2(:) = h2in
  h2(nmols+1) = 1d-5
  idh(:) = 1d0 / 2e0 / h2(:)
  tt = 0d0

  lb(:) = 0d0
  ub(:) = 1d20
  lb(:) = Tgas_min
  ub(:) = Tgas_max

  y(:) = x(:)

  do i=1,nmols+1
    y(i) = max(x(i) - h2(i), lb(i))
    dxl(:) = f(y(:))
    y(i) = min(x(i) + h2(i), ub(i))
    dxu(:) = f(y(:))
    jex(:, i) = (dxu(:) - dxl(:)) * idh(i)
    y(i) = x(i)
  end do

end function jex


end module prizmo_equilibrium

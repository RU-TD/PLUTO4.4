module prizmo_linear_solver
contains

  ! *******************************
  ! solve Ax=B linear system
  function linear_solver(a, b) result(c)
    implicit none
    real*8,intent(in)::a(:, :), b(:)
    integer::n
    real*8::c(size(b))

    n = size(b)

    if(n==0) then
      print *, "ERROR: linear solver called with zero-order matix!"
      stop
    elseif(n==2) then
       c(:) = linear_solver_n2(a(:, :), b(:))
    elseif(n==3) then
       c(:) = linear_solver_n3(a(:, :), b(:))
    else
       c(:) = linear_solver_nn(a(:, :), b(:), n)
    end if

  end function linear_solver

  ! ***********************
  ! solve Ax=B analytically for a 2-levels system with first row of ones
  ! (  1,   1) * (c1) = (b1)
  ! (a21, a22)   (c2)   (0 )
  function linear_solver_n2(a, b) result(c)
    implicit none
    integer,parameter::n=2
    real*8,intent(in)::a(n, n), b(n)
    real*8::c(n), iab

    c(1) = a(2, 2) * b(1) / (a(2, 2) - a(2, 1))
    c(2) = b(1) - c(1)

  end function linear_solver_n2

  ! ************************
  ! solve Ax=B analytically for a 3-levels system with first row of ones
  ! (  1,   1,   1) * (c1) = (b1)
  ! (a21, a22, a23)   (c2)   (0 )
  ! (a31, a32, a33)   (c3)   (0 )
  function linear_solver_n3(a, b) result(c)
    implicit none
    integer,parameter::n=3
    real*8,intent(in)::a(n, n), b(n)
    real*8::iab, c(n)

    iab = b(1) / (a(2, 1) * (a(3, 3) - a(3, 2)) &
      + a(2, 2) * (a(3, 1) - a(3, 3)) &
      + a(2, 3) * (a(3, 2) - a(3, 1)))

    c(1) = (a(2, 3) * a(3, 2) - a(2, 2) * a(3, 3)) * iab
    c(2) = (a(2, 1) * a(3, 3) - a(2, 3) * a(3, 1)) * iab
    c(3) = b(1) - c(1) - c(2)

  end function linear_solver_n3

  ! ************************
  ! solve Ax=B for a n-levels system
  function linear_solver_nn(a, b, n) result(c)
    implicit none
    integer,intent(in)::n
    real*8,intent(in)::a(n, n), b(n)
    real*8::c(n)
    integer::ipiv(n), ierr, i, j

    c(:) = b(:)
    print *, "LAPACK DGESV not set!"
    stop

    ! call dgesv(n, 1, a(:, :), n, ipiv(:), c(:), n, ierr)
    !
    ! if(ierr /= 0) then
    !    print *, "ERROR: dgesv error", ierr
    !    if(ierr < 0) then
    !      print *, "the ith argument had an illegal value, ith=", ierr
    !    else
    !      print *, "U(i, i) is exactly zero (singular matrix), i=", ierr
    !      do i=1,n
    !        do j=1,n
    !          if(a(i,j) /= 0d0) then
    !            print '(a10,2I5,E17.8e3)', "A(i, j):", i, j, a(i, j)
    !          end if
    !        end do
    !      end do
    !    end if
    !    print *, "NOTE: the next floating point exception is triggered to get backtrace"
    !    c = -1d0
    !    c = log10(c)
    !    stop
    ! end if

  end function linear_solver_nn

  ! *************************
  ! solve a linear system as
  ! (a11, -a12,    0,    0)   (c1)   (0)
  ! (  0,  a22, -a23,    0)   (c2)   (0)
  ! (  0,    0,  a33, -a34)   (c3)   (0)
  ! (  1,    1,    1,    1) * (c4) = (b)
  function bidiagonal_solver(a, b, n) result(c)
    implicit none
    integer,intent(in)::n
    real*8,intent(in)::a(n, n), b
    real*8::c(n), w(n)
    integer::i

    ! compute coefficients
    w(2) = a(1, 1) / a(1, 2)
    do i=3,n
      w(i) = w(i-1) * a(i-1, i-1) / a(i-1, i)
    end do

    ! find unknowns
    c(1) = b / (1d0 + sum(w(2:n)))
    c(2:n) = w(2:n) * c(1)

  end function bidiagonal_solver

end module prizmo_linear_solver

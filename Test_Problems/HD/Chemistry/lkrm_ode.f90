module lkrm_ode
contains

  !************************
  !evolve chemistry for a time-step dt (s)
  ! n(:) are species number densities
  subroutine dochem(n, Tgas, jflux, dt)
    use lkrm_commons
    use lkrm_rates
    implicit none
    real*8,intent(inout)::n(nmols)
    real*8,intent(in)::dt, Tgas, jflux(nphoto)

    !DLSODES variables
    integer,parameter::meth=2 !1=adam, 2=BDF
    integer,parameter::lrw=20+3*nmols**2+20*nmols
    integer,parameter::liw=30
    !tolerances
    real*8,parameter::rtol(nmols) = 1d-8
    real*8,parameter::atol(nmols) = 1d-20
    integer::neqa(1),itol,itask,istate,iopt,mf
    integer::iwork(liw)
    real*8::rwork(lrw),tloc

    iwork(:) = 0
    rwork(:) = 0d0
    itol = 4 !both tolerances are arrays

    !number of equations is a single sized array
    neqa(:) = nmols

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
    iwork(6) = int(1e5) !maximum number of iteration before warning
    iwork(5) = 0 !maximum integration order, zero means default
    MF = 222 !internally generated Jacobian and sparsity

    istate = 1
    tloc = 0d0
    do
       !call solver, see DLSODES documentation
       CALL DLSODES(fex, neqa(:), n(:), tloc, dt, &
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
       elseif(istate==-5) then
          !problem with sparsity, need to recompute
          istate = 3
          cycle
       else
          !unknonw problem stop program
          print *,istate
          stop
       end if
    end do

  end subroutine dochem

  !*********************
  !differential equations, returns dn(:)
  ! see DLSODES documentation
  subroutine fex(neq, tt, n, dn)
    use lkrm_commons
    implicit none
    integer::neq
    real*8::n(neq), tt, dn(neq)

    !!BEGIN_ODE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    dn(idx_H2) = -kall(1)*n(idx_H2) &
        +kall(2)*n(idx_H)*n(idx_H)

    dn(idx_H) = +kall(1)*n(idx_H2) &
        +kall(1)*n(idx_H2) &
        -kall(2)*n(idx_H)*n(idx_H) &
        -kall(2)*n(idx_H)*n(idx_H)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ODE

  end subroutine fex

  !***************************
  !Jacobian, pd(i,j)=df(i)/dx(j), see DLSODES documentation
  subroutine jes(neq, tt, n, j, ian, jan, pdj)
    use lkrm_commons
    implicit none
    integer::neq, j, ian, jan
    real*8::tt, n(neq), pdj(neq)

  end subroutine jes

end module lkrm_ode

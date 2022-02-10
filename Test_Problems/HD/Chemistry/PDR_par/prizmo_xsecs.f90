module prizmo_xsecs
contains
  ! **************************
  subroutine load_xsecs_all()
    use prizmo_commons, only: xsecs
    implicit none

    call load_energy_grid()

    ! initialize xsecs to zero
    xsecs(:, :) = 0d0

    !!BEGIN_LOAD_XSECS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    !H2 -> H + H
    call load_xsecs(258, "xsecs_runtime/H2__H_H.dat")

    !CO -> C + O
    call load_xsecs(259, "xsecs_runtime/CO__C_O.dat")

    !H2+ -> H+ + H
    call load_xsecs(260, "xsecs_runtime/H2j__Hj_H.dat")

    !H3+ -> H2+ + H
    call load_xsecs(261, "xsecs_runtime/H3j__H2j_H.dat")

    !H3+ -> H2 + H+
    call load_xsecs(262, "xsecs_runtime/H3j__H2_Hj.dat")

    !C -> C+ + E
    call load_xsecs(263, "xsecs_runtime/C__Cj_E.dat")

    !CH -> CH+ + E
    call load_xsecs(264, "xsecs_runtime/CH__CHj_E.dat")

    !CH -> C + H
    call load_xsecs(265, "xsecs_runtime/CH__C_H.dat")

    !CH+ -> C+ + H
    call load_xsecs(266, "xsecs_runtime/CHj__Cj_H.dat")

    !CH2 -> CH2+ + E
    call load_xsecs(267, "xsecs_runtime/CH2__CH2j_E.dat")

    !CH2 -> CH + H
    call load_xsecs(268, "xsecs_runtime/CH2__CH_H.dat")

    !CH2+ -> CH+ + H
    call load_xsecs(269, "xsecs_runtime/CH2j__CHj_H.dat")

    !CH3 -> CH3+ + E
    call load_xsecs(270, "xsecs_runtime/CH3__CH3j_E.dat")

    !CH3 -> CH2 + H
    call load_xsecs(271, "xsecs_runtime/CH3__CH2_H.dat")

    !CH3 -> CH + H2
    call load_xsecs(272, "xsecs_runtime/CH3__CH_H2.dat")

    !CH4 -> CH3 + H
    call load_xsecs(273, "xsecs_runtime/CH4__CH3_H.dat")

    !CH4 -> CH2 + H2
    call load_xsecs(274, "xsecs_runtime/CH4__CH2_H2.dat")

    !CH4 -> CH + H2 + H
    call load_xsecs(275, "xsecs_runtime/CH4__CH_H2_H.dat")

    !OH -> OH+ + E
    call load_xsecs(276, "xsecs_runtime/OH__OHj_E.dat")

    !OH -> O + H
    call load_xsecs(277, "xsecs_runtime/OH__O_H.dat")

    !OH+ -> O + H+
    call load_xsecs(278, "xsecs_runtime/OHj__O_Hj.dat")

    !H2O -> OH + H
    call load_xsecs(279, "xsecs_runtime/H2O__OH_H.dat")

    !H2O -> H2O+ + E
    call load_xsecs(280, "xsecs_runtime/H2O__H2Oj_E.dat")

    !O2 -> O + O
    call load_xsecs(281, "xsecs_runtime/O2__O_O.dat")

    !O2 -> O2+ + E
    call load_xsecs(282, "xsecs_runtime/O2__O2j_E.dat")

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_LOAD_XSECS

  end subroutine load_xsecs_all

  ! ************************
  ! load xsecs from file_name
  subroutine load_energy_grid()
    use prizmo_commons
    implicit none
    integer::unit, ios, i


    ! open file to read
    open(newunit=unit, file="runtime_data/energy_grid.dat", status="old", iostat=ios)

    ! check if file exists
    if(ios/=0) then
      print *, "ERROR: missing energy grid file"
      stop
    end if

    ! loop to read energy grid
    do i=1,nphoto
      read(unit, *) energy_grid(i)
      inv_energy_grid(i) = 1d0 / energy_grid(i)
    end do
    close(unit)

  end subroutine load_energy_grid

  ! ************************
  ! load xsecs from file_name
  subroutine load_xsecs(idx, file_name, verbose)
    use prizmo_commons
    implicit none
    integer,intent(in)::idx
    integer,intent(in),optional::verbose
    character(len=*),intent(in)::file_name
    integer::unit, ios, i
    real*8::delta_trapz(nphoto), energy


    if(present(verbose)) then
      if(verbose>0) then
        print *, "loading ", trim(file_name)
      end if
    end if

    ! open file to read
    open(newunit=unit, file=trim(file_name), status="old", iostat=ios)

    ! check if file exists
    if(ios/=0) then
      print *, "ERROR: missing file", trim(file_name)
      stop
    end if

    ! loop to read xsecs
    do i=1,nphoto
      read(unit, *) energy, xsecs(i, idx)
    end do
    close(unit)

    ! store energy difference
    delta_energy_grid(:) = energy_grid(2:nphoto) - energy_grid(1:nphoto-1)

    delta_trapz(1) = energy_grid(2) - energy_grid(1)
    ! FIXME delta_trapz(2:nphoto-1) = delta_energy_grid(2:nphoto-1) + delta_energy_grid(1:nphoto-2)
    delta_trapz(2:nphoto-1) = energy_grid(3:nphoto) - energy_grid(1:nphoto-2)
    delta_trapz(nphoto) = energy_grid(nphoto) - energy_grid(nphoto-1)

    xsecs_trapz(:, idx) = delta_trapz(:) * xsecs(:, idx) / energy_grid(:)
    xsecs_trapz_heat(:, idx) = delta_trapz(:) * xsecs(:, idx)

    ! search grid points indexes for Draine flux limits (default is full range)
    imin_fDraine = 1
    imax_fDraine = nphoto
    do i=1,nphoto-1
      if(energy_grid(i) <= 6d0 .and. energy_grid(i+1) >= 6d0) then
        imin_fDraine = i
      end if

      if(energy_grid(i) <= 13.6d0 .and. energy_grid(i+1) >= 13.6d0) then
        imax_fDraine = i
        exit  ! energy is ascending, so here no need to continue the loop
      end if
    end do

  end subroutine load_xsecs
end module prizmo_xsecs

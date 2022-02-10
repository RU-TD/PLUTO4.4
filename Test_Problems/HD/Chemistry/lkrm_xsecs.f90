module lkrm_xsecs
contains
  ! **************************
  subroutine load_xsecs_all()
    use lkrm_commons, only: xsecs
    implicit none

    ! initialize xsecs to zero
    xsecs(:, :) = 0d0

    !!BEGIN_LOAD_XSECS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_LOAD_XSECS

  end subroutine load_xsecs_all

  ! ************************
  ! load xsecs from file_name
  subroutine load_xsecs(idx, file_name)
    use lkrm_commons, only: xsecs, nphoto, energy_grid, inv_energy_grid
    implicit none
    integer,intent(in)::idx
    character(len=*),intent(in)::file_name
    integer::unit, ios, i

    print *, "loading ", trim(file_name)

    ! open file to read
    open(newunit=unit, file=trim(file_name), status="old", iostat=ios)

    ! check if file exists
    if(ios/=0) then
       print *, "ERROR: missing file", trim(file_name)
       stop
    end if

    ! loop to read
    do i=1,nphoto
       read(unit, *) energy_grid(i), xsecs(i, idx)
       inv_energy_grid(i) = 1d0 / energy_grid(i)
    end do
    close(unit)

  end subroutine load_xsecs
end module lkrm_xsecs

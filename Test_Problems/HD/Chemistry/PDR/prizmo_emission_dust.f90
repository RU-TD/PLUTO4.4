module prizmo_emission_dust
  use prizmo_commons
!!BEGIN_EMISSION_DUST_COMMONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2021-01-20 16:00:19
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    real*8,parameter::em_xmin=0.000000e+00
    real*8,parameter::em_invdx=8.627827e+00
    real*8,parameter::em_fact=8.340233e+00
    integer,parameter::nsteps=30

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!END_EMISSION_DUST_COMMONS

  real*8::em_data(nphoto, nsteps), xdata(nsteps)

contains
  ! ************************
  ! load emission data from file as a function of photo bins and Tdust
  subroutine load_data()
    implicit none
    integer::i, j, unit, idx ,ios
    real*8::energy

    open(newunit=unit, file="runtime_data/dust_emission_bins.dat", status="old", iostat=ios)
    if(ios /= 0) then
       print *, "ERROR: problems while loading dust emission bins!"
       stop
    end if
    ! loop on temperature steps
    do j=1,nsteps
       ! loop on photobins
       do i=1,nphoto
         read(unit, *) xdata(j), idx, energy, em_data(i, j)
       end do
    end do
    close(unit)

    xdata(:) = log10(xdata(:))
    em_data(:, :) = log10(em_data(:, :))

  end subroutine load_data

  ! ***************************
  ! add dust emission to the emission array
  subroutine add_emission(emission, Tdust)
    use prizmo_commons
    implicit none
    real*8,intent(inout)::emission(nphoto)
    real*8,intent(in)::Tdust
    real*8::logTdust, pre
    integer::idx

    logTdust = log10(Tdust)
    idx = int((logTdust - em_xmin) * em_fact) + 1

    pre = (logTdust - xdata(idx)) * em_invdx
    emission(:) = emission(:) + 1d1**(pre * (em_data(:, idx+1) - em_data(:, idx)) + em_data(:, idx))

  end subroutine add_emission

end module prizmo_emission_dust

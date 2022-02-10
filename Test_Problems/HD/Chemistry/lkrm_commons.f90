module lkrm_commons

  !constants
  real*8,parameter::spy=365d0*24d0*3600d0  ! seconds per year
  real*8,parameter::pi=acos(-1d0)
  real*8,parameter::pi43=4d0*pi/3d0
  real*8,parameter::ev2erg=1.602176487d-12  ! eV -> erg
  real*8,parameter::hplanck=6.6260755d-27  ! erg*s
  real*8,parameter::hplanck_eV=hplanck/ev2erg  ! eV*s

  !!BEGIN_ARRAYSIZE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! number of species
    integer,parameter::nmols=2
    ! number of reactions
    integer,parameter::nrea=2

    ! number of energy bins
    integer,parameter::nphoto=10000

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_ARRAYSIZE

  ! rate coefficients
  real*8::kall(nrea)

  ! cross-section, cm2
  real*8::xsecs(nphoto, nrea)

  ! energy grid points, eV
  real*8::energy_grid(nphoto)

  ! inverse energy grid points, 1/eV
  real*8::inv_energy_grid(nphoto)

  !!BEGIN_IDXLIST
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    integer,parameter::idx_H2=1
    integer,parameter::idx_H=2

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_IDXLIST

  !!BEGIN_SPECIES_MASS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2018-03-08 13:34:08
    ! CHANGESET: 3f0f625
    ! URL: https://GiovanniPicogna@bitbucket.org/tgrassi/mocassin_xray_chem.git
    ! BY: picogna@picogna-laptop

    ! species mass, g
    real*8,parameter::mass(nmols) = (/&
        3.3452438d-24, &
        1.6726219d-24/)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_SPECIES_MASS

end module lkrm_commons

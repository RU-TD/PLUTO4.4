module prizmo_commons

  !constants
  real*8,parameter::spy=365d0*24d0*3600d0  ! seconds per year
  real*8,parameter::pi=acos(-1d0)
  real*8,parameter::pi43=4d0*pi/3d0
  real*8,parameter::ev2erg=1.602176487d-12  ! eV -> erg
  real*8,parameter::hplanck=6.6260755d-27  ! erg*s
  real*8,parameter::hplanck_eV=hplanck/ev2erg  ! eV*s
  real*8,parameter::kboltzmann=1.38064852d-16  ! erg/K
  real*8,parameter::kboltzmann_eV=kboltzmann/ev2erg  ! eV/K
  real*8,parameter::clight=2.99792458d10  ! cm/s
  real*8,parameter::stefan_boltzmann=5.6704d-5  ! erg/cm2/s/K4
  real*8,parameter::stefan_boltzmann_eV=stefan_boltzmann/ev2erg  ! eV/cm2/s/K4
  real*8,parameter::proton_mass=1.6726219d-24  ! g
  integer,parameter::max_character_len=30

  !!BEGIN_ARRAYSIZE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! number of species
    integer,parameter::nmols=31
    ! number of reactions
    integer,parameter::nrea=289
    ! number of energy bins
    integer,parameter::nphoto=1000

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_ARRAYSIZE

  !!BEGIN_USER_VARIABLES
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! user variables found in the network or in context object
    real*8::variable_Av
    real*8::variable_G0
    real*8::variable_crflux
    real*8::variable_NCO_incoming
    real*8::variable_NCO_escaping
    real*8::variable_NH2_incoming
    real*8::variable_NH2_escaping
    real*8::variable_NH2O_escaping

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_USER_VARIABLES

  ! cooling and heating mode flags, to change cooling configurations, defualt is zero
  integer::cooling_mode
  integer::heating_mode
  integer::tdust_mode
  integer::chemistry_mode
  integer::beta_escape_mode

  ! service variable, ignore for standard use
  integer::debug

  ! dust variables
  real*8,parameter::d2g = 1d-2  ! dust/gas mass ratio
  real*8,parameter::dust_d2g = d2g  ! alias
  real*8,parameter::dust_rho_bulk = 3d0  ! bulk density, g/cm3
  real*8,parameter::dust_amin = 5d-7  ! minimum size, cm
  real*8,parameter::dust_amax = 2.5d-5  ! maximum size, cm
  real*8,parameter::dust_pexp = -3.5  ! power-law dust distribution exponent
  real*8,parameter::dust_p1 = dust_pexp + 1d0  ! power-law exponent for number density integrals
  real*8,parameter::dust_p3 = dust_pexp + 3d0  ! power-law exponent for area integrals
  real*8,parameter::dust_p4 = dust_pexp + 4d0  ! power-law exponent for mass integrals
  real*8,parameter::dust_app2 = (3d-8)**2  !sites separation area, cm2  ! exponents
  real*8,parameter::dust_debye_nu = 1d12  ! Debye frequency, 1/s
  ! integrals with a^0, a^2, and a^3 as kernel
  real*8,parameter::dust_fact1 = (dust_amax**dust_p1 - dust_amin**dust_p1) / dust_p1
  real*8,parameter::dust_fact3 = (dust_amax**dust_p3 - dust_amin**dust_p3) / dust_p3
  real*8,parameter::dust_fact4 = (dust_amax**dust_p4 - dust_amin**dust_p4) / dust_p4

  !!BEGIN_DUST_CHARGE_COMMONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! min dust charge
    integer,parameter::dust_zmin=-4
    ! max dust charge
    integer,parameter::dust_zmax=4
    ! dust charge bins
    ! integer,parameter::dust_ncharge=9

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_DUST_CHARGE_COMMONS

  ! dust charge bins
  integer,parameter::dust_ncharge = dust_zmax - dust_zmin + 1


  ! line velocity b-factor, cm/s
  real*8,parameter::b_line = 1d5

  ! velocity gradient escape probability, 1/s
  real*8,parameter::dvdz = 1d-5

  ! ortho-para ratio H2
  real*8,parameter::opratio_H2 = 1d0 / 3d0

  ! runtime folder where data required at runtime are stored
  character(len=13)::runtime_folder="runtime_data/"

  ! forced temperature limits, K
  ! default are 1e99, -1e99
  real*8::Tgas_min_forced
  real*8::Tgas_max_forced

  !!BEGIN_TEMPERATURE_LIMITS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! min, max clipping Tgas for current problem, K
    real*8,parameter::Tgas_min = 3.000000e+00
    real*8,parameter::Tgas_max = 1.000000e+08

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_TEMPERATURE_LIMITS

  ! rate coefficients
  real*8::kall(nrea), krate_evaluate_once(nrea)

  ! cross-section, cm2
  real*8::xsecs(nphoto, nrea)

  ! integration kernels of trapezoidal method for photo rates integral
  ! xsecs * dE / E
  real*8::xsecs_trapz(nphoto, nrea)
  ! analogous for heating xsecs * dE
  real*8::xsecs_trapz_heat(nphoto, nrea)

  ! energy grid points, eV
  real*8::energy_grid(nphoto)

  ! indexes in energy_grid of Draine flux validity
  integer::imin_fDraine, imax_fDraine

  ! inverse energy grid points, 1/eV
  real*8::inv_energy_grid(nphoto)

  ! difference between i and i+1 energy grid points, eV
  real*8::delta_energy_grid(nphoto-1)

  ! extra variables index
  integer,parameter::idx_Tgas = nmols + 1

  !!BEGIN_IDXLIST
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    integer,parameter::idx_H=1
    integer,parameter::idx_CH=2
    integer,parameter::idx_C=3
    integer,parameter::idx_H2=4
    integer,parameter::idx_CH3=5
    integer,parameter::idx_CH2=6
    integer,parameter::idx_CH4=7
    integer,parameter::idx_OH=8
    integer,parameter::idx_O=9
    integer,parameter::idx_H2O=10
    integer,parameter::idx_CO=11
    integer,parameter::idx_O2=12
    integer,parameter::idx_CH2j=13
    integer,parameter::idx_CHj=14
    integer,parameter::idx_CH3j=15
    integer,parameter::idx_Hej=16
    integer,parameter::idx_He=17
    integer,parameter::idx_Hj=18
    integer,parameter::idx_Cj=19
    integer,parameter::idx_Oj=20
    integer,parameter::idx_H2j=21
    integer,parameter::idx_COj=22
    integer,parameter::idx_E=23
    integer,parameter::idx_H3j=24
    integer,parameter::idx_CH4j=25
    integer,parameter::idx_OHj=26
    integer,parameter::idx_CH5j=27
    integer,parameter::idx_H2Oj=28
    integer,parameter::idx_H3Oj=29
    integer,parameter::idx_HCOj=30
    integer,parameter::idx_O2j=31

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_IDXLIST

  !!BEGIN_SPECIES_MASS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    ! species mass, g
    real*8,parameter::mass(nmols) = (/&
        1.6726219d-24, &
        2.17440847d-23, &
        2.00714628d-23, &
        3.3452438d-24, &
        2.50893285d-23, &
        2.34167066d-23, &
        2.67619504d-23, &
        2.84345723d-23, &
        2.67619504d-23, &
        3.01071942d-23, &
        4.68334132d-23, &
        5.35239008d-23, &
        2.34157957d-23, &
        2.17431738d-23, &
        2.50884176d-23, &
        6.68957666d-24, &
        6.6904876d-24, &
        1.67171096d-24, &
        2.00705519d-23, &
        2.67610395d-23, &
        3.34433286d-24, &
        4.68325023d-23, &
        9.10442514d-28, &
        5.01695476d-24, &
        2.67610395d-23, &
        2.84336614d-23, &
        2.84336614d-23, &
        3.01062833d-23, &
        3.17789052d-23, &
        4.85051242d-23, &
        5.35229899d-23/)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_SPECIES_MASS

end module prizmo_commons

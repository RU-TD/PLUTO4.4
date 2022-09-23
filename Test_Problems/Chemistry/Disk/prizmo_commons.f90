module prizmo_commons
use prizmo_fit
implicit none

!! PREPROCESS_COMMON_VARS

integer,parameter::nspecies=40
integer,parameter::nphoto=1000
integer,parameter::nreactions=315

integer,parameter::idx_C=1
integer,parameter::idx_CH=2
integer,parameter::idx_CH2=3
integer,parameter::idx_CH2j=4
integer,parameter::idx_CH3=5
integer,parameter::idx_CH3j=6
integer,parameter::idx_CH4=7
integer,parameter::idx_CH4j=8
integer,parameter::idx_CH5j=9
integer,parameter::idx_CHj=10
integer,parameter::idx_CO=11
integer,parameter::idx_CO_DUST=12
integer,parameter::idx_COj=13
integer,parameter::idx_Cj=14
integer,parameter::idx_Cjj=15
integer,parameter::idx_Cjjj=16
integer,parameter::idx_Cjjjj=17
integer,parameter::idx_E=18
integer,parameter::idx_H=19
integer,parameter::idx_H2=20
integer,parameter::idx_H2O=21
integer,parameter::idx_H2O_DUST=22
integer,parameter::idx_H2Oj=23
integer,parameter::idx_H2j=24
integer,parameter::idx_H3Oj=25
integer,parameter::idx_H3j=26
integer,parameter::idx_HCOj=27
integer,parameter::idx_He=28
integer,parameter::idx_Hej=29
integer,parameter::idx_Hejj=30
integer,parameter::idx_Hj=31
integer,parameter::idx_O=32
integer,parameter::idx_O2=33
integer,parameter::idx_O2j=34
integer,parameter::idx_OH=35
integer,parameter::idx_OHj=36
integer,parameter::idx_Oj=37
integer,parameter::idx_Ojj=38
integer,parameter::idx_Ojjj=39
integer,parameter::idx_Ojjjj=40

!! PREPROCESS_END

integer,parameter::idx_Tgas=nspecies + 1
integer,parameter::zmin=-2, zmax=2

integer,parameter::jtab_fit_nt=30
real*8::jion_fit_data(zmin:zmax, jtab_fit_nt), jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx
real*8::jele_fit_data(zmin:zmax, jtab_fit_nt)
real*8::jion_cool_fit_data(zmin:zmax, jtab_fit_nt)
real*8::jele_cool_fit_data(zmin:zmax, jtab_fit_nt)

!! PREPROCESS_ATOMIC_COOLING_NVEC
!! PREPROCESS_END

integer,parameter::atomic_cooling_n1=1000  ! tgas, K
integer,parameter::atomic_cooling_3lev_nvec=6  ! 3 levels
integer,parameter::atomic_cooling_2lev_nvec=2  ! 2 levels

!! PREPROCESS_ATOMIC_COOLING_COMMONS

type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_C_Hj
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_C_H
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_C_e
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_C_H2or
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_C_H2pa
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_O_H
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_O_Hj
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_O_e
type(fit1d_data_vec(nv=atomic_cooling_2lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_Cj_e
type(fit1d_data_vec(nv=atomic_cooling_2lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_Cj_H
type(fit1d_data_vec(nv=atomic_cooling_3lev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_Oj_e

!! PREPROCESS_END

! ! atomic cooling table, 2D
! integer,parameter::atomic_cooling2d_n1=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling2d_n2=100  ! tgas, K
! type(fit2d_data_vec(nv=atomic_cooling_nvec2d, n1=atomic_cooling2d_n1, n2=atomic_cooling2d_n2))::atomic_cooling_table_2d
!
! ! atomic cooling table, 3D
! integer,parameter::atomic_cooling3d_n1=30  ! H abundance, cm-3
! integer,parameter::atomic_cooling3d_n2=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling3d_n3=100  ! tgas, K
! type(fit3d_data_vec(nv=atomic_cooling_nvec3d, n1=atomic_cooling3d_n1, n2=atomic_cooling3d_n2, &
! n3=atomic_cooling3d_n3))::atomic_cooling_table_3d
!
! ! atomic cooling table, 4D
! integer,parameter::atomic_cooling4d_n1=30  ! H abundance, cm-3
! integer,parameter::atomic_cooling4d_n2=30  ! H+ abundance, cm-3
! integer,parameter::atomic_cooling4d_n3=30  ! e- abundance, cm-3
! integer,parameter::atomic_cooling4d_n4=100  ! tgas, K
! type(fit4d_data_vec(nv=atomic_cooling_nvec4d, n1=atomic_cooling4d_n1, n2=atomic_cooling4d_n2, &
!                     n3=atomic_cooling4d_n3, n4=atomic_cooling4d_n4))::atomic_cooling_table_4d

! dust cooling table, 3D
integer,parameter::dust_cooling_table_n1=100
integer,parameter::dust_cooling_table_n2=100
integer,parameter::dust_cooling_table_n3=50
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::dust_cooling_table_data
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::dust_heating_table_data
type(fit3d_data(n1=dust_cooling_table_n1, n2=dust_cooling_table_n2, n3=dust_cooling_table_n3))::tdust_table_data

! dust cooling pre-jflux
real*8::pre_dust_cooling_table(nphoto)

! photoelectic heating Jpe table
real*8::jpe_table(nphoto, zmin:zmax)
real*8::jpe_heating_table(nphoto, zmin:zmax)
real*8::phterm_jpe(zmin:zmax), phterm_jpe_heating(zmin:zmax)

! H2 cooling tables
integer,parameter::cool_H2_vec=5
integer,parameter::cool_H2_nx=1000
type(fit1d_data_vec(nv=cool_H2_vec, n1=cool_H2_nx))::cooling_H2data

! dust kabs table, cm2/g
real*8::dust_kappa_opacity(nphoto)

! shielding H2 tables
integer,parameter::shielding_H2_n1=30, shielding_H2_n2=100
type(fit2d_data(n1=shielding_H2_n1, n2=shielding_H2_n2))::shielding_H2_data

! shielding CO tables
integer,parameter::shielding_CO_n1=50, shielding_CO_n2=50
type(fit2d_data(n1=shielding_CO_n1, n2=shielding_CO_n2))::shielding_CO_data

! cooling CO tables
integer,parameter::cool_CO_tab_n1=40, cool_CO_tab_n2=40, cool_CO_tab_n3=40
type(fit3d_data(n1=cool_CO_tab_n1, n2=cool_CO_tab_n2, n3=cool_CO_tab_n3))::cool_CO_tab_data


real*8::kall(nreactions)
real*8::kall_heat(nreactions)
real*8::photo_xsecs(nphoto, nreactions)
real*8::energy_threshold(nphoto)

real*8::radial_Ncol_H2
real*8::radial_Ncol_CO
real*8::vertical_Ncol_CO

real*8::energy(nphoto)
real*8::delta_energy(nphoto-1)

real*8::ode_atol(nspecies+1)
real*8::ode_rtol(nspecies+1)

real*8::gamma_ad, d2g, user_Av, user_cr, ortho_to_para
real*8::chi_FUV  ! habing flux in range 912-1100 AA
real*8::rho_gas, rho_dust  ! gas and dust mass densities, g/cm3, do not change during integration


!! PREPROCESS_RADIATION_CONSTANTS

! FUV interval photobins indexes
integer,parameter::fuv_idx1=2
integer,parameter::fuv_idx2=159
real*8,parameter::habing_flux=1.273584d-18

!! PREPROCESS_END

!! PREPROCESS_KABS_INTEGRAL

! kabs integral in FUV range
real*8,parameter::kabs_integral=6.308145d+04

!! PREPROCESS_END

real*8,parameter::d2g_min=1d-8  ! below this limit no dust cooling/heating are calculated

real*8,parameter::hplanck=6.6260755d-27  ! erg * s
real*8,parameter::kboltzmann=1.380658d-16  ! erg / K
real*8,parameter::erg2ev=6.24150647996e11  ! 1 erg in eV
real*8,parameter::ev2erg=1d0 / erg2ev  ! 1 eV in erg
real*8,parameter::pmass=1.6726219d-24  ! proton mass in g
real*8,parameter::pi=acos(-1d0)

logical::solve_thermo
logical::solve_chemistry

character(len=1024)::runtime_data_folder
character(len=128)::reactions_verbatim(nreactions)

! FIXME: could cause problems with OPENMP
real*8::jflux_common(nphoto), log_Eabsorption

!! PREPROCESS_MASSES

real*8,parameter::masses(nspecies) = (/2.007146d-23, &
	2.341671d-23, &
	2.676195d-23, &
	2.676104d-23, &
	3.010719d-23, &
	3.010628d-23, &
	3.345244d-23, &
	3.345153d-23, &
	3.679677d-23, &
	2.341580d-23, &
	4.683341d-23, &
	4.683341d-23, &
	4.683250d-23, &
	2.007055d-23, &
	2.006964d-23, &
	2.006873d-23, &
	2.006782d-23, &
	9.109380d-28, &
	3.345244d-24, &
	6.690488d-24, &
	3.345244d-23, &
	3.345244d-23, &
	3.345153d-23, &
	6.689577d-24, &
	3.679677d-23, &
	1.003482d-23, &
	5.017775d-23, &
	6.690488d-24, &
	6.689577d-24, &
	6.688666d-24, &
	3.344333d-24, &
	2.676195d-23, &
	5.352390d-23, &
	5.352299d-23, &
	3.010719d-23, &
	3.010628d-23, &
	2.676104d-23, &
	2.676013d-23, &
	2.675922d-23, &
	2.675831d-23/)

!! PREPROCESS_END

end module prizmo_commons

module prizmo_loaders
  use prizmo_commons
  use prizmo_fit
contains

  ! ******************
  subroutine load_energy()
    implicit none
    integer::i, unit

    print *, "loading energy..."

    open(newunit=unit, file=trim(runtime_data_folder)//"energy.dat", status="old")
    do i=1,nphoto
      read(unit, *) energy(i)
    end do
    close(unit)

    delta_energy(:) = energy(2:nphoto) - energy(1:nphoto-1)

  end subroutine

  ! ******************
  subroutine load_all_photo_xsecs()

    print *, "loading photo xsecs..."

    !! PREPROCESS_LOAD_XSECS

    photo_xsecs(:, 283) = load_photo_xsecs("photo_xsecs_H__H+_E.dat")
    photo_xsecs(:, 284) = load_photo_xsecs("photo_xsecs_He__He+_E.dat")
    photo_xsecs(:, 285) = load_photo_xsecs("photo_xsecs_He+__He++_E.dat")
    photo_xsecs(:, 286) = load_photo_xsecs("photo_xsecs_O__O+_E.dat")
    photo_xsecs(:, 287) = load_photo_xsecs("photo_xsecs_O+__O++_E.dat")
    photo_xsecs(:, 288) = load_photo_xsecs("photo_xsecs_O++__O+++_E.dat")
    photo_xsecs(:, 289) = load_photo_xsecs("photo_xsecs_O+++__O++++_E.dat")
    photo_xsecs(:, 290) = load_photo_xsecs("photo_xsecs_C__C+_E.dat")
    photo_xsecs(:, 291) = load_photo_xsecs("photo_xsecs_C+__C++_E.dat")
    photo_xsecs(:, 292) = load_photo_xsecs("photo_xsecs_C++__C+++_E.dat")
    photo_xsecs(:, 293) = load_photo_xsecs("photo_xsecs_C+++__C++++_E.dat")
    photo_xsecs(:, 294) = load_photo_xsecs("photo_xsecs_H2__H_H.dat")
    photo_xsecs(:, 295) = load_photo_xsecs("photo_xsecs_CO__C_O.dat")
    photo_xsecs(:, 296) = load_photo_xsecs("photo_xsecs_H2+__H+_H.dat")
    photo_xsecs(:, 297) = load_photo_xsecs("photo_xsecs_CH__CH+_E.dat")
    photo_xsecs(:, 298) = load_photo_xsecs("photo_xsecs_CH__C_H.dat")
    photo_xsecs(:, 299) = load_photo_xsecs("photo_xsecs_CH+__C+_H.dat")
    photo_xsecs(:, 300) = load_photo_xsecs("photo_xsecs_CH2__CH_H.dat")
    photo_xsecs(:, 301) = load_photo_xsecs("photo_xsecs_CH2+__CH+_H.dat")
    photo_xsecs(:, 302) = load_photo_xsecs("photo_xsecs_CH3__CH3+_E.dat")
    photo_xsecs(:, 303) = load_photo_xsecs("photo_xsecs_CH3__CH2_H.dat")
    photo_xsecs(:, 304) = load_photo_xsecs("photo_xsecs_CH4__CH3_H.dat")
    photo_xsecs(:, 305) = load_photo_xsecs("photo_xsecs_OH__O_H.dat")
    photo_xsecs(:, 306) = load_photo_xsecs("photo_xsecs_OH+__O_H+.dat")
    photo_xsecs(:, 307) = load_photo_xsecs("photo_xsecs_H2O__OH_H.dat")
    photo_xsecs(:, 308) = load_photo_xsecs("photo_xsecs_H2O__H2O+_E.dat")
    photo_xsecs(:, 309) = load_photo_xsecs("photo_xsecs_O2__O_O.dat")
    photo_xsecs(:, 310) = load_photo_xsecs("photo_xsecs_O2__O2+_E.dat")

    !! PREPROCESS_END

  end subroutine load_all_photo_xsecs

  ! ******************
  subroutine load_energy_thresholds()
    implicit none
    integer::i, unit

    print *, "loading energy thresholds..."

    energy_threshold(:) = 0d0

    open(newunit=unit, file=trim(runtime_data_folder)//"energy_thresholds.dat", status="old")
    do i=1,nreactions
      read(unit, *) energy_threshold(i)
    end do
    close(unit)

    ! eV -> erg
    energy_threshold = energy_threshold * ev2erg

  end subroutine load_energy_thresholds

  ! *****************
  function load_photo_xsecs(fname) result(xsecs)
    implicit none
    character(len=*)::fname
    integer::i, unit
    real*8::xsecs(nphoto)

    open(newunit=unit, file=trim(runtime_data_folder)//trim(fname), status="old")
    do i=1,nphoto
      read(unit, *) xsecs(i)
    end do
    close(unit)

  end function load_photo_xsecs

  ! ****************************
  subroutine load_all_atomic_cooling_tables()
    implicit none
    ! character(len=1024)::fnames2d(atomic_cooling_nvec2d)
    ! character(len=1024)::fnames3d(atomic_cooling_nvec3d)
    ! character(len=1024)::fnames4d(atomic_cooling_nvec4d)
    character(len=1024)::fnames_3lev(6), fnames_2lev(2)

    print *, "loading atomic cooling..."

    !! PREPROCESS_LOAD_ATOMIC_COOLING

    fnames_3lev = (/"runtime_data/cool_C_H+_k01k02.dat", &
    "runtime_data/cool_C_H+_k10.dat", &
    "runtime_data/cool_C_H+_k20.dat", &
    "runtime_data/cool_C_H+_k02.dat", &
    "runtime_data/cool_C_H+_k12.dat", &
    "runtime_data/cool_C_H+_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_C_Hj, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_C_H_k01k02.dat", &
    "runtime_data/cool_C_H_k10.dat", &
    "runtime_data/cool_C_H_k20.dat", &
    "runtime_data/cool_C_H_k02.dat", &
    "runtime_data/cool_C_H_k12.dat", &
    "runtime_data/cool_C_H_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_C_H, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_C_e_k01k02.dat", &
    "runtime_data/cool_C_e_k10.dat", &
    "runtime_data/cool_C_e_k20.dat", &
    "runtime_data/cool_C_e_k02.dat", &
    "runtime_data/cool_C_e_k12.dat", &
    "runtime_data/cool_C_e_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_C_e, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_C_H2or_k01k02.dat", &
    "runtime_data/cool_C_H2or_k10.dat", &
    "runtime_data/cool_C_H2or_k20.dat", &
    "runtime_data/cool_C_H2or_k02.dat", &
    "runtime_data/cool_C_H2or_k12.dat", &
    "runtime_data/cool_C_H2or_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_C_H2or, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_C_H2pa_k01k02.dat", &
    "runtime_data/cool_C_H2pa_k10.dat", &
    "runtime_data/cool_C_H2pa_k20.dat", &
    "runtime_data/cool_C_H2pa_k02.dat", &
    "runtime_data/cool_C_H2pa_k12.dat", &
    "runtime_data/cool_C_H2pa_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_C_H2pa, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_O_H_k01k02.dat", &
    "runtime_data/cool_O_H_k10.dat", &
    "runtime_data/cool_O_H_k20.dat", &
    "runtime_data/cool_O_H_k02.dat", &
    "runtime_data/cool_O_H_k12.dat", &
    "runtime_data/cool_O_H_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_O_H, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_O_H+_k01k02.dat", &
    "runtime_data/cool_O_H+_k10.dat", &
    "runtime_data/cool_O_H+_k20.dat", &
    "runtime_data/cool_O_H+_k02.dat", &
    "runtime_data/cool_O_H+_k12.dat", &
    "runtime_data/cool_O_H+_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_O_Hj, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_O_e_k01k02.dat", &
    "runtime_data/cool_O_e_k10.dat", &
    "runtime_data/cool_O_e_k20.dat", &
    "runtime_data/cool_O_e_k02.dat", &
    "runtime_data/cool_O_e_k12.dat", &
    "runtime_data/cool_O_e_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_O_e, do_log=.false.)
    
    fnames_2lev = (/"runtime_data/cool_C+_e_k01.dat", &
    "runtime_data/cool_C+_e_k10.dat"/)
    
    call load_1d_fit_vec(fnames_2lev, atomic_cooling_2lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_Cj_e, do_log=.false.)
    
    fnames_2lev = (/"runtime_data/cool_C+_H_k01.dat", &
    "runtime_data/cool_C+_H_k10.dat"/)
    
    call load_1d_fit_vec(fnames_2lev, atomic_cooling_2lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_Cj_H, do_log=.false.)
    
    fnames_3lev = (/"runtime_data/cool_O+_e_k01k02.dat", &
    "runtime_data/cool_O+_e_k10.dat", &
    "runtime_data/cool_O+_e_k20.dat", &
    "runtime_data/cool_O+_e_k02.dat", &
    "runtime_data/cool_O+_e_k12.dat", &
    "runtime_data/cool_O+_e_k21k20.dat"/)
    
    call load_1d_fit_vec(fnames_3lev, atomic_cooling_3lev_nvec, atomic_cooling_n1, &
    atomic_cooling_table_Oj_e, do_log=.false.)

    !! PREPROCESS_END

    ! call load_2d_fit_vec(fnames2d, atomic_cooling_nvec2d, atomic_cooling2d_n1, atomic_cooling2d_n2, &
    !   atomic_cooling_table_2d, do_log=.true.)
    !
    ! call load_3d_fit_vec(fnames3d, atomic_cooling_nvec3d, atomic_cooling3d_n1, atomic_cooling3d_n2, atomic_cooling3d_n3, &
    !   atomic_cooling_table_3d, do_log=.true.)
    !
    ! call load_4d_fit_vec(fnames4d, atomic_cooling_nvec4d, atomic_cooling4d_n1, atomic_cooling4d_n2, atomic_cooling4d_n3, atomic_cooling4d_n4, &
    !   atomic_cooling_table_4d, do_log=.true.)

  end subroutine load_all_atomic_cooling_tables

  ! ! ****************************
  ! subroutine load_atomic_cooling_table(fname, idx_atom)
  !   implicit none
  !   character(len=*),intent(in)::fname
  !   integer,intent(in)::idx_atom
  !   !real*8::fdata(nxp_cooling_table, nxe_cooling_table, ntgas_cooling_table)
  !
  !   ! print *, "loading atomic cooling table from "//trim(fname)
  !   ! call load_data_3d(trim(fname), fdata, &
  !   !    atomic_cooling_xp, nxp_cooling_table, atomic_cooling_xp_min, atomic_cooling_xp_fact, atomic_cooling_xp_invdx, &
  !   !    atomic_cooling_xe, nxe_cooling_table, atomic_cooling_xe_min, atomic_cooling_xe_fact, atomic_cooling_xe_invdx, &
  !   !    atomic_cooling_tgas, ntgas_cooling_table, atomic_cooling_tgas_min, atomic_cooling_tgas_fact, atomic_cooling_tgas_invdx, &
  !   !    do_log=.true.)
  !   ! atomic_cooling_table(:, :, :, idx_atom) = fdata
  !
  ! end subroutine load_atomic_cooling_table

  ! ***************************
  subroutine load_photoelectric_tables()
    implicit none
    integer::unit, i, j, z

    print *, "load phe tables..."
    open(newunit=unit, file=trim(runtime_data_folder)//"jptot.dat", status="old")
    do i=zmin, zmax
      do j=1,nphoto
        read(unit, *) z, jpe_table(j, i)
      end do
    end do
    close(unit)

    open(newunit=unit, file=trim(runtime_data_folder)//"jptot_heating.dat", status="old")
    do i=zmin, zmax
      do j=1,nphoto
        read(unit, *) z, jpe_heating_table(j, i)
      end do
    end do
    close(unit)

    ! ion-grain rates
    call load_jtab(trim(runtime_data_folder)//"jion.dat", jtab_fit_nt, jion_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx)
    call load_jtab(trim(runtime_data_folder)//"jelectron.dat", jtab_fit_nt, jele_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx)

    ! ion-grain recombination cooling
    call load_jtab(trim(runtime_data_folder)//"jion_cooling.dat", jtab_fit_nt, jion_cool_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx)
    call load_jtab(trim(runtime_data_folder)//"jelectron_cooling.dat", jtab_fit_nt, jele_cool_fit_data, jtab_fit_xmin, jtab_fit_dx, jtab_fit_invdx)

  end subroutine load_photoelectric_tables

  ! ***************************
  subroutine load_jtab(fname, nt, j_table, xmin, dx, invdx)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nt
    real*8,intent(inout)::j_table(zmin:zmax, nt), xmin, dx, invdx
    integer::z, unit, i, j
    real*8::tgas(nt)

    open(newunit=unit, file=trim(fname), status="old")
    do i=zmin, zmax
      do j=1,nt
        read(unit, *) z, tgas(j), j_table(i, j)
      end do
    end do
    close(unit)

    tgas = log10(tgas)
    j_table = log10(j_table + 1d-40)

    xmin = tgas(1)
    dx = (tgas(nt) - tgas(1)) / (nt - 1)
    invdx = 1d0 / dx

  end subroutine load_jtab

  ! ************************
  subroutine load_dust_cooling_table()
    implicit none
    integer::i, unit
    character(len=1024)::fname
    real*8:: dummy

    print *, "loading atomic tdust table..."

    fname = trim(runtime_data_folder)//"dust_cooling_grid.dat"
    call load_3d_fit(trim(fname), dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3, dust_cooling_table_data, &
      do_log=.true., skip_lines=1)

    fname = trim(runtime_data_folder)//"dust_heating_grid.dat"
    call load_3d_fit(trim(fname), dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3, dust_heating_table_data, &
      do_log=.true., skip_lines=1)

    fname = trim(runtime_data_folder)//"tdust_grid.dat"
    call load_3d_fit(trim(fname), dust_cooling_table_n1, dust_cooling_table_n2, dust_cooling_table_n3, tdust_table_data, &
      do_log=.true., skip_lines=1)

    ! call load_data_3d(trim(fname), dust_cooling_table, &
    !    dust_cooling_emiss, nemiss_dust_cooling_table, dust_cooling_emiss_min, dust_cooling_emiss_fact, dust_cooling_emiss_invdx, &
    !    dust_cooling_tgas, ntgas_dust_cooling_table, dust_cooling_tgas_min, dust_cooling_tgas_fact, dust_cooling_tgas_invdx, &
    !    dust_cooling_ngas, nngas_dust_cooling_table, dust_cooling_ngas_min, dust_cooling_ngas_fact, dust_cooling_ngas_invdx, &
    !    do_log=.true., skip_lines=1)

    open(newunit=unit, file=trim(runtime_data_folder)//"tdust_Ea_preint.dat", status="old")
    read(unit, *)
    do i=1,nphoto
      read(unit, *) dummy, pre_dust_cooling_table(i)
    end do
    close(unit)

  end subroutine load_dust_cooling_table

  ! ****************************
  subroutine load_dust_kappa_opacity()
    implicit none
    integer::i, unit

    print *, "loading dust opacity..."

    open(newunit=unit, file=trim(runtime_data_folder)//"kappa_dust.dat", status="old")
    do i=1,nphoto
      read(unit, *) dust_kappa_opacity(i)
    end do
    close(unit)

  end subroutine load_dust_kappa_opacity

  ! ***************************
  subroutine load_H2_cooling_tabs()
    implicit none
    integer::i, unit
    character(len=1024)::fnames(cool_H2_vec)

    print *, "loading H2 cooling..."

    fnames = (/trim(runtime_data_folder)//"cool_H2_H.dat" , &
      trim(runtime_data_folder)//"cool_H2_Hj.dat", &
      trim(runtime_data_folder)//"cool_H2_H2.dat", &
      trim(runtime_data_folder)//"cool_H2_e.dat", &
      trim(runtime_data_folder)//"cool_H2_HDL.dat"/)

      !load_1d_fit_vec(fname, nvec, nx, data, do_log, skip_lines)
    call load_1d_fit_vec(fnames, cool_H2_vec, cool_H2_nx, cooling_H2data, do_log=.true.)

    ! fname = trim(runtime_data_folder)//"cool_H2_H.dat"
    ! call load_data_1d(trim(fname), cooling_H2_H, &
    !   cool_H2_tgas, cool_H2_ntgas, cool_H2_min, cool_H2_fact, cool_H2_invdx, do_log=.true.)
    !
    ! fname = trim(runtime_data_folder)//"cool_H2_Hj.dat"
    ! call load_data_1d(trim(fname), cooling_H2_Hj, &
    !   cool_H2_tgas, cool_H2_ntgas, cool_H2_min, cool_H2_fact, cool_H2_invdx, do_log=.true.)
    !
    ! fname = trim(runtime_data_folder)//"cool_H2_H2.dat"
    ! call load_data_1d(trim(fname), cooling_H2_H2, &
    !   cool_H2_tgas, cool_H2_ntgas, cool_H2_min, cool_H2_fact, cool_H2_invdx, do_log=.true.)
    !
    ! fname = trim(runtime_data_folder)//"cool_H2_e.dat"
    ! call load_data_1d(trim(fname), cooling_H2_e, &
    !   cool_H2_tgas, cool_H2_ntgas, cool_H2_min, cool_H2_fact, cool_H2_invdx, do_log=.true.)
    !
    ! fname = trim(runtime_data_folder)//"cool_H2_HDL.dat"
    ! call load_data_1d(trim(fname), cooling_H2_HDL, &
    !   cool_H2_tgas, cool_H2_ntgas, cool_H2_min, cool_H2_fact, cool_H2_invdx, do_log=.true.)

  end subroutine load_H2_cooling_tabs

  ! ! ************************
  ! subroutine load_shielding_H2_table()
  !   implicit none
  !   integer::i, unit
  !   character(len=1024)::fname
  !   real*8:: dummy
  !
  !   print *, "loading shielding H2 table"
  !
  !   fname = trim(runtime_data_folder)//"shielding_H2.dat"
  !
  !   call load_data_2d(trim(fname), shielding_H2_table, &
  !      shielding_H2_ncol, nncol_shielding_H2, shielding_H2_ncol_min, shielding_H2_ncol_fact, shielding_H2_ncol_invdx, &
  !      shielding_H2_tgas, ntgas_shielding_H2, shielding_H2_tgas_min, shielding_H2_tgas_fact, shielding_H2_tgas_invdx, &
  !      do_log=.true.)
  !
  !   ! store to avoid numerical upper value problem, but physically OK
  !   shielding_H2_ncol_max = maxval(shielding_H2_ncol) * 0.99
  !
  ! end subroutine load_shielding_H2_table

  ! ************************
  subroutine load_shielding_H2_table()
    implicit none
    integer::i, unit
    character(len=1024)::fname
    real*8:: dummy

    print *, "loading shielding H2 table..."

    fname = trim(runtime_data_folder)//"shielding_H2.dat"

    call load_2d_fit(trim(fname), shielding_H2_n1, shielding_H2_n2, shielding_H2_data, do_log=.true.)

  end subroutine load_shielding_H2_table

  ! ************************
  subroutine load_shielding_CO_table()
    implicit none
    integer::i, unit
    character(len=1024)::fname
    real*8:: dummy

    print *, "loading shielding CO table..."

    fname = trim(runtime_data_folder)//"shielding_CO.dat"

    call load_2d_fit(trim(fname), shielding_CO_n1, shielding_CO_n2, shielding_CO_data, do_log=.true.)

  end subroutine load_shielding_CO_table

  ! ****************
  subroutine load_CO_cooling()
    implicit none
    character(len=1024)::fname

    print *, "loading CO cooling..."

    fname = trim(runtime_data_folder)//"CO_cooling.dat"

    call load_3d_fit(trim(fname), cool_CO_tab_n1, cool_CO_tab_n2, cool_CO_tab_n2, cool_CO_tab_data)

  end subroutine load_CO_cooling

  ! *****************
  subroutine load_verbatim_reactions()
    implicit none
    integer::i, unit

    print *, "loading verbatim reactions..."

    open(newunit=unit, file=trim(runtime_data_folder)//"reactions.dat", status="old")
    do i=1,nreactions
      read(unit, '(a100)') reactions_verbatim(i)
    end do
    close(unit)

  end subroutine load_verbatim_reactions

end module prizmo_loaders

module prizmo_emission

  !!BEGIN_EMISSION_COMMONS
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    real*8,parameter::kcoll_xmin=4.771213e-01
    real*8,parameter::kcoll_invdx=1.329278e+02
    real*8,parameter::kcoll_fact=1.327949e+02
    integer,parameter::nsteps=1000
    integer,parameter::nlevels_H=6
    real*8::Aij_H(nlevels_H, nlevels_H)
    real*8::kcoll_data_H_E(nlevels_H, nlevels_H, nsteps)
    integer,parameter::nlevels_C=3
    real*8::Aij_C(nlevels_C, nlevels_C)
    real*8::kcoll_data_C_H(nlevels_C, nlevels_C, nsteps)
    real*8::kcoll_data_C_Hj(nlevels_C, nlevels_C, nsteps)
    real*8::kcoll_data_C_E(nlevels_C, nlevels_C, nsteps)
    real*8::kcoll_data_C_oH2(nlevels_C, nlevels_C, nsteps)
    real*8::kcoll_data_C_pH2(nlevels_C, nlevels_C, nsteps)
    integer,parameter::nlevels_O=3
    real*8::Aij_O(nlevels_O, nlevels_O)
    real*8::kcoll_data_O_H(nlevels_O, nlevels_O, nsteps)
    real*8::kcoll_data_O_Hj(nlevels_O, nlevels_O, nsteps)
    real*8::kcoll_data_O_E(nlevels_O, nlevels_O, nsteps)
    integer,parameter::nlevels_CO=41
    real*8::Aij_CO(nlevels_CO, nlevels_CO)
    real*8::kcoll_data_CO_pH2(nlevels_CO, nlevels_CO, nsteps)
    real*8::kcoll_data_CO_oH2(nlevels_CO, nlevels_CO, nsteps)
    integer,parameter::nlevels_Hej=6
    real*8::Aij_Hej(nlevels_Hej, nlevels_Hej)
    real*8::kcoll_data_Hej_E(nlevels_Hej, nlevels_Hej, nsteps)
    integer,parameter::nlevels_He=6
    real*8::Aij_He(nlevels_He, nlevels_He)
    real*8::kcoll_data_He_E(nlevels_He, nlevels_He, nsteps)
    integer,parameter::nlevels_Cj=2
    real*8::Aij_Cj(nlevels_Cj, nlevels_Cj)
    real*8::kcoll_data_Cj_E(nlevels_Cj, nlevels_Cj, nsteps)
    real*8::kcoll_data_Cj_H(nlevels_Cj, nlevels_Cj, nsteps)
    integer,parameter::nlevels_Oj=3
    real*8::Aij_Oj(nlevels_Oj, nlevels_Oj)
    real*8::kcoll_data_Oj_E(nlevels_Oj, nlevels_Oj, nsteps)
    integer,parameter::ntransitions=70

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  !!END_EMISSION_COMMONS

  real*8,parameter::nlimit_frac = 1d-20
  real*8::xdata(nsteps)
  real*8::emission_array_energy(ntransitions)
  real*8::emission_array_flux(ntransitions)
  character(len=50)::emission_array_names(ntransitions)

contains

  ! **************************
  subroutine load_emission_data()
    implicit none
    !!BEGIN_EMISSION_LOAD
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    call load_data_Aij("runtime_data//rate_Aij_H.dat", nlevels_H, Aij_H)
    call load_data_kcoll("runtime_data//rate_collisions_H_E.dat", nlevels_H, kcoll_data_H_E)
    call load_data_Aij("runtime_data//rate_Aij_C.dat", nlevels_C, Aij_C)
    call load_data_kcoll("runtime_data//rate_collisions_C_H.dat", nlevels_C, kcoll_data_C_H)
    call load_data_kcoll("runtime_data//rate_collisions_C_Hj.dat", nlevels_C, kcoll_data_C_Hj)
    call load_data_kcoll("runtime_data//rate_collisions_C_E.dat", nlevels_C, kcoll_data_C_E)
    call load_data_kcoll("runtime_data//rate_collisions_C_oH2.dat", nlevels_C, kcoll_data_C_oH2)
    call load_data_kcoll("runtime_data//rate_collisions_C_pH2.dat", nlevels_C, kcoll_data_C_pH2)
    call load_data_Aij("runtime_data//rate_Aij_O.dat", nlevels_O, Aij_O)
    call load_data_kcoll("runtime_data//rate_collisions_O_H.dat", nlevels_O, kcoll_data_O_H)
    call load_data_kcoll("runtime_data//rate_collisions_O_Hj.dat", nlevels_O, kcoll_data_O_Hj)
    call load_data_kcoll("runtime_data//rate_collisions_O_E.dat", nlevels_O, kcoll_data_O_E)
    call load_data_Aij("runtime_data//rate_Aij_CO.dat", nlevels_CO, Aij_CO)
    call load_data_kcoll("runtime_data//rate_collisions_CO_pH2.dat", nlevels_CO, kcoll_data_CO_pH2)
    call load_data_kcoll("runtime_data//rate_collisions_CO_oH2.dat", nlevels_CO, kcoll_data_CO_oH2)
    call load_data_Aij("runtime_data//rate_Aij_Hej.dat", nlevels_Hej, Aij_Hej)
    call load_data_kcoll("runtime_data//rate_collisions_Hej_E.dat", nlevels_Hej, kcoll_data_Hej_E)
    call load_data_Aij("runtime_data//rate_Aij_He.dat", nlevels_He, Aij_He)
    call load_data_kcoll("runtime_data//rate_collisions_He_E.dat", nlevels_He, kcoll_data_He_E)
    call load_data_Aij("runtime_data//rate_Aij_Cj.dat", nlevels_Cj, Aij_Cj)
    call load_data_kcoll("runtime_data//rate_collisions_Cj_E.dat", nlevels_Cj, kcoll_data_Cj_E)
    call load_data_kcoll("runtime_data//rate_collisions_Cj_H.dat", nlevels_Cj, kcoll_data_Cj_H)
    call load_data_Aij("runtime_data//rate_Aij_Oj.dat", nlevels_Oj, Aij_Oj)
    call load_data_kcoll("runtime_data//rate_collisions_Oj_E.dat", nlevels_Oj, kcoll_data_Oj_E)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_EMISSION_LOAD

    call load_emission_names()

  end subroutine load_emission_data

  ! **************************
  ! load emission names as string
  subroutine load_emission_names()
    implicit none
    integer::i, unit

    open(newunit=unit, file="runtime_data/emission_names_list.dat", status="old")
    do i=1,ntransitions
      read(unit, '(a50)') emission_array_names(i)
    end do
    close(unit)

  end subroutine load_emission_names

  ! **************************
  ! load collisional rates from file, k(i, j, Tgas), cm-3/s
  subroutine load_data_kcoll(fname, nlevels, kdata)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nlevels
    integer::unit, j, ilev, jlev, ios
    real*8::ydata(nsteps)
    real*8,intent(out)::kdata(nlevels, nlevels, nsteps)

    ! i_level, j_level, Tgas_index
    ! large negative because fit is in log space
    kdata(:, :, :) = -1d99

    ! open file to read
    open(newunit=unit, file=trim(fname), status="old", iostat=ios)
    ! check if file exists
    if(ios /= 0) then
      print *, "ERROR: problems with "//trim(fname)
      stop
    end if

    ! loop to load, files are i_level, j_level blocks
    ! with nsteps temperature steps each
    do
      !read levels number
      read(unit, *, iostat=ios) ilev, jlev
      ! read until it can
      if (ios /= 0) then
        exit
      end if
      ! loop on temperature steps
      do j=1,nsteps
        read(unit, *) xdata(j), ydata(j)
      end do
      ! store log of data
      kdata(jlev, ilev, :) = log10(ydata(:) + 1d-40)
    end do
    close(unit)

    ! temperature grid is the same for all the collisional rates (commons)
    xdata(:) = log10(xdata)

  end subroutine load_data_kcoll

  ! ******************************
  ! load Aij coefficients, 1/s
  subroutine load_data_Aij(fname, nlevels, Aij)
    implicit none
    character(len=*),intent(in)::fname
    integer,intent(in)::nlevels
    integer::unit, ilev, jlev, ios
    real*8::aval
    real*8,intent(out)::Aij(nlevels, nlevels)

    ! i_level, j_level
    Aij(:, :) = 0d0

    ! read from file
    open(newunit=unit, file=fname, status="old")
    ! loop on transitions
    do
      read(unit, *, iostat=ios) ilev, jlev, aval
      ! read until it can
      if (ios /= 0) then
        exit
      end if
      ! transitions are stored into Aij=Aul matrix
      Aij(jlev, ilev) = aval
    end do
    close(unit)

  end subroutine load_data_Aij

  ! *********************
  ! get atomic cooling from emission, erg/s/cm3
  function get_atomic_cooling(n, Tgas_in) result(cool)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas_in
    real*8::n_pH2, n_oH2, logTgas, cool, Tgas, nlimit

    !!BEGIN_ATOMIC_COOLING_POPULATION_DECLARE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    real*8::Cij_H(nlevels_H, nlevels_H)
    real*8::npop_H(nlevels_H)
    real*8::Cij_C(nlevels_C, nlevels_C)
    real*8::npop_C(nlevels_C)
    real*8::Cij_O(nlevels_O, nlevels_O)
    real*8::npop_O(nlevels_O)
    real*8::Cij_Hej(nlevels_Hej, nlevels_Hej)
    real*8::npop_Hej(nlevels_Hej)
    real*8::Cij_He(nlevels_He, nlevels_He)
    real*8::npop_He(nlevels_He)
    real*8::Cij_Cj(nlevels_Cj, nlevels_Cj)
    real*8::npop_Cj(nlevels_Cj)
    real*8::Cij_Oj(nlevels_Oj, nlevels_Oj)
    real*8::npop_Oj(nlevels_Oj)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ATOMIC_COOLING_POPULATION_DECLARE

    ! density limit to compute cooling, cm-3
    nlimit = nlimit_frac * get_ntot(n(:))

    ! compute ortho/para ratio
    n_pH2 = n(idx_H2) / (opratio_H2 + 1d0)
    n_oH2 = n(idx_H2) - n_pH2

    ! check Tgas range
    Tgas = max(min(Tgas_in, Tgas_max), Tgas_min)

    ! log of Tgas for rate interpolation
    logTgas = log10(Tgas)

    !!BEGIN_ATOMIC_COOLING_POPULATION
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    if(n(idx_H) > nlimit) then
        Cij_H(:, :) = 0d0
        call add_Cij(Cij_H(:, :), kcoll_data_H_E(:, :, :), logTgas, nlevels_H, n(idx_E))
        npop_H(:) = population(Aij_H(:, :), Cij_H(:, :), nlevels_H, n(idx_H))
    end if

    if(n(idx_C) > nlimit) then
        Cij_C(:, :) = 0d0
        call add_Cij(Cij_C(:, :), kcoll_data_C_H(:, :, :), logTgas, nlevels_C, n(idx_H))
        call add_Cij(Cij_C(:, :), kcoll_data_C_Hj(:, :, :), logTgas, nlevels_C, n(idx_Hj))
        call add_Cij(Cij_C(:, :), kcoll_data_C_E(:, :, :), logTgas, nlevels_C, n(idx_E))
        call add_Cij(Cij_C(:, :), kcoll_data_C_oH2(:, :, :), logTgas, nlevels_C, n_oH2)
        call add_Cij(Cij_C(:, :), kcoll_data_C_pH2(:, :, :), logTgas, nlevels_C, n_pH2)
        npop_C(:) = population(Aij_C(:, :), Cij_C(:, :), nlevels_C, n(idx_C))
    end if

    if(n(idx_O) > nlimit) then
        Cij_O(:, :) = 0d0
        call add_Cij(Cij_O(:, :), kcoll_data_O_H(:, :, :), logTgas, nlevels_O, n(idx_H))
        call add_Cij(Cij_O(:, :), kcoll_data_O_Hj(:, :, :), logTgas, nlevels_O, n(idx_Hj))
        call add_Cij(Cij_O(:, :), kcoll_data_O_E(:, :, :), logTgas, nlevels_O, n(idx_E))
        npop_O(:) = population(Aij_O(:, :), Cij_O(:, :), nlevels_O, n(idx_O))
    end if

    if(n(idx_Hej) > nlimit) then
        Cij_Hej(:, :) = 0d0
        call add_Cij(Cij_Hej(:, :), kcoll_data_Hej_E(:, :, :), logTgas, nlevels_Hej, n(idx_E))
        npop_Hej(:) = population(Aij_Hej(:, :), Cij_Hej(:, :), nlevels_Hej, n(idx_Hej))
    end if

    if(n(idx_He) > nlimit) then
        Cij_He(:, :) = 0d0
        call add_Cij(Cij_He(:, :), kcoll_data_He_E(:, :, :), logTgas, nlevels_He, n(idx_E))
        npop_He(:) = population(Aij_He(:, :), Cij_He(:, :), nlevels_He, n(idx_He))
    end if

    if(n(idx_Cj) > nlimit) then
        Cij_Cj(:, :) = 0d0
        call add_Cij(Cij_Cj(:, :), kcoll_data_Cj_E(:, :, :), logTgas, nlevels_Cj, n(idx_E))
        call add_Cij(Cij_Cj(:, :), kcoll_data_Cj_H(:, :, :), logTgas, nlevels_Cj, n(idx_H))
        npop_Cj(:) = population(Aij_Cj(:, :), Cij_Cj(:, :), nlevels_Cj, n(idx_Cj))
    end if

    if(n(idx_Oj) > nlimit) then
        Cij_Oj(:, :) = 0d0
        call add_Cij(Cij_Oj(:, :), kcoll_data_Oj_E(:, :, :), logTgas, nlevels_Oj, n(idx_E))
        npop_Oj(:) = population(Aij_Oj(:, :), Cij_Oj(:, :), nlevels_Oj, n(idx_Oj))
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ATOMIC_COOLING_POPULATION

    cool = 0d0

    !!BEGIN_ATOMIC_COOLING_TOTAL
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    if(n(idx_H) > nlimit) then
        cool = cool + npop_H(2) * 1.634030e-11 * 2.496000e-06 &
            * beta_escape(2.496d-6, 1.63402955d-11, dvdz, npop_H(2), npop_H(1), 2d0, 2d0)
        cool = cool + npop_H(3) * 1.634029e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63402886d-11, dvdz, npop_H(3), npop_H(1), 2d0, 2d0)
        cool = cool + npop_H(4) * 1.634036e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63403613d-11, dvdz, npop_H(4), npop_H(1), 4d0, 2d0)
        cool = cool + npop_H(5) * 3.026013e-12 * 2.100000e+06 &
            * beta_escape(2.1d6, 3.02601303d-12, dvdz, npop_H(5), npop_H(3), 2d0, 2d0)
        cool = cool + npop_H(5) * 3.025940e-12 * 4.210000e+06 &
            * beta_escape(4.21d6, 3.02594033d-12, dvdz, npop_H(5), npop_H(4), 2d0, 4d0)
        cool = cool + npop_H(6) * 1.936630e-11 * 1.670000e+08 &
            * beta_escape(1.67d8, 1.93662994d-11, dvdz, npop_H(6), npop_H(1), 2d0, 2d0)
        cool = cool + npop_H(6) * 3.026004e-12 * 2.240000e+07 &
            * beta_escape(2.24d7, 3.02600389d-12, dvdz, npop_H(6), npop_H(2), 2d0, 2d0)
    end if

    if(n(idx_C) > nlimit) then
        cool = cool + npop_C(2) * 3.261092e-15 * 7.900000e-08 &
            * beta_escape(7.9d-8, 3.2610918d-15, dvdz, npop_C(2), npop_C(1), 3d0, 1d0)
        cool = cool + npop_C(3) * 8.623531e-15 * 2.100000e-14 &
            * beta_escape(2.1d-14, 8.62353066d-15, dvdz, npop_C(3), npop_C(1), 5d0, 1d0)
        cool = cool + npop_C(3) * 5.362439e-15 * 2.700000e-07 &
            * beta_escape(2.7d-7, 5.36243885d-15, dvdz, npop_C(3), npop_C(2), 5d0, 3d0)
    end if

    if(n(idx_O) > nlimit) then
        cool = cool + npop_O(2) * 3.143861e-14 * 8.900000e-05 &
            * beta_escape(8.9d-5, 3.14386094d-14, dvdz, npop_O(2), npop_O(1), 3d0, 5d0)
        cool = cool + npop_O(3) * 4.508784e-14 * 1.300000e-10 &
            * beta_escape(1.3d-10, 4.50878387d-14, dvdz, npop_O(3), npop_O(1), 1d0, 5d0)
        cool = cool + npop_O(3) * 1.364923e-14 * 1.800000e-05 &
            * beta_escape(1.8d-5, 1.36492293d-14, dvdz, npop_O(3), npop_O(2), 1d0, 3d0)
    end if

    if(n(idx_Hej) > nlimit) then
        cool = cool + npop_Hej(2) * 6.538978e-11 * 5.266000e+02 &
            * beta_escape(5.266d2, 6.53897756d-11, dvdz, npop_Hej(2), npop_Hej(1), 2d0, 2d0)
        cool = cool + npop_Hej(3) * 6.538968e-11 * 1.002000e+10 &
            * beta_escape(1.002d+10, 6.53896826d-11, dvdz, npop_Hej(3), npop_Hej(1), 2d0, 2d0)
        cool = cool + npop_Hej(4) * 6.539085e-11 * 1.003000e+10 &
            * beta_escape(1.003d+10, 6.53908461d-11, dvdz, npop_Hej(4), npop_Hej(1), 4d0, 2d0)
        cool = cool + npop_Hej(5) * 1.210971e-11 * 3.369000e+07 &
            * beta_escape(3.369d7, 1.21097056d-11, dvdz, npop_Hej(5), npop_Hej(3), 2d0, 2d0)
        cool = cool + npop_Hej(5) * 1.210854e-11 * 6.737000e+07 &
            * beta_escape(6.737d7, 1.21085421d-11, dvdz, npop_Hej(5), npop_Hej(4), 2d0, 4d0)
        cool = cool + npop_Hej(6) * 7.749936e-11 * 2.675000e+09 &
            * beta_escape(2.675d9, 7.74993604d-11, dvdz, npop_Hej(6), npop_Hej(1), 2d0, 2d0)
        cool = cool + npop_Hej(6) * 1.210958e-11 * 3.590000e+08 &
            * beta_escape(3.59d8, 1.21095848d-11, dvdz, npop_Hej(6), npop_Hej(2), 2d0, 2d0)
    end if

    if(n(idx_He) > nlimit) then
        cool = cool + npop_He(2) * 3.175454e-11 * 1.730000e-04 &
            * beta_escape(1.73d-4, 3.1754543d-11, dvdz, npop_He(2), npop_He(1), 3d0, 1d0)
        cool = cool + npop_He(3) * 3.303013e-11 * 5.094000e+01 &
            * beta_escape(5.094d1, 3.3030132d-11, dvdz, npop_He(3), npop_He(1), 1d0, 1d0)
        cool = cool + npop_He(4) * 1.833647e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83364677d-12, dvdz, npop_He(4), npop_He(2), 5d0, 3d0)
        cool = cool + npop_He(5) * 3.358820e-11 * 2.330000e+02 &
            * beta_escape(2.33d2, 3.35882049d-11, dvdz, npop_He(5), npop_He(1), 3d0, 1d0)
        cool = cool + npop_He(5) * 1.833662e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83366187d-12, dvdz, npop_He(5), npop_He(2), 3d0, 3d0)
        cool = cool + npop_He(6) * 1.833858e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83385813d-12, dvdz, npop_He(6), npop_He(2), 1d0, 3d0)
    end if

    if(n(idx_Cj) > nlimit) then
        cool = cool + npop_Cj(2) * 1.259842e-14 * 2.400000e-06 &
            * beta_escape(2.4d-6, 1.25984177d-14, dvdz, npop_Cj(2), npop_Cj(1), 4d0, 2d0)
    end if

    if(n(idx_Oj) > nlimit) then
        cool = cool + npop_Oj(2) * 5.325769e-12 * 5.100000e-05 &
            * beta_escape(5.1d-5, 5.32576883d-12, dvdz, npop_Oj(2), npop_Oj(1), 6d0, 1d0)
        cool = cool + npop_Oj(3) * 5.329745e-12 * 1.700000e-04 &
            * beta_escape(1.7d-4, 5.32974509d-12, dvdz, npop_Oj(3), npop_Oj(1), 4d0, 1d0)
        cool = cool + npop_Oj(3) * 3.976268e-15 * 1.300000e-07 &
            * beta_escape(1.3d-7, 3.97626774d-15, dvdz, npop_Oj(3), npop_Oj(2), 4d0, 6d0)
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_ATOMIC_COOLING_TOTAL

  end function get_atomic_cooling

  ! *********************
  ! add emission from given chemical species to the emission array
  subroutine add_emission(emission, n, Tgas_in)
    use prizmo_commons
    use prizmo_utils
    implicit none
    real*8,intent(in)::n(nmols), Tgas_in
    real*8,intent(inout)::emission(nphoto)
    real*8::Tgas, logTgas, n_pH2, n_oH2, emission_lost, nlimit
    !!BEGIN_EMISSION_POPULATION_DECLARE
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    real*8::Cij_H(nlevels_H, nlevels_H)
    real*8::npop_H(nlevels_H)
    real*8::Cij_C(nlevels_C, nlevels_C)
    real*8::npop_C(nlevels_C)
    real*8::Cij_O(nlevels_O, nlevels_O)
    real*8::npop_O(nlevels_O)
    real*8::Cij_CO(nlevels_CO, nlevels_CO)
    real*8::npop_CO(nlevels_CO)
    real*8::Cij_Hej(nlevels_Hej, nlevels_Hej)
    real*8::npop_Hej(nlevels_Hej)
    real*8::Cij_He(nlevels_He, nlevels_He)
    real*8::npop_He(nlevels_He)
    real*8::Cij_Cj(nlevels_Cj, nlevels_Cj)
    real*8::npop_Cj(nlevels_Cj)
    real*8::Cij_Oj(nlevels_Oj, nlevels_Oj)
    real*8::npop_Oj(nlevels_Oj)

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_EMISSION_POPULATION_DECLARE

    ! density limit to compute cooling, cm-3
    nlimit = nlimit_frac * get_ntot(n(:))

    ! compute ortho/para ratio
    n_pH2 = n(idx_H2) / (opratio_H2 + 1d0)
    n_oH2 = n(idx_H2) - n_pH2

    ! check Tgas range
    Tgas = max(min(Tgas_in, Tgas_max), Tgas_min)

    ! log of Tgas for rate interpolation
    logTgas = log10(Tgas)

    !!BEGIN_EMISSION_POPULATION
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    if(n(idx_H) > nlimit) then
        Cij_H(:, :) = 0d0
        call add_Cij(Cij_H(:, :), kcoll_data_H_E(:, :, :), logTgas, nlevels_H, n(idx_E))
        npop_H(:) = population(Aij_H(:, :), Cij_H(:, :), nlevels_H, n(idx_H))
    end if

    if(n(idx_C) > nlimit) then
        Cij_C(:, :) = 0d0
        call add_Cij(Cij_C(:, :), kcoll_data_C_H(:, :, :), logTgas, nlevels_C, n(idx_H))
        call add_Cij(Cij_C(:, :), kcoll_data_C_Hj(:, :, :), logTgas, nlevels_C, n(idx_Hj))
        call add_Cij(Cij_C(:, :), kcoll_data_C_E(:, :, :), logTgas, nlevels_C, n(idx_E))
        call add_Cij(Cij_C(:, :), kcoll_data_C_oH2(:, :, :), logTgas, nlevels_C, n_oH2)
        call add_Cij(Cij_C(:, :), kcoll_data_C_pH2(:, :, :), logTgas, nlevels_C, n_pH2)
        npop_C(:) = population(Aij_C(:, :), Cij_C(:, :), nlevels_C, n(idx_C))
    end if

    if(n(idx_O) > nlimit) then
        Cij_O(:, :) = 0d0
        call add_Cij(Cij_O(:, :), kcoll_data_O_H(:, :, :), logTgas, nlevels_O, n(idx_H))
        call add_Cij(Cij_O(:, :), kcoll_data_O_Hj(:, :, :), logTgas, nlevels_O, n(idx_Hj))
        call add_Cij(Cij_O(:, :), kcoll_data_O_E(:, :, :), logTgas, nlevels_O, n(idx_E))
        npop_O(:) = population(Aij_O(:, :), Cij_O(:, :), nlevels_O, n(idx_O))
    end if

    if(n(idx_CO) > nlimit) then
        Cij_CO(:, :) = 0d0
        call add_Cij(Cij_CO(:, :), kcoll_data_CO_pH2(:, :, :), logTgas, nlevels_CO, n_pH2)
        call add_Cij(Cij_CO(:, :), kcoll_data_CO_oH2(:, :, :), logTgas, nlevels_CO, n_oH2)
        npop_CO(:) = population(Aij_CO(:, :), Cij_CO(:, :), nlevels_CO, n(idx_CO))
    end if

    if(n(idx_Hej) > nlimit) then
        Cij_Hej(:, :) = 0d0
        call add_Cij(Cij_Hej(:, :), kcoll_data_Hej_E(:, :, :), logTgas, nlevels_Hej, n(idx_E))
        npop_Hej(:) = population(Aij_Hej(:, :), Cij_Hej(:, :), nlevels_Hej, n(idx_Hej))
    end if

    if(n(idx_He) > nlimit) then
        Cij_He(:, :) = 0d0
        call add_Cij(Cij_He(:, :), kcoll_data_He_E(:, :, :), logTgas, nlevels_He, n(idx_E))
        npop_He(:) = population(Aij_He(:, :), Cij_He(:, :), nlevels_He, n(idx_He))
    end if

    if(n(idx_Cj) > nlimit) then
        Cij_Cj(:, :) = 0d0
        call add_Cij(Cij_Cj(:, :), kcoll_data_Cj_E(:, :, :), logTgas, nlevels_Cj, n(idx_E))
        call add_Cij(Cij_Cj(:, :), kcoll_data_Cj_H(:, :, :), logTgas, nlevels_Cj, n(idx_H))
        npop_Cj(:) = population(Aij_Cj(:, :), Cij_Cj(:, :), nlevels_Cj, n(idx_Cj))
    end if

    if(n(idx_Oj) > nlimit) then
        Cij_Oj(:, :) = 0d0
        call add_Cij(Cij_Oj(:, :), kcoll_data_Oj_E(:, :, :), logTgas, nlevels_Oj, n(idx_E))
        npop_Oj(:) = population(Aij_Oj(:, :), Cij_Oj(:, :), nlevels_Oj, n(idx_Oj))
    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_EMISSION_POPULATION


    !!BEGIN_EMISSION_TOTAL
    ! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! NOTE: This block is auto-generated
    ! WHEN: 2020-12-15 15:59:15
    ! CHANGESET: xxxxxxx
    ! URL:
    ! BY: picogna@errc1

    if(n(idx_H) > nlimit) then
        ! wl = 1.215673e-01 micron, E = 1.019881e+01 eV
        emission_array_flux(1) = npop_H(2) * 1.634030e-11 * 2.496000e-06 &
            * beta_escape(2.496d-6, 1.63402955d-11, dvdz, npop_H(2), npop_H(1), 2d0, 2d0)
        emission_array_energy(1) = 1.019881e+01
        emission(634) = emission(634) + npop_H(2) * 1.634030e-11 * 2.496000e-06 &
            * beta_escape(2.496d-6, 1.63402955d-11, dvdz, npop_H(2), npop_H(1), 2d0, 2d0)

        ! wl = 1.215674e-01 micron, E = 1.019881e+01 eV
        emission_array_flux(2) = npop_H(3) * 1.634029e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63402886d-11, dvdz, npop_H(3), npop_H(1), 2d0, 2d0)
        emission_array_energy(2) = 1.019881e+01
        emission(634) = emission(634) + npop_H(3) * 1.634029e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63402886d-11, dvdz, npop_H(3), npop_H(1), 2d0, 2d0)

        ! wl = 1.215668e-01 micron, E = 1.019885e+01 eV
        emission_array_flux(3) = npop_H(4) * 1.634036e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63403613d-11, dvdz, npop_H(4), npop_H(1), 4d0, 2d0)
        emission_array_energy(3) = 1.019885e+01
        emission(634) = emission(634) + npop_H(4) * 1.634036e-11 * 6.260000e+08 &
            * beta_escape(6.26d8, 1.63403613d-11, dvdz, npop_H(4), npop_H(1), 4d0, 2d0)

        ! wl = 6.564564e-01 micron, E = 1.888689e+00 eV
        emission_array_flux(4) = npop_H(5) * 3.026013e-12 * 2.100000e+06 &
            * beta_escape(2.1d6, 3.02601303d-12, dvdz, npop_H(5), npop_H(3), 2d0, 2d0)
        emission_array_energy(4) = 1.888689e+00

        ! wl = 6.564722e-01 micron, E = 1.888643e+00 eV
        emission_array_flux(5) = npop_H(5) * 3.025940e-12 * 4.210000e+06 &
            * beta_escape(4.21d6, 3.02594033d-12, dvdz, npop_H(5), npop_H(4), 2d0, 4d0)
        emission_array_energy(5) = 1.888643e+00

        ! wl = 1.025723e-01 micron, E = 1.208749e+01 eV
        emission_array_flux(6) = npop_H(6) * 1.936630e-11 * 1.670000e+08 &
            * beta_escape(1.67d8, 1.93662994d-11, dvdz, npop_H(6), npop_H(1), 2d0, 2d0)
        emission_array_energy(6) = 1.208749e+01
        emission(851) = emission(851) + npop_H(6) * 1.936630e-11 * 1.670000e+08 &
            * beta_escape(1.67d8, 1.93662994d-11, dvdz, npop_H(6), npop_H(1), 2d0, 2d0)

        ! wl = 6.564584e-01 micron, E = 1.888683e+00 eV
        emission_array_flux(7) = npop_H(6) * 3.026004e-12 * 2.240000e+07 &
            * beta_escape(2.24d7, 3.02600389d-12, dvdz, npop_H(6), npop_H(2), 2d0, 2d0)
        emission_array_energy(7) = 1.888683e+00

    end if

    if(n(idx_C) > nlimit) then
        ! wl = 6.091352e+02 micron, E = 2.035413e-03 eV
        emission_array_flux(8) = npop_C(2) * 3.261092e-15 * 7.900000e-08 &
            * beta_escape(7.9d-8, 3.2610918d-15, dvdz, npop_C(2), npop_C(1), 3d0, 1d0)
        emission_array_energy(8) = 2.035413e-03

        ! wl = 2.303518e+02 micron, E = 5.382384e-03 eV
        emission_array_flux(9) = npop_C(3) * 8.623531e-15 * 2.100000e-14 &
            * beta_escape(2.1d-14, 8.62353066d-15, dvdz, npop_C(3), npop_C(1), 5d0, 1d0)
        emission_array_energy(9) = 5.382384e-03

        ! wl = 3.704370e+02 micron, E = 3.346971e-03 eV
        emission_array_flux(10) = npop_C(3) * 5.362439e-15 * 2.700000e-07 &
            * beta_escape(2.7d-7, 5.36243885d-15, dvdz, npop_C(3), npop_C(2), 5d0, 3d0)
        emission_array_energy(10) = 3.346971e-03

    end if

    if(n(idx_O) > nlimit) then
        ! wl = 6.318491e+01 micron, E = 1.962244e-02 eV
        emission_array_flux(11) = npop_O(2) * 3.143861e-14 * 8.900000e-05 &
            * beta_escape(8.9d-5, 3.14386094d-14, dvdz, npop_O(2), npop_O(1), 3d0, 5d0)
        emission_array_energy(11) = 1.962244e-02

        ! wl = 4.405724e+01 micron, E = 2.814162e-02 eV
        emission_array_flux(12) = npop_O(3) * 4.508784e-14 * 1.300000e-10 &
            * beta_escape(1.3d-10, 4.50878387d-14, dvdz, npop_O(3), npop_O(1), 1d0, 5d0)
        emission_array_energy(12) = 2.814162e-02

        ! wl = 1.455354e+02 micron, E = 8.519179e-03 eV
        emission_array_flux(13) = npop_O(3) * 1.364923e-14 * 1.800000e-05 &
            * beta_escape(1.8d-5, 1.36492293d-14, dvdz, npop_O(3), npop_O(2), 1d0, 3d0)
        emission_array_energy(13) = 8.519179e-03

    end if

    if(n(idx_CO) > nlimit) then
        ! wl = 2.600758e+03 micron, E = 4.767234e-04 eV
        emission_array_flux(14) = npop_CO(2) * 7.637950e-16 * 7.203000e-08 &
            * beta_escape(7.203d-8, 7.63795036d-16, dvdz, npop_CO(2), npop_CO(1), 3d0, 1d0)
        emission_array_energy(14) = 4.767234e-04

        ! wl = 1.300409e+03 micron, E = 9.534244e-04 eV
        emission_array_flux(15) = npop_CO(3) * 1.527554e-15 * 6.910000e-07 &
            * beta_escape(6.91d-7, 1.52755426d-15, dvdz, npop_CO(3), npop_CO(2), 5d0, 3d0)
        emission_array_energy(15) = 9.534244e-04

        ! wl = 8.669574e+02 micron, E = 1.430107e-03 eV
        emission_array_flux(16) = npop_CO(4) * 2.291284e-15 * 2.497000e-06 &
            * beta_escape(2.497d-6, 2.29128429d-15, dvdz, npop_CO(4), npop_CO(3), 7d0, 5d0)
        emission_array_energy(16) = 1.430107e-03

        ! wl = 6.502520e+02 micron, E = 1.906710e-03 eV
        emission_array_flux(17) = npop_CO(5) * 3.054886e-15 * 6.126000e-06 &
            * beta_escape(6.126d-6, 3.05488583d-15, dvdz, npop_CO(5), npop_CO(4), 9d0, 7d0)
        emission_array_energy(17) = 1.906710e-03

        ! wl = 5.202328e+02 micron, E = 2.383245e-03 eV
        emission_array_flux(18) = npop_CO(6) * 3.818379e-15 * 1.221000e-05 &
            * beta_escape(1.221d-5, 3.81837872d-15, dvdz, npop_CO(6), npop_CO(5), 1.1d1, 9d0)
        emission_array_energy(18) = 2.383245e-03

        ! wl = 4.335549e+02 micron, E = 2.859712e-03 eV
        emission_array_flux(19) = npop_CO(7) * 4.581763e-15 * 2.137000e-05 &
            * beta_escape(2.137d-5, 4.581763d-15, dvdz, npop_CO(7), npop_CO(6), 1.3d1, 1.1d1)
        emission_array_energy(19) = 2.859712e-03

        ! wl = 3.716512e+02 micron, E = 3.336036e-03 eV
        emission_array_flux(20) = npop_CO(8) * 5.344919e-15 * 3.422000e-05 &
            * beta_escape(3.422d-5, 5.34491948d-15, dvdz, npop_CO(8), npop_CO(7), 1.5d1, 1.3d1)
        emission_array_energy(20) = 3.336036e-03

        ! wl = 3.252252e+02 micron, E = 3.812256e-03 eV
        emission_array_flux(21) = npop_CO(9) * 6.107908e-15 * 5.134000e-05 &
            * beta_escape(5.134d-5, 6.10790763d-15, dvdz, npop_CO(9), npop_CO(8), 1.7d1, 1.5d1)
        emission_array_energy(21) = 3.812256e-03

        ! wl = 2.891197e+02 micron, E = 4.288334e-03 eV
        emission_array_flux(22) = npop_CO(10) * 6.870668e-15 * 7.330000e-05 &
            * beta_escape(7.33d-5, 6.8706682d-15, dvdz, npop_CO(10), npop_CO(9), 1.9d1, 1.7d1)
        emission_array_energy(22) = 4.288334e-03

        ! wl = 2.602403e+02 micron, E = 4.764220e-03 eV
        emission_array_flux(23) = npop_CO(11) * 7.633121e-15 * 1.006000e-04 &
            * beta_escape(1.006d-4, 7.6331213d-15, dvdz, npop_CO(11), npop_CO(10), 2.1d1, 1.9d1)
        emission_array_energy(23) = 4.764220e-03

        ! wl = 2.366133e+02 micron, E = 5.239951e-03 eV
        emission_array_flux(24) = npop_CO(12) * 8.395327e-15 * 1.339000e-04 &
            * beta_escape(1.339d-4, 8.3953268d-15, dvdz, npop_CO(12), npop_CO(11), 2.3d1, 2.1d1)
        emission_array_energy(24) = 5.239951e-03

        ! wl = 2.169271e+02 micron, E = 5.715478e-03 eV
        emission_array_flux(25) = npop_CO(13) * 9.157205e-15 * 1.735000e-04 &
            * beta_escape(1.735d-4, 9.15720524d-15, dvdz, npop_CO(13), npop_CO(12), 2.5d1, 2.3d1)
        emission_array_energy(25) = 5.715478e-03

        ! wl = 2.002725e+02 micron, E = 6.190776e-03 eV
        emission_array_flux(26) = npop_CO(14) * 9.918717e-15 * 2.200000e-04 &
            * beta_escape(2.2d-4, 9.91871683d-15, dvdz, npop_CO(14), npop_CO(13), 2.7d1, 2.5d1)
        emission_array_energy(26) = 6.190776e-03

        ! wl = 1.859995e+02 micron, E = 6.665833e-03 eV
        emission_array_flux(27) = npop_CO(15) * 1.067984e-14 * 2.739000e-04 &
            * beta_escape(2.739d-4, 1.06798418d-14, dvdz, npop_CO(15), npop_CO(14), 2.9d1, 2.7d1)
        emission_array_energy(27) = 6.665833e-03

        ! wl = 1.736313e+02 micron, E = 7.140661e-03 eV
        emission_array_flux(28) = npop_CO(16) * 1.144060e-14 * 3.354000e-04 &
            * beta_escape(3.354d-4, 1.14405998d-14, dvdz, npop_CO(16), npop_CO(15), 3.1d1, 2.9d1)
        emission_array_energy(28) = 7.140661e-03

        ! wl = 1.628115e+02 micron, E = 7.615198e-03 eV
        emission_array_flux(29) = npop_CO(17) * 1.220089e-14 * 4.050000e-04 &
            * beta_escape(4.05d-4, 1.22008919d-14, dvdz, npop_CO(17), npop_CO(16), 3.3d1, 3.1d1)
        emission_array_energy(29) = 7.615198e-03

        ! wl = 1.532669e+02 micron, E = 8.089431e-03 eV
        emission_array_flux(30) = npop_CO(18) * 1.296070e-14 * 4.829000e-04 &
            * beta_escape(4.829d-4, 1.2960698d-14, dvdz, npop_CO(18), npop_CO(17), 3.5d1, 3.3d1)
        emission_array_energy(30) = 8.089431e-03

        ! wl = 1.447841e+02 micron, E = 8.563386e-03 eV
        emission_array_flux(31) = npop_CO(19) * 1.372006e-14 * 5.695000e-04 &
            * beta_escape(5.695d-4, 1.37200577d-14, dvdz, npop_CO(19), npop_CO(18), 3.7d1, 3.5d1)
        emission_array_energy(31) = 8.563386e-03

        ! wl = 1.371964e+02 micron, E = 9.036989e-03 eV
        emission_array_flux(32) = npop_CO(20) * 1.447885e-14 * 6.650000e-04 &
            * beta_escape(6.65d-4, 1.44788521d-14, dvdz, npop_CO(20), npop_CO(19), 3.9d1, 3.7d1)
        emission_array_energy(32) = 9.036989e-03

        ! wl = 1.303690e+02 micron, E = 9.510250e-03 eV
        emission_array_flux(33) = npop_CO(21) * 1.523710e-14 * 7.695000e-04 &
            * beta_escape(7.695d-4, 1.5237101d-14, dvdz, npop_CO(21), npop_CO(20), 4.1d1, 3.9d1)
        emission_array_energy(33) = 9.510250e-03

        ! wl = 1.241933e+02 micron, E = 9.983159e-03 eV
        emission_array_flux(34) = npop_CO(22) * 1.599478e-14 * 8.833000e-04 &
            * beta_escape(8.833d-4, 1.59947848d-14, dvdz, npop_CO(22), npop_CO(21), 4.3d1, 4.1d1)
        emission_array_energy(34) = 9.983159e-03

        ! wl = 1.185807e+02 micron, E = 1.045568e-02 eV
        emission_array_flux(35) = npop_CO(23) * 1.675184e-14 * 1.006000e-03 &
            * beta_escape(1.006d-3, 1.67518432d-14, dvdz, npop_CO(23), npop_CO(22), 4.5d1, 4.3d1)
        emission_array_energy(35) = 1.045568e-02

        ! wl = 1.134575e+02 micron, E = 1.092781e-02 eV
        emission_array_flux(36) = npop_CO(24) * 1.750828e-14 * 1.139000e-03 &
            * beta_escape(1.139d-3, 1.75082769d-14, dvdz, npop_CO(24), npop_CO(23), 4.7d1, 4.5d1)
        emission_array_energy(36) = 1.092781e-02

        ! wl = 1.087629e+02 micron, E = 1.139950e-02 eV
        emission_array_flux(37) = npop_CO(25) * 1.826401e-14 * 1.281000e-03 &
            * beta_escape(1.281d-3, 1.82640062d-14, dvdz, npop_CO(25), npop_CO(24), 4.9d1, 4.7d1)
        emission_array_energy(37) = 1.139950e-02

        ! wl = 1.044449e+02 micron, E = 1.187077e-02 eV
        emission_array_flux(38) = npop_CO(26) * 1.901907e-14 * 1.432000e-03 &
            * beta_escape(1.432d-3, 1.90190708d-14, dvdz, npop_CO(26), npop_CO(25), 5.1d1, 4.9d1)
        emission_array_energy(38) = 1.187077e-02

        ! wl = 1.004605e+02 micron, E = 1.234158e-02 eV
        emission_array_flux(39) = npop_CO(27) * 1.977339e-14 * 1.592000e-03 &
            * beta_escape(1.592d-3, 1.97733915d-14, dvdz, npop_CO(27), npop_CO(26), 5.3d1, 5.1d1)
        emission_array_energy(39) = 1.234158e-02

        ! wl = 9.677249e+01 micron, E = 1.281193e-02 eV
        emission_array_flux(40) = npop_CO(28) * 2.052697e-14 * 1.761000e-03 &
            * beta_escape(1.761d-3, 2.05269681d-14, dvdz, npop_CO(28), npop_CO(27), 5.5d1, 5.3d1)
        emission_array_energy(40) = 1.281193e-02

        ! wl = 9.334906e+01 micron, E = 1.328178e-02 eV
        emission_array_flux(41) = npop_CO(29) * 2.127976e-14 * 1.940000e-03 &
            * beta_escape(1.94d-3, 2.12797609d-14, dvdz, npop_CO(29), npop_CO(28), 5.7d1, 5.5d1)
        emission_array_energy(41) = 1.328178e-02

        ! wl = 9.016303e+01 micron, E = 1.375111e-02 eV
        emission_array_flux(42) = npop_CO(30) * 2.203171e-14 * 2.126000e-03 &
            * beta_escape(2.126d-3, 2.20317104d-14, dvdz, npop_CO(30), npop_CO(29), 5.9d1, 5.7d1)
        emission_array_energy(42) = 1.375111e-02

        ! wl = 8.719045e+01 micron, E = 1.421993e-02 eV
        emission_array_flux(43) = npop_CO(31) * 2.278284e-14 * 2.321000e-03 &
            * beta_escape(2.321d-3, 2.27828366d-14, dvdz, npop_CO(31), npop_CO(30), 6.1d1, 5.9d1)
        emission_array_energy(43) = 1.421993e-02

        ! wl = 8.441072e+01 micron, E = 1.468821e-02 eV
        emission_array_flux(44) = npop_CO(32) * 2.353310e-14 * 2.524000e-03 &
            * beta_escape(2.524d-3, 2.35330995d-14, dvdz, npop_CO(32), npop_CO(31), 6.3d1, 6.1d1)
        emission_array_energy(44) = 1.468821e-02

        ! wl = 8.180579e+01 micron, E = 1.515592e-02 eV
        emission_array_flux(45) = npop_CO(33) * 2.428246e-14 * 2.735000e-03 &
            * beta_escape(2.735d-3, 2.42824595d-14, dvdz, npop_CO(33), npop_CO(32), 6.5d1, 6.3d1)
        emission_array_energy(45) = 1.515592e-02

        ! wl = 7.935982e+01 micron, E = 1.562304e-02 eV
        emission_array_flux(46) = npop_CO(34) * 2.503088e-14 * 2.952000e-03 &
            * beta_escape(2.952d-3, 2.5030877d-14, dvdz, npop_CO(34), npop_CO(33), 6.7d1, 6.5d1)
        emission_array_energy(46) = 1.562304e-02

        ! wl = 7.705868e+01 micron, E = 1.608958e-02 eV
        emission_array_flux(47) = npop_CO(35) * 2.577835e-14 * 3.175000e-03 &
            * beta_escape(3.175d-3, 2.57783523d-14, dvdz, npop_CO(35), npop_CO(34), 6.9d1, 6.7d1)
        emission_array_energy(47) = 1.608958e-02

        ! wl = 7.489006e+01 micron, E = 1.655549e-02 eV
        emission_array_flux(48) = npop_CO(36) * 2.652483e-14 * 3.404000e-03 &
            * beta_escape(3.404d-3, 2.65248254d-14, dvdz, npop_CO(36), npop_CO(35), 7.1d1, 6.9d1)
        emission_array_energy(48) = 1.655549e-02

        ! wl = 7.284289e+01 micron, E = 1.702077e-02 eV
        emission_array_flux(49) = npop_CO(37) * 2.727028e-14 * 3.638000e-03 &
            * beta_escape(3.638d-3, 2.72702764d-14, dvdz, npop_CO(37), npop_CO(36), 7.3d1, 7.1d1)
        emission_array_energy(49) = 1.702077e-02

        ! wl = 7.090720e+01 micron, E = 1.748542e-02 eV
        emission_array_flux(50) = npop_CO(38) * 2.801473e-14 * 3.878000e-03 &
            * beta_escape(3.878d-3, 2.80147253d-14, dvdz, npop_CO(38), npop_CO(37), 7.5d1, 7.3d1)
        emission_array_energy(50) = 1.748542e-02

        ! wl = 6.907442e+01 micron, E = 1.794937e-02 eV
        emission_array_flux(51) = npop_CO(39) * 2.875805e-14 * 4.120000e-03 &
            * beta_escape(4.12d-3, 2.87580535d-14, dvdz, npop_CO(39), npop_CO(38), 7.7d1, 7.5d1)
        emission_array_energy(51) = 1.794937e-02

        ! wl = 6.733646e+01 micron, E = 1.841264e-02 eV
        emission_array_flux(52) = npop_CO(40) * 2.950030e-14 * 4.365000e-03 &
            * beta_escape(4.365d-3, 2.95002993d-14, dvdz, npop_CO(40), npop_CO(39), 7.9d1, 7.7d1)
        emission_array_energy(52) = 1.841264e-02

        ! wl = 6.568625e+01 micron, E = 1.887521e-02 eV
        emission_array_flux(53) = npop_CO(41) * 3.024142e-14 * 4.613000e-03 &
            * beta_escape(4.613d-3, 3.0241425d-14, dvdz, npop_CO(41), npop_CO(40), 8.1d1, 7.9d1)
        emission_array_energy(53) = 1.887521e-02

    end if

    if(n(idx_Hej) > nlimit) then
        ! wl = 3.037854e-02 micron, E = 4.081309e+01 eV
        emission_array_flux(54) = npop_Hej(2) * 6.538978e-11 * 5.266000e+02 &
            * beta_escape(5.266d2, 6.53897756d-11, dvdz, npop_Hej(2), npop_Hej(1), 2d0, 2d0)
        emission_array_energy(54) = 4.081309e+01

        ! wl = 3.037858e-02 micron, E = 4.081303e+01 eV
        emission_array_flux(55) = npop_Hej(3) * 6.538968e-11 * 1.002000e+10 &
            * beta_escape(1.002d+10, 6.53896826d-11, dvdz, npop_Hej(3), npop_Hej(1), 2d0, 2d0)
        emission_array_energy(55) = 4.081303e+01

        ! wl = 3.037804e-02 micron, E = 4.081376e+01 eV
        emission_array_flux(56) = npop_Hej(4) * 6.539085e-11 * 1.003000e+10 &
            * beta_escape(1.003d+10, 6.53908461d-11, dvdz, npop_Hej(4), npop_Hej(1), 4d0, 2d0)
        emission_array_energy(56) = 4.081376e+01

        ! wl = 1.640375e-01 micron, E = 7.558284e+00 eV
        emission_array_flux(57) = npop_Hej(5) * 1.210971e-11 * 3.369000e+07 &
            * beta_escape(3.369d7, 1.21097056d-11, dvdz, npop_Hej(5), npop_Hej(3), 2d0, 2d0)
        emission_array_energy(57) = 7.558284e+00
        emission(277) = emission(277) + npop_Hej(5) * 1.210971e-11 * 3.369000e+07 &
            * beta_escape(3.369d7, 1.21097056d-11, dvdz, npop_Hej(5), npop_Hej(3), 2d0, 2d0)

        ! wl = 1.640533e-01 micron, E = 7.557558e+00 eV
        emission_array_flux(58) = npop_Hej(5) * 1.210854e-11 * 6.737000e+07 &
            * beta_escape(6.737d7, 1.21085421d-11, dvdz, npop_Hej(5), npop_Hej(4), 2d0, 4d0)
        emission_array_energy(58) = 7.557558e+00
        emission(277) = emission(277) + npop_Hej(5) * 1.210854e-11 * 6.737000e+07 &
            * beta_escape(6.737d7, 1.21085421d-11, dvdz, npop_Hej(5), npop_Hej(4), 2d0, 4d0)

        ! wl = 2.563177e-02 micron, E = 4.837130e+01 eV
        emission_array_flux(59) = npop_Hej(6) * 7.749936e-11 * 2.675000e+09 &
            * beta_escape(2.675d9, 7.74993604d-11, dvdz, npop_Hej(6), npop_Hej(1), 2d0, 2d0)
        emission_array_energy(59) = 4.837130e+01

        ! wl = 1.640391e-01 micron, E = 7.558208e+00 eV
        emission_array_flux(60) = npop_Hej(6) * 1.210958e-11 * 3.590000e+08 &
            * beta_escape(3.59d8, 1.21095848d-11, dvdz, npop_Hej(6), npop_Hej(2), 2d0, 2d0)
        emission_array_energy(60) = 7.558208e+00
        emission(277) = emission(277) + npop_Hej(6) * 1.210958e-11 * 3.590000e+08 &
            * beta_escape(3.59d8, 1.21095848d-11, dvdz, npop_Hej(6), npop_Hej(2), 2d0, 2d0)

    end if

    if(n(idx_He) > nlimit) then
        ! wl = 6.255627e-02 micron, E = 1.981963e+01 eV
        emission_array_flux(61) = npop_He(2) * 3.175454e-11 * 1.730000e-04 &
            * beta_escape(1.73d-4, 3.1754543d-11, dvdz, npop_He(2), npop_He(1), 3d0, 1d0)
        emission_array_energy(61) = 1.981963e+01

        ! wl = 6.014041e-02 micron, E = 2.061579e+01 eV
        emission_array_flux(62) = npop_He(3) * 3.303013e-11 * 5.094000e+01 &
            * beta_escape(5.094d1, 3.3030132d-11, dvdz, npop_He(3), npop_He(1), 1d0, 1d0)
        emission_array_energy(62) = 2.061579e+01

        ! wl = 1.083331e+00 micron, E = 1.144472e+00 eV
        emission_array_flux(63) = npop_He(4) * 1.833647e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83364677d-12, dvdz, npop_He(4), npop_He(2), 5d0, 3d0)
        emission_array_energy(63) = 1.144472e+00

        ! wl = 5.914117e-02 micron, E = 2.096411e+01 eV
        emission_array_flux(64) = npop_He(5) * 3.358820e-11 * 2.330000e+02 &
            * beta_escape(2.33d2, 3.35882049d-11, dvdz, npop_He(5), npop_He(1), 3d0, 1d0)
        emission_array_energy(64) = 2.096411e+01

        ! wl = 1.083322e+00 micron, E = 1.144482e+00 eV
        emission_array_flux(65) = npop_He(5) * 1.833662e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83366187d-12, dvdz, npop_He(5), npop_He(2), 3d0, 3d0)
        emission_array_energy(65) = 1.144482e+00

        ! wl = 1.083206e+00 micron, E = 1.144604e+00 eV
        emission_array_flux(66) = npop_He(6) * 1.833858e-12 * 1.150000e+07 &
            * beta_escape(1.15d7, 1.83385813d-12, dvdz, npop_He(6), npop_He(2), 1d0, 3d0)
        emission_array_energy(66) = 1.144604e+00

    end if

    if(n(idx_Cj) > nlimit) then
        ! wl = 1.576742e+02 micron, E = 7.863314e-03 eV
        emission_array_flux(67) = npop_Cj(2) * 1.259842e-14 * 2.400000e-06 &
            * beta_escape(2.4d-6, 1.25984177d-14, dvdz, npop_Cj(2), npop_Cj(1), 4d0, 2d0)
        emission_array_energy(67) = 7.863314e-03

    end if

    if(n(idx_Oj) > nlimit) then
        ! wl = 3.729876e-01 micron, E = 3.324083e+00 eV
        emission_array_flux(68) = npop_Oj(2) * 5.325769e-12 * 5.100000e-05 &
            * beta_escape(5.1d-5, 5.32576883d-12, dvdz, npop_Oj(2), npop_Oj(1), 6d0, 1d0)
        emission_array_energy(68) = 3.324083e+00

        ! wl = 3.727093e-01 micron, E = 3.326565e+00 eV
        emission_array_flux(69) = npop_Oj(3) * 5.329745e-12 * 1.700000e-04 &
            * beta_escape(1.7d-4, 5.32974509d-12, dvdz, npop_Oj(3), npop_Oj(1), 4d0, 1d0)
        emission_array_energy(69) = 3.326565e+00

        ! wl = 4.995755e+02 micron, E = 2.481791e-03 eV
        emission_array_flux(70) = npop_Oj(3) * 3.976268e-15 * 1.300000e-07 &
            * beta_escape(1.3d-7, 3.97626774d-15, dvdz, npop_Oj(3), npop_Oj(2), 4d0, 6d0)
        emission_array_energy(70) = 2.481791e-03

    end if

    ! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    !!END_EMISSION_TOTAL

  end subroutine add_emission

  ! **************
  ! add collisional rates to the Cij matrix for a given collider
  ! using linear interpolation in temperature
  ! kcoll, matrix contating the fit data for each transition
  ! logTgas, log10(temperature), log10(K)
  ! nlevels, number of levels
  ! ncollider, abundance of the collider, cm-3
  subroutine add_Cij(Cij, kcoll, logTgas, nlevels, ncollider)
    implicit none
    integer,intent(in)::nlevels
    real*8,intent(in)::logTgas, kcoll(nlevels, nlevels, nsteps), ncollider
    real*8::k(nlevels, nlevels), Cij(nlevels, nlevels), pre
    integer::idx

    ! find index
    idx = floor((logTgas - kcoll_xmin) * kcoll_fact) + 1

    if(idx < 1) then
      !print *, idx, logTgas, kcoll_xmin
      idx = max(idx, 1)
    end if

    ! precompute fit factor
    pre = (logTgas - xdata(idx)) * kcoll_invdx
    ! do fit
    k(:, :) = 1d1**(pre * (kcoll(:, :, idx+1) - kcoll(:, :, idx)) + kcoll(:, :, idx))

    ! add rate to the matrix
    Cij(:, :) = Cij(:, :) + k(:, :) * ncollider

  end subroutine add_Cij

  ! *********************
  ! compute level population with the given atomic data
  ! Aij, matrix of dexctiation coefficients, 1/s
  ! Cij, matrix of (de)exctiation collisiona rate coefficients, cm-3/s
  ! nlevels, number of energy levels
  ! ncool, coolant number density, cm-3
  ! returns: levels array of population number density, cm-3
  ! see https://www.aanda.org/articles/aa/pdf/2009/25/aa11821-09.pdf
  function population(Aij, Cij, nlevels, ncool) result(n)
    use prizmo_commons
    use prizmo_linear_solver
    implicit none
    integer,intent(in)::nlevels
    real*8,intent(in)::Aij(nlevels, nlevels), ncool
    real*8::Cij(nlevels, nlevels), Rij(nlevels, nlevels), pre, n(nlevels), logT
    integer::idx, i, j, ierr, ipiv(nlevels)

    ! Rij is a matrix where indexes are levels (up, low)
    Rij(:, :) = Aij(:, :) + Cij(:, :)

    ! loop on the diagonal to construct the negative terms of the matrix.
    ! summing Rij(i, i) is for consistency, but it should be already zero
    ! since Aii and Cii are both not transitions (hence they are both zero)
    do i=1,nlevels
      Rij(i, i) = -sum(Rij(:, i)) + Rij(i, i)
    end do

    ! first row is set to one for mass conservation
    Rij(1, :) = 1d0
    ! RHS of the system Rij * x = n
    n(1) = ncool
    n(2:nlevels) = 0d0

    ! use analytical linear solver for system of size 2 and 3
    ! dgesv otherwise
    n(:) = linear_solver(Rij(: ,:), n(:))

    ! avoid negative levels
    do i=1,nlevels
      n(i) = max(n(i), 0d0)
    end do

    ! normalize to ensure mass conservation
    n(:) = n(:) / (sum(n) + 1d-40) * ncool

  end function population

  ! ***********************
  ! line optical thickness (escape probability)
  ! Tielens+2005, pag.50, eqns. 2.43-2.45
  ! input arguments are cgs
  ! Aul, Einstein coefficient
  ! Eul, energy difference
  ! dvdz_in, velocity gradient
  ! n_up, n_low, population of the levels
  ! g_up, g_low, statistical weight
  function beta_escape(Aul, Eul, dvdz_in, n_up, n_low, g_up, g_low) result(beta)
    use prizmo_commons
    implicit none
    real*8,intent(in)::Aul, Eul, dvdz_in, n_up, n_low, g_up, g_low
    real*8,parameter::pre=(clight * hplanck)**3 / 8e0 / pi
    real*8::tau_ul, beta

    if(beta_escape_mode == -1) then
      beta = 1d0
      return
    end if

    tau_ul = pre * Aul / Eul**3 * n_up / dvdz_in * (n_low * g_up / (n_up + 1d-40) / g_low - 1d0)
    ! this avoids numerical instability
    tau_ul = max(tau_ul, 1d-10)

    if(tau_ul > 7d0) then
       beta = 1d0 / (4d0 * tau_ul * sqrt(log(tau_ul / sqrt(pi))))
    else
       beta = (1d0 - exp(-2.34 * tau_ul)) / 4.68d0 / tau_ul
    end if

    ! check limits
    if(beta < 0d0 .or. beta > 1d0) then
       print *, "ERROR: problem with beta escape!"
       print *, " beta < 0 or beta > 1:", beta
       stop
    end if

  end function beta_escape


end module prizmo_emission

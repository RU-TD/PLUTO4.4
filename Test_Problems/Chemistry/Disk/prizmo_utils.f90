module prizmo_utils
  use prizmo_commons
contains

  ! ************************
  function get_electrons(x) result(ne)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::ne

    !! PREPROCESS_ELECTRONS

    ne = x(idx_CH2j) + x(idx_CH3j) + x(idx_CH4j) + x(idx_CH5j) + x(idx_CHj) + x(idx_COj) + x(idx_Cj) + 2 * x(idx_Cjj) + 3 * x(idx_Cjjj) + 4 * x(idx_Cjjjj) + x(idx_H2Oj) + x(idx_H2j) + x(idx_H3Oj) + x(idx_H3j) + x(idx_HCOj) + x(idx_Hej) + 2 * x(idx_Hejj) + x(idx_Hj) + x(idx_O2j) + x(idx_OHj) + x(idx_Oj) + 2 * x(idx_Ojj) + 3 * x(idx_Ojjj) + 4 * x(idx_Ojjjj)

    !! PREPROCESS_END

    if(ne < 0d0) then
      print *, "ERROR: negative electrons!"
      stop
    end if

    ne = max(ne, 0d0)

  end function get_electrons

  ! ************************
  function get_rho(x) result(rho)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::rho

    !! PREPROCESS_GET_RHO

    rho = x(idx_C) * masses(idx_C)&
     + x(idx_CH) * masses(idx_CH)&
     + x(idx_CH2) * masses(idx_CH2)&
     + x(idx_CH2j) * masses(idx_CH2j)&
     + x(idx_CH3) * masses(idx_CH3)&
     + x(idx_CH3j) * masses(idx_CH3j)&
     + x(idx_CH4) * masses(idx_CH4)&
     + x(idx_CH4j) * masses(idx_CH4j)&
     + x(idx_CH5j) * masses(idx_CH5j)&
     + x(idx_CHj) * masses(idx_CHj)&
     + x(idx_CO) * masses(idx_CO)&
     + x(idx_COj) * masses(idx_COj)&
     + x(idx_Cj) * masses(idx_Cj)&
     + x(idx_Cjj) * masses(idx_Cjj)&
     + x(idx_Cjjj) * masses(idx_Cjjj)&
     + x(idx_Cjjjj) * masses(idx_Cjjjj)&
     + x(idx_E) * masses(idx_E)&
     + x(idx_H) * masses(idx_H)&
     + x(idx_H2) * masses(idx_H2)&
     + x(idx_H2O) * masses(idx_H2O)&
     + x(idx_H2Oj) * masses(idx_H2Oj)&
     + x(idx_H2j) * masses(idx_H2j)&
     + x(idx_H3Oj) * masses(idx_H3Oj)&
     + x(idx_H3j) * masses(idx_H3j)&
     + x(idx_HCOj) * masses(idx_HCOj)&
     + x(idx_He) * masses(idx_He)&
     + x(idx_Hej) * masses(idx_Hej)&
     + x(idx_Hejj) * masses(idx_Hejj)&
     + x(idx_Hj) * masses(idx_Hj)&
     + x(idx_O) * masses(idx_O)&
     + x(idx_O2) * masses(idx_O2)&
     + x(idx_O2j) * masses(idx_O2j)&
     + x(idx_OH) * masses(idx_OH)&
     + x(idx_OHj) * masses(idx_OHj)&
     + x(idx_Oj) * masses(idx_Oj)&
     + x(idx_Ojj) * masses(idx_Ojj)&
     + x(idx_Ojjj) * masses(idx_Ojjj)&
     + x(idx_Ojjjj) * masses(idx_Ojjjj)

    !! PREPROCESS_END

  end function get_rho

  ! ************************
  function get_Hnuclei(x) result(nH)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nH

    !! PREPROCESS_HNUCLEI

    nH = x(idx_CH) + 2 * x(idx_CH2) + 2 * x(idx_CH2j) + 3 * x(idx_CH3) + 3 * x(idx_CH3j) + 4 * x(idx_CH4) + 4 * x(idx_CH4j) + 5 * x(idx_CH5j) + x(idx_CHj) + x(idx_H) + 2 * x(idx_H2) + 2 * x(idx_H2O) + 2 * x(idx_H2O_DUST) + 2 * x(idx_H2Oj) + 2 * x(idx_H2j) + 3 * x(idx_H3Oj) + 3 * x(idx_H3j) + x(idx_HCOj) + x(idx_Hj) + x(idx_OH) + x(idx_OHj)

    !! PREPROCESS_END

    if(nH < 0d0) then
      print *, "ERROR: negative H nuclei!"
      print *, nH
      !stop
    end if

    nH = max(nH, 0d0)

  end function get_Hnuclei

  ! ************************
  function get_Cnuclei(x) result(nC)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nC

    !! PREPROCESS_CNUCLEI

    nC = x(idx_C) + x(idx_CH) + x(idx_CH2) + x(idx_CH2j) + x(idx_CH3) + x(idx_CH3j) + x(idx_CH4) + x(idx_CH4j) + x(idx_CH5j) + x(idx_CHj) + x(idx_CO) + x(idx_CO_DUST) + x(idx_COj) + x(idx_Cj) + x(idx_Cjj) + x(idx_Cjjj) + x(idx_Cjjjj) + x(idx_HCOj)

    !! PREPROCESS_END

    if(nC < 0d0) then
      print *, "ERROR: negative C nuclei!"
      print *, nC
      !stop
    end if

    nC = max(nC, 0d0)

  end function get_Cnuclei

  ! ************************
  function get_Onuclei(x) result(nO)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nO

    !! PREPROCESS_ONUCLEI

    nO = x(idx_CO) + x(idx_CO_DUST) + x(idx_COj) + x(idx_H2O) + x(idx_H2O_DUST) + x(idx_H2Oj) + x(idx_H3Oj) + x(idx_HCOj) + x(idx_O) + 2 * x(idx_O2) + 2 * x(idx_O2j) + x(idx_OH) + x(idx_OHj) + x(idx_Oj) + x(idx_Ojj) + x(idx_Ojjj) + x(idx_Ojjjj)

    !! PREPROCESS_END

    if(nO < 0d0) then
      print *, "ERROR: negative O nuclei!"
      print *, nO
      !stop
    end if

    nO = max(nO, 0d0)

  end function get_Onuclei

  ! ************************
  function get_HEnuclei(x) result(nHE)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::nHE

    !! PREPROCESS_HENUCLEI

    nHe = x(idx_He) + x(idx_Hej) + x(idx_Hejj)

    !! PREPROCESS_END

    if(nHE < 0d0) then
      print *, "ERROR: negative HE nuclei!"
      print *, nHE
      !stop
    end if

    nHE = max(nHE, 0d0)

  end function get_HEnuclei

  ! **************************
  subroutine ranker(array, nbest, array_string)
    implicit none
    real*8,intent(in):: array(:)
    integer,intent(in)::nbest
    character(len=*),intent(in):: array_string(:)
    integer::idxs(size(array)), i
    real*8::amax

    idxs = argsort_r(array)

    amax = array(idxs(1))

    print *, "**********************"
    do i=1,min(nbest, size(array))
      print '(2I5,2E17.8e3,x,a30)', i, idxs(i), array(idxs(i)), array(idxs(i)) / (amax + 1d-40), trim(array_string(idxs(i)))
    end do

  end subroutine ranker

  ! **********************
  function argsort_r(array) result(idxs)
    implicit none
    real*8,intent(in)::array(:)
    integer::idxs(size(array)), na

    na = size(array)

    idxs = argsort(array)
    idxs = idxs(na:1:-1)

  end function

  ! **********************
  function argsort(array) result(idxs)
    implicit none
    real*8,intent(in)::array(:)
    real*8::a(size(array)), rtmp
    integer::idxs(size(array)), i, na, itmp
    logical::swap

    na = size(array)
    do i=1,na
      idxs(i) = i
    end do

    a = array

    do
      swap = .false.
      do i=2,na
        if(a(i-1) > a(i)) then
          rtmp = a(i)
          a(i) = a(i-1)
          a(i-1) = rtmp
          itmp = idxs(i)
          idxs(i) = idxs(i-1)
          idxs(i-1) = itmp
          swap = .true.
        end if
      end do
      if(.not.swap) exit
    end do

  end function argsort

end module prizmo_utils

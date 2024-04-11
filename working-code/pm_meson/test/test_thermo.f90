module test_thermo
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds , only : dp
    !> imports go here
    implicit none
    private
  
    public :: collect_thermo
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_thermo(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_calc_lnq_sys", test_calc_lnq_sys), &
      new_unittest("test_add_lnq_indi", test_add_lnq_indi), &
      new_unittest("test_calc_helmholtz_energy", test_calc_helmholtz_energy), &
      new_unittest("test_calc_gibbs_enthalpy", test_calc_gibbs_enthalpy), &
      new_unittest("test_calc_internal_energy", test_calc_internal_energy) &
      ]
  
  end subroutine collect_thermo

  ! --------------------------------
  ! Unit test for calculate_lnq_sys
  ! --------------------------------
  ! Calculation of the system partition function  
  subroutine test_calc_lnq_sys(error)
  
    ! Dependencies
    use thermo, only : calculate_lnq_sys

    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: ntemp = 5
    integer :: nclust = 4
    real(dp), dimension(5, 4) :: pop = reshape([0.5_dp, 0.4_dp, 0.3_dp, 0.2_dp, 0.1_dp, &
                                                0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp, &
                                                0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
                                                0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], [5, 4])
    real(dp), dimension(5, 4) :: ln_clust = reshape([-5.2_dp, -4.2_dp, -3.2_dp, -2.2_dp, -1.2_dp, &
                                                      3.8_dp,  3.3_dp,  3.3_dp,  3.3_dp,  3.3_dp, &
                                                      2.3_dp,  2.8_dp,  3.3_dp,  4.3_dp,  5.3_dp, &
                                                      1.3_dp,  1.3_dp,  1.3_dp,  1.3_dp,  1.3_dp], [5, 4])
    real(dp), dimension(5) :: lnq_sys

    ! Expected values
    real(dp) :: expected(5) = [-1.10000000e+00_dp, -3.33066907e-16_dp, &
                                1.15000000e+00_dp,  2.40000000e+00_dp,  3.65000000e+00_dp]
    
    ! Call the subroutine
    call calculate_lnq_sys(ntemp, nclust, pop, ln_clust, lnq_sys)

    ! Check the results
    call check(error, lnq_sys(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq_sys(2), expected(2), thr=thr, rel=.false.)
    call check(error, lnq_sys(3), expected(3), thr=thr, rel=.false.)
    call check(error, lnq_sys(4), expected(4), thr=thr, rel=.false.)
    call check(error, lnq_sys(5), expected(5), thr=thr, rel=.false.)
  
  end subroutine test_calc_lnq_sys


  ! --------------------------------  
  ! Unit test for add_lnq_indi
  ! --------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_add_lnq_indi(error)
  
    ! Dependencies
    use thermo, only : add_lnq_indi

    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: ntemp = 5
    integer :: nclust = 4
    real(dp), dimension(5, 4) :: pop = reshape([0.5_dp, 0.4_dp, 0.3_dp, 0.2_dp, 0.1_dp, &
                                                0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp, &
                                                0.1_dp, 0.2_dp, 0.3_dp, 0.4_dp, 0.5_dp, &
                                                0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp, 0.1_dp], [5, 4])
    real(dp), dimension(5) :: lnq_sys = [-1.10000000e+00_dp, -3.33066907e-16_dp, &
                                         1.15000000e+00_dp,  2.40000000e+00_dp,  3.65000000e+00_dp]

    ! Expected values
    real(dp) :: expected(5) = [-0.771298070337215_dp, 0.363034254943387_dp, 1.52439686978342_dp, &
                                2.76303425494339_dp, 3.97870192966279_dp]
    
    ! Call the subroutine
    call add_lnq_indi(ntemp, nclust, pop, lnq_sys)

    ! Check the results
    call check(error, lnq_sys(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq_sys(2), expected(2), thr=thr, rel=.false.)
    call check(error, lnq_sys(3), expected(3), thr=thr, rel=.false.)
    call check(error, lnq_sys(4), expected(4), thr=thr, rel=.false.)
    call check(error, lnq_sys(5), expected(5), thr=thr, rel=.false.)
  
  end subroutine test_add_lnq_indi

  ! -----------------------------------------
  ! Unit test for calculate_helmholtz_energy
  ! -----------------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_calc_helmholtz_energy(error)
  
    ! Dependencies
    use thermo, only : calculate_helmholtz_energy

    ! Precision
    real(dp) :: thr = 1.0e-25_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: ntemp = 8
    real(dp), dimension(8) :: temp = [100.15_dp, 200.0_dp, 300.0_dp, 406.05_dp, &
                                      500.0_dp, 604.0_dp, 730.4_dp, 850.0_dp]
    real(dp), dimension(8) :: lnq = [100.03_dp, 2.45_dp, 3.2_dp, 4.99_dp, &
                                     5.04_dp, 6.0_dp, 7.089_dp, 8.54_dp]
    real(dp), dimension(8) :: a

    ! Expected values
    real(dp) :: expected(8) = [-1.38313459e-19_dp, -6.76517912e-21_dp, -1.32542285e-20_dp, &
                                -2.79745610e-20_dp, -3.47923498e-20_dp, -5.00347125e-20_dp, &
                                -7.14873109e-20_dp, -1.00221296e-19_dp]
    
    ! Call the subroutine
    call calculate_helmholtz_energy(ntemp, temp, lnq, a)

    ! Check the results
    call check(error, a(1), expected(1), thr=thr, rel=.false.)
    call check(error, a(2), expected(2), thr=thr, rel=.false.)
    call check(error, a(3), expected(3), thr=thr, rel=.false.)
    call check(error, a(4), expected(4), thr=thr, rel=.false.)
    call check(error, a(5), expected(5), thr=thr, rel=.false.)
    call check(error, a(6), expected(6), thr=thr, rel=.false.)
    call check(error, a(7), expected(7), thr=thr, rel=.false.)
    call check(error, a(8), expected(8), thr=thr, rel=.false.)
  
  end subroutine test_calc_helmholtz_energy

  ! -----------------------------------------
  ! Unit test for calculate_helmholtz_energy
  ! -----------------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_calc_gibbs_enthalpy(error)
  
    ! Dependencies
    use thermo, only : calculate_gibbs_enthalpy
    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: ntemp = 12
    real(dp), dimension(12) :: temp = [100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp, 423.0_dp, 1100.0_dp, 1200.8_dp]
    real(dp), dimension(12) :: lnq = [1.03_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                     5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, &
                                     9.0_dp, 10.0_dp, 11.0_dp, -12.0_dp]
    real(dp), dimension(12) :: vol = [0.5_dp, 0.6_dp, 4.7_dp, 0.8_dp, &
                                      0.9_dp, 1.0_dp, 1.1_dp, 61.2_dp, &
                                      1.3_dp, 1.4_dp, 1.5_dp, 12.6_dp]
    real(dp) :: press = 1.01325_dp
    real(dp), dimension(12) :: g

    ! Expected values
    real(dp) :: expected(12) = [0.506625_dp, 0.60795_dp, 4.762275_dp, 0.8106_dp, &
                                0.911925_dp, 1.01325_dp, 1.114575_dp, 62.0109_dp, &
                                1.317225_dp, 1.41855_dp, 1.519875_dp, 12.76695_dp]
    
    ! Call the subroutine
    call calculate_gibbs_enthalpy(ntemp, temp, lnq, vol, press, g)

    ! Check the results
    call check(error, g(1), expected(1), thr=thr, rel=.false.)
    call check(error, g(2), expected(2), thr=thr, rel=.false.)
    call check(error, g(3), expected(3), thr=thr, rel=.false.)
    call check(error, g(4), expected(4), thr=thr, rel=.false.)
    call check(error, g(5), expected(5), thr=thr, rel=.false.)
    call check(error, g(6), expected(6), thr=thr, rel=.false.)
    call check(error, g(7), expected(7), thr=thr, rel=.false.)
    call check(error, g(8), expected(8), thr=thr, rel=.false.)
    call check(error, g(9), expected(9), thr=thr, rel=.false.)
    call check(error, g(10), expected(10), thr=thr, rel=.false.)
    call check(error, g(11), expected(11), thr=thr, rel=.false.)
    call check(error, g(12), expected(12), thr=thr, rel=.false.)
  
  end subroutine test_calc_gibbs_enthalpy

  ! -----------------------------------------
  ! Unit test for calculate_internal_energy
  ! -----------------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_calc_internal_energy(error)
  
    ! Dependencies
    use thermo, only : calculate_internal_energy
    ! Precision
    real(dp) :: thr = 1.0e-24_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    integer :: ntemp = 21
    real(dp), dimension(21) :: temp = [100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp, 423.0_dp, 1100.0_dp, 1200.8_dp, &
                                      100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp]
    real(dp), dimension(21) :: dlnq = [1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp, 10.0_dp, 11.0_dp, -12.0_dp, &
                                       1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, 5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp]
    real(dp), dimension(21) :: u


    ! Expected values
    real(dp) :: expected(21) = [1.42633767e-21_dp, -1.84268461e-18_dp, 1.15817216e-17_dp, 4.52997000e-18_dp, &
                                1.95463421e-17_dp, 2.66951869e-16_dp, -5.52415699e-17_dp, 8.51881019e-17_dp, &
                                1.00649298e-16_dp, 2.47038109e-17_dp, 1.83764355e-16_dp, -2.38894320e-16_dp, &
                                1.42633767e-21_dp, -1.84268461e-18_dp, 1.15817216e-17_dp, 4.52997000e-18_dp, &
                                1.95463421e-17_dp, 2.66951869e-16_dp, -5.52415699e-17_dp, 8.51881019e-17_dp, &
                                1.00649298e-16_dp]
    
    ! Call the subroutine
    call calculate_internal_energy(ntemp, temp, dlnq, u)

    ! Check the results
    do i = 1, 21
      call check(error, u(i), expected(i), thr=thr, rel=.false.)
    end do
  
  end subroutine test_calc_internal_energy
  
  end module test_thermo
  
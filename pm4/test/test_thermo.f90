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
      new_unittest("test_calc_internal_energy", test_calc_internal_energy), &
      new_unittest("test_calc_enthalpy", test_calc_enthalpy), &
      new_unittest("test_calc_entropy", test_calc_entropy), &
      new_unittest("test_calc_expansion_coefficient", test_calc_expansion_coefficient), &
      new_unittest("test_calc_cv", test_calc_cv), &
      new_unittest("test_calc_cp", test_calc_cp) &
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

  ! -----------------------------------------
  ! Unit test for calculate_enthalpy
  ! -----------------------------------------
  subroutine test_calc_enthalpy(error)
  
    ! Dependencies
    use thermo, only : calculate_enthalpy
    ! Precision
    real(dp) :: thr = 1.0e-24_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    integer :: ntemp = 18
    real(dp), dimension(18) :: temp = [412.25_dp, 123.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp, 423.0_dp, 1100.0_dp, 1200.8_dp, &
                                      100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp]
    real(dp), dimension(18) :: dlnq = [-2.4_dp, 38.5_dp, 20.43_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp, 10.0_dp, 11.0_dp, -12.0_dp, &
                                       1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, 5.04_dp, 53.0_dp]
    real(dp), dimension(18) :: vol = [0.33_dp, 23.5_dp, 17.32_dp, 0.43_dp, &
                                      8.3_dp, 1.0_dp, 1.1_dp, 61.2_dp, 1.3_dp, 1.4_dp, 1.5_dp, 12.6_dp, &
                                      0.33_dp, 23.5_dp, 17.32_dp, 0.43_dp, 8.3_dp, 1.0_dp]
    real(dp) :: press = 1e-20_dp
    real(dp), dimension(18) :: h


    ! Expected values
    real(dp) :: expected(18) = [-5.62809240e-18_dp, 8.32920629e-18_dp, 7.41152536e-17_dp, 4.53427000e-18_dp, &
                                1.96293421e-17_dp, 2.66961869e-16_dp, -5.52305699e-17_dp, 8.58001019e-17_dp, &
                                1.00662298e-16_dp, 2.47178109e-17_dp, 1.83779355e-16_dp, -2.38768320e-16_dp, &
                                4.72633767e-21_dp, -1.60768461e-18_dp, 1.17549216e-17_dp, 4.53427000e-18_dp, &
                                1.96293421e-17_dp, 2.66961869e-16_dp]
    
    ! Call the subroutine
    call calculate_enthalpy(ntemp, temp, dlnq, vol, press, h)

    ! Check the results
    do i = 1, 18
      call check(error, h(i), expected(i), thr=thr, rel=.false.)
    end do
  
  end subroutine test_calc_enthalpy

  ! -----------------------------------------
  ! Unit test for calculate_entropy
  ! -----------------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_calc_entropy(error)
  
    ! Dependencies
    use thermo, only : calculate_entropy
    ! Precision
    real(dp) :: thr = 1.0e-26_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    integer :: ntemp = 33
    real(dp), dimension(33) :: temp = [412.25_dp, 123.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp, 423.0_dp, 1100.0_dp, 1200.8_dp, &
                                      100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp, 423.0_dp, 1100.0_dp, 1200.8_dp, &
                                      100.15_dp, 233.4_dp, 512.0_dp, 406.05_dp, &
                                      530.0_dp, 604.0_dp, 730.4_dp, 850.0_dp, &
                                      900.0_dp]
    real(dp), dimension(33) :: lnq = [-2.4_dp, 38.5_dp, 20.43_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp,  &
                                       9.0_dp, 10.0_dp, 11.0_dp, -12.0_dp, &
                                       1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp, &
                                       10.0_dp, 11.0_dp, -12.0_dp, 1.03e-2_dp, &
                                       -2.45_dp, 3.2_dp, 1.99_dp, 5.04_dp, 53.0_dp, &
                                       -7.5_dp, 8.54_dp, 9.0_dp]
    real(dp), dimension(33) :: dlnq = [1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, &
                                       9.0_dp, 10.0_dp, 11.0_dp, -12.0_dp, &
                                       1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp, &
                                       10.0_dp, 11.0_dp, -12.0_dp, &
                                       1.03e-2_dp, -2.45_dp, 3.2_dp, 1.99_dp, &
                                       5.04_dp, 53.0_dp, -7.5_dp, 8.54_dp, 9.0_dp]
    real(dp), dimension(33) :: s


    ! Expected values
    real(dp) :: expected(33) = [2.54891930e-23_dp, -3.64256573e-21_dp,  2.29026165e-20_dp,  1.11836626e-20_dp, &
                                3.69494754e-20_dp,  4.42705038e-19_dp, -7.57354899e-20_dp,  1.00339204e-19_dp, &
                                1.11956811e-19_dp,  5.85395091e-20_dp,  1.67210376e-19_dp, -1.99111647e-19_dp, &
                                1.43842205e-23_dp, -7.92878993e-21_dp,  2.26647307e-20_dp,  1.11836626e-20_dp, &
                                3.69494754e-20_dp,  4.42705038e-19_dp, -7.57354899e-20_dp,  1.00339204e-19_dp, &
                                1.11956811e-19_dp,  5.85395091e-20_dp,  1.67210376e-19_dp, -1.99111647e-19_dp, &
                                1.43842205e-23_dp, -7.92878993e-21_dp,  2.26647307e-20_dp,  1.11836626e-20_dp, &
                                3.69494754e-20_dp,  4.42705038e-19_dp, -7.57354899e-20_dp,  1.00339204e-19_dp, &
                                1.11956811e-19_dp]
    
    ! Call the subroutine
    call calculate_entropy(ntemp, temp, lnq, dlnq, s)

    ! Check the results
    do i = 1, 33
      call check(error, s(i), expected(i), thr=thr, rel=.false.)
    end do
  
  end subroutine test_calc_entropy

  ! -----------------------------------------
  ! Unit test for  calculate_expansion_coefficient
  ! -----------------------------------------
  ! Adds the part of the system partition function arising from particle
  ! indistinguishability. 
  subroutine test_calc_expansion_coefficient(error)
  
    ! Dependencies
    use thermo, only :  calculate_expansion_coefficient
    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    integer :: ntemp = 10
    real(dp), dimension(10) :: vol = [12.42_dp, 1.93_dp, 5.44_dp, 3.23_dp, 4.23_dp, &
                                       5.23_dp, 6.23_dp, 7.23_dp, 8.23_dp, 9.23_dp]
    real(dp), dimension(10) :: dvol = [-1.3_dp, 300.4_dp, 23.3_dp, -12.44_dp, 25.32_dp, &
                                      -212.42_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
    
    real(dp), dimension(10) :: alpha


    ! Expected values
    real(dp) :: expected(10) = [-1.04669887e-01_dp,  1.55647668e+02_dp,  4.28308824e+00_dp, -3.85139319e+00_dp, &
                                 5.98581560e+00_dp, -4.06156788e+01_dp,  1.60513644e-01_dp,  2.76625173e-01_dp, &
                                 3.64520049e-01_dp,  4.33369447e-01_dp]
    
    ! Call the subroutine
    call  calculate_expansion_coefficient(ntemp, vol, dvol, alpha)

    ! Check the results
    do i = 1, 10
      call check(error, alpha(i), expected(i), thr=thr, rel=.false.)
    end do
  
  end subroutine test_calc_expansion_coefficient

  ! -----------------------------------------
  ! Unit test for  calculate_cv
  ! -----------------------------------------
  ! Calculates the heat capacity at constant volume.
  ! BUT : The volume changes at every temperature.
  !       It is not constant and thus it is not possible to calculate the heat capacity at constant volume.
  subroutine test_calc_cv(error)
    
      ! Dependencies
      use thermo, only : calculate_cv
      ! Precision
      real(dp) :: thr = 1.0e-25_dp
      
      ! Arguments
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      integer :: ntemp = 10
      real(dp), dimension(10) :: temp = [12.42_dp, 1.93_dp, 5.44_dp, 3.23_dp, 4.23_dp, &
                                         5.23_dp, 6.23_dp, 7.23_dp, 8.23_dp, 9.23_dp]
      real(dp), dimension(10) :: dlnq = [-1.3_dp, 300.4_dp, 23.3_dp, -12.44_dp, 25.32_dp, &
                                        -212.42_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      real(dp), dimension(10) :: ddlnq = [38.2_dp, -0.2_dp, 0.0_dp, 100.3_dp, -4.2_dp, &
                                         3.2_dp, 2.3_dp, 1.3_dp, 0.3_dp, 0.1_dp]
      real(dp), dimension(10) :: cv

      ! Expected values
      real(dp) :: expected(10) = [8.09101959e-20_dp, 1.59989448e-20_dp, 3.49999993e-21_dp, 1.33378609e-20_dp, &
                                   1.91988908e-21_dp, -2.94683460e-20_dp, 1.40452947e-21_dp, 1.33750035e-21_dp, &
                                   9.62310419e-22_dp, 1.13709255e-21_dp]

      ! Call the subroutine
      call calculate_cv(ntemp, temp, dlnq, ddlnq, cv)

      ! Check the results
      do i = 1, 10
        call check(error, cv(i), expected(i), thr=thr, rel=.false.)
      end do

  end subroutine test_calc_cv

  !-----------------------------------------
  ! Unit test for calculate_cp
  !-----------------------------------------
  ! Calculates the heat capacity at constant pressure.
  ! BUT : Is calculated using the heat capacity at constant volume, which is not possible to calculate.

  subroutine test_calc_cp(error)
    
      ! Dependencies
      use thermo, only : calculate_cp
      ! Precision
      real(dp) :: thr = 1.0e-13_dp
      
      ! Arguments
      type(error_type), allocatable, intent(out) :: error
      integer :: i
      integer :: ntemp = 10
      real(dp), dimension(10) :: temp = [-1.3_dp, 300.4_dp, 23.3_dp, -12.44_dp, 25.32_dp, &
                                         -212.42_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      real(dp), dimension(10) :: dlnq = [5.2_dp, 3.8_dp, 2.3_dp, 1.3_dp, 4.2_dp, &
                                        3.3_dp, 2.8_dp, 1.3_dp, 3.2_dp, 3.3_dp]
      real(dp), dimension(10) :: ddlnq = [1.3_dp, 4.2_dp, 28.3_dp, 2.8_dp, -1.3_dp, &
                                         3.2_dp, 33.3_dp, 1.3_dp, 8.3_dp, 1.3_dp]
      real(dp), dimension(10) :: dvol = [5.3e-8_dp, 2.1e-7_dp, -2.6e-9_dp, 1.0e-8_dp, 1.3e-6_dp, &
                                        2.3e-8_dp, 3.3e-8_dp, 4.3e-6_dp, 5.3e-9_dp, 6.3e-8_dp]
      real(dp) :: press = 1.01325_dp
      real(dp), dimension(10) :: cp

      ! Expected values
      real(dp) :: expected(10) = [5.370225e-08_dp, 2.127825e-07_dp, -2.634450e-09_dp, 1.013250e-08_dp, 1.317225e-06_dp, &
                                   2.330475e-08_dp, 3.343725e-08_dp, 4.356975e-06_dp, 5.370225e-09_dp, 6.383475e-08_dp]

      ! Call the subroutine
      call calculate_cp(ntemp, temp, dlnq, ddlnq, dvol, press, cp)

      ! Check the results
      do i = 1, 10
        call check(error, cp(i), expected(i), thr=thr, rel=.false.)
      end do

    end subroutine test_calc_cp
 
  
  end module test_thermo
  
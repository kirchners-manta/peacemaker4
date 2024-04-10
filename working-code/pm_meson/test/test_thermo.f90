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
      new_unittest("test_calc_lnq_sys", test_calc_lnq_sys) &
      ]
  
  end subroutine collect_thermo
  
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

  
  end module test_thermo
  
! This module tests the subroutine solve_polynomial3 which is used to solve the volume polynomial.

module test_vol_polynomial
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds, only : dp
    implicit none
    private
  
    public :: collect_vol_polynomial
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_vol_polynomial(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_solve_polynomial3", test_solve_polynomial3), &
      new_unittest("test_solve_polynomial3_2", test_solve_polynomial3_2) &
      ]
  
  end subroutine collect_vol_polynomial
  
  !------------------------------------------
  ! Unit tests for solve_polynomial3
  !------------------------------------------
  subroutine test_solve_polynomial3(error)
    ! 2x³ - 3x² - 3x + 2 = 0
    ! x = -1, 2, 1/2

    ! Dependencies
    use polynomial, only : solve_polynomial3

    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(4) :: coeffs = [2.0_dp, -3.0_dp, -3.0_dp, 2.0_dp]
    complex(dp), dimension(3) :: roots

    ! Expected values
    real(dp), dimension(3) :: expected_real = [2.0_dp, -1.0_dp, 0.5_dp]
    real(dp), dimension(3) :: expected_imag = [0.0_dp, 0.0_dp, 0.0_dp]
    
    ! Call the subroutine
    call solve_polynomial3(coeffs, roots)

    ! Check the results
    call check(error, realpart(roots(1)), expected_real(1), thr=thr, rel=.false.)
    call check(error, realpart(roots(2)), expected_real(2), thr=thr, rel=.false.)
    call check(error, realpart(roots(3)), expected_real(3), thr=thr, rel=.false.)
    call check(error, imagpart(roots(1)), expected_imag(1), thr=thr, rel=.false.)
    call check(error, imagpart(roots(2)), expected_imag(2), thr=thr, rel=.false.)
    call check(error, imagpart(roots(3)), expected_imag(3), thr=thr, rel=.false.)

  end subroutine test_solve_polynomial3

  subroutine test_solve_polynomial3_2(error)
    ! 5x³ + x² - 2x + 4 = 0
    ! x = -0.84566, 0.32283 + 0.91749 i, 0.32283 - 0.91749 i

    ! Dependencies
    use polynomial, only : solve_polynomial3

    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(4) :: coeffs = [4.0_dp, 2.0_dp, 1.0_dp, 5.0_dp]
    complex(dp), dimension(3) :: roots

    ! Expected values
    real(dp), dimension(3) :: expected_real = [-0.84566_dp, 0.32283_dp, 0.32283_dp]
    real(dp), dimension(3) :: expected_imag = [0.0_dp, 0.91749_dp, -0.91749_dp]
    
    ! Call the subroutine
    call solve_polynomial3(coeffs, roots)

    ! Check the results
    call check(error, realpart(roots(1)), expected_real(1), thr=thr, rel=.false.)
    call check(error, realpart(roots(2)), expected_real(2), thr=thr, rel=.false.)
    call check(error, realpart(roots(3)), expected_real(3), thr=thr, rel=.false.)
    call check(error, imagpart(roots(1)), expected_imag(1), thr=thr, rel=.false.)
    call check(error, imagpart(roots(2)), expected_imag(2), thr=thr, rel=.false.)
    call check(error, imagpart(roots(3)), expected_imag(3), thr=thr, rel=.false.)

  end subroutine test_solve_polynomial3_2
  
  
  end module test_vol_polynomial
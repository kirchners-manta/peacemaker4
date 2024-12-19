! This module tests the subroutine solve_polynomial3 which is used to solve the volume polynomial.
! The volume polynomial is a cubic polynomial of the form a_3 x³ + a_2 x² + a_1 x + a_0 = 0.
! The soubroutine fails if the coefficient of the cubic term is zero since we devide by it.
! The subroutine fails if a_1 and a_2 are both zero.

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
      new_unittest("test_solve_polynomial3_2", test_solve_polynomial3_2), &
      new_unittest("test_solve_polynomial3_3", test_solve_polynomial3_3) &
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

  subroutine test_solve_polynomial3_3(error)
    ! 3.4x³ + 2.8x² - 5.1x + 2.1 = 0

    ! Dependencies
    use polynomial, only : solve_polynomial3

    ! Precision
    real(dp) :: thr = 1.0e-5_dp
    
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(4) :: coeffs = [2.1_dp, -5.1_dp, 2.8_dp, 3.4_dp]
    complex(dp), dimension(3) :: roots

    ! Expected values
    real(dp), dimension(3) :: expected_real = [0.502518460445945_dp, -1.82856633265660_dp, 0.502518460445945_dp]
    real(dp), dimension(3) :: expected_imag = [-0.291979234759594_dp, 0.0_dp, 0.291979234759594_dp]
    
    ! Call the subroutine
    call solve_polynomial3(coeffs, roots)

    ! Check the results
    call check(error, realpart(roots(1)), expected_real(1), thr=thr, rel=.false.)
    call check(error, realpart(roots(2)), expected_real(2), thr=thr, rel=.false.)
    call check(error, realpart(roots(3)), expected_real(3), thr=thr, rel=.false.)
    call check(error, imagpart(roots(1)), expected_imag(1), thr=thr, rel=.false.)
    call check(error, imagpart(roots(2)), expected_imag(2), thr=thr, rel=.false.)
    call check(error, imagpart(roots(3)), expected_imag(3), thr=thr, rel=.false.)

  end subroutine test_solve_polynomial3_3
    
  end module test_vol_polynomial
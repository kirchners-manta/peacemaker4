module test_auxiliary
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds, only : dp
    !> imports go here
    implicit none
    private
  
    public :: collect_auxiliary
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_auxiliary(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_ln_factorial", test_ln_factorial), &
      new_unittest("test_ln_factorial_2", test_ln_factorial_2, should_fail=.true.), &
      new_unittest("test_ln_factorial_3", test_ln_factorial_3), &
      new_unittest("test_ln_factorial_4", test_ln_factorial_4) &
      ]
  
  end subroutine collect_auxiliary
  

!--------------------------------------
! Unit tests for calculation of ln(n!)
!--------------------------------------
  ! Gives infinity for negative input.
  ! This is no problem since it is only called with non-negative integers. (Populations)
  subroutine test_ln_factorial(error)
    use auxiliary, only : ln_factorial
    !> Precision of the test
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x = 0.0_dp

    ! Expected result
    real(dp) :: expected = 0.0_dp

    ! Computed result
    real(dp) :: computed

    ! Call the function
    computed = ln_factorial(x)

    call check(error, computed, expected, thr=thr, rel=.false.)
  end subroutine test_ln_factorial

  subroutine test_ln_factorial_2(error)
    use auxiliary, only : ln_factorial
    !> Precision of the test
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x = -1.0_dp

    ! Expected result
    real(dp) :: expected = 0.0_dp

    ! Computed result
    real(dp) :: computed

    ! Call the function
    computed = ln_factorial(x)

    call check(error, computed, expected, thr=thr, rel=.false.)
  end subroutine test_ln_factorial_2
  
  subroutine test_ln_factorial_3(error)
    use auxiliary, only : ln_factorial
    !> Precision of the test
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x = 300.0_dp

    ! Expected result
    real(dp) :: expected = 1414.9058499450679885_dp

    ! Computed result
    real(dp) :: computed

    ! Call the function
    computed = ln_factorial(x)

    call check(error, computed, expected, thr=thr, rel=.false.)
  end subroutine test_ln_factorial_3

  subroutine test_ln_factorial_4(error)
    use auxiliary, only : ln_factorial
    !> Precision of the test
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: x = 0.021_dp

    ! Expected result
    real(dp) :: expected = -0.0117625_dp

    ! Computed result
    real(dp) :: computed

    ! Call the function
    computed = ln_factorial(x)

    call check(error, computed, expected, thr=thr, rel=.false.)
  end subroutine test_ln_factorial_4
  
  end module test_auxiliary
  
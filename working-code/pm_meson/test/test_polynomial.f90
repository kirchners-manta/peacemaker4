module test_polynomial
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use kinds, only : dp
  implicit none
  private

  public :: collect_polynomial

contains

!> Collect all exported unit tests
subroutine collect_polynomial(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test_newton", test_newton), &
    new_unittest("test1.2", test_invalid, should_fail=.true.) &
    ]

end subroutine collect_polynomial

!----------------------------------------
! Unit test for Newton-Raphson algorithm
!----------------------------------------
subroutine test_newton(error)
  ! Example system:
  ! f_1(x,y) = x^2 + y^2 - 1
  ! f_2(x,y) = 2x - y
  ! Solution: (1/sqrt(5), 2/sqrt(5))

  ! f(x) = 1 - x 

  ! Dependency
  use polynomial, only : newton
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 1                               ! Number of the polynomials
  integer, dimension(n) :: d = [1]                          ! Degree of the polynomials
  integer :: l = 1                                          ! Length of the coefficient array
  real(dp), dimension(0:1) :: coeffs = [1.0, -1.0]          ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [1.0]                      ! Initial guess
  integer :: iter = 500                                     ! Maximum number of iterations
  logical :: success                                        ! Success flag

  ! Expected result
  !real(dp), dimension(n) :: expected = [1.0/sqrt(5.0), 2.0/sqrt(5.0)]
  real(dp) :: expected = 1.0

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success)
  write(*,*) "Success: ", success
  write(*,*) "x0: ", x0
  write(*,*) "expected: ", expected

  ! Check the result
  call check(error, x0(1), expected, thr=0.01_dp, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  
  call check(error, 1, 0)
  if (allocated(error)) return
end subroutine test_invalid

end module test_polynomial

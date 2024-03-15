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
    new_unittest("test_newton_2", test_newton_2) &
    !new_unittest("test_newton_3", test_newton_3) &
    ]

end subroutine collect_polynomial

!----------------------------------------
! Unit test for Newton-Raphson algorithm
!----------------------------------------
subroutine test_newton(error)
  ! Example polynomial:
  ! f(x) = 2 - 9x 
  ! Solution: 2/9

  ! Dependency
  use polynomial, only : newton
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 1                               ! Number of the polynomials
  integer, dimension(n) :: d = [1]                          ! Degree of the polynomials
  integer :: l = 2                                          ! Length of the coefficient array
  real(dp), dimension(2) :: coeffs = [2.0, -9.0]            ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [0.5]                      ! Initial guess
  integer :: iter = 500                                     ! Maximum number of iterations
  logical :: success                                        ! Success flag

  ! Expected result
  real(dp) :: expected = 2.0/9.0

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success)
  write(*,*) "Success: ", success

  ! Check the result
  call check(error, x0(1), expected, thr=0.001_dp, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton

subroutine test_newton_2(error)
  ! Example system:
  ! f_1(x,y) = -1 + x + y
  ! f_2(x,y) = -3 +2x +4y
  ! Solution: (0.5, 0.5)

  ! Dependency
  use polynomial, only : newton
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                                                       ! Number of the polynomials
  integer, dimension(n) :: d = [2, 2]                                                ! Degree of each dimension
  integer :: l = 8                                                                  ! Length of the coefficient array
  real(dp), dimension(8) :: coeffs = [-1.0, 1.0, 1.0, 0.0, -3.0, 2.0, 4.0, 0.0]           ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [0.54, 0.55]                                   ! Initial guess
  integer :: iter = 500                                                             ! Maximum number of iterations
  logical :: success                                                                ! Success flag

  ! Expected result
  real(dp), dimension(n) :: expected = [0.5, 0.5]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success)
  write(*,*) "Success: ", success
  write(*,*) "x0: ", x0

  ! Check the result
  call check(error, x0(1), expected(1), thr=0.01_dp, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_2

subroutine test_newton_3(error)
  ! Example system:
  ! f_1(x,y) = -1 + 0x + x^2 + 0xy + 0y + y^2
  ! f_2(x,y) = 0 + 2x - y + 0xy 
  ! Solution: (1/sqrt(5), 2/sqrt(5))

  ! Dependency
  use polynomial, only : newton
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                               ! Number of the polynomials
  integer, dimension(n) :: d = [2,2]                        ! Degree of the polynomials
  integer :: l = 18                                         ! Length of the coefficient array
  real(dp), dimension(18) :: coeffs = [-1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0]         ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [1.0/sqrt(5.0), 2.0/sqrt(5.0)]         ! Initial guess
  integer :: iter = 500                                     ! Maximum number of iterations
  logical :: success                                        ! Success flag

  ! Expected result
  real(dp), dimension(n) :: expected = [1.0/sqrt(5.0), 2.0/sqrt(5.0)]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success)
  write(*,*) "Success: ", success
  write(*,*) "x0: ", x0

  ! Check the result
  call check(error, x0(1), expected(1), thr=0.01_dp, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_3

end module test_polynomial

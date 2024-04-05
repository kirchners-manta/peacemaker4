! In this module the Newton and Horner algortihms for solving the population polynomial 
! for obtaining the monomer populations are tested.
! Additionally the subroutine calculating remaining populations is tested.

module test_pop_polynomial
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use kinds, only : dp
  implicit none
  private

  public :: collect_pop_polynomial

contains

!> Collect all exported unit tests
subroutine collect_pop_polynomial(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test_horner", test_horner), &
    new_unittest("test_newton", test_newton), &
    new_unittest("test_newton_2", test_newton_2), &
    new_unittest("test_newton_3", test_newton_3), &
    new_unittest("test_newton_4", test_newton_4), &
    new_unittest("test_newton_5", test_newton_5), &
    new_unittest("test_newton_6", test_newton_6), &
    new_unittest("test_newton_7", test_newton_7), &
    new_unittest("test_calc_remaining_populations", test_calc_remaining_populations), &
    new_unittest("test_calc_remaining_populations_2", test_calc_remaining_populations_2) &
    ]

end subroutine collect_pop_polynomial

!----------------------------------------
! Unit test for Horner's algorithm
!----------------------------------------
subroutine test_horner(error)
  ! Example polynomial:
  ! f(x,y) = -2 + x + x² + 7x³ + 2y + xy + 8x²y - 2x³y - 5y² +xy² - x²y² + 3x³y² - y³ + xy³  + x²y² + 12x³y³

  ! Dependency
  use polynomial, only : horner
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                                                 ! Number of the polynomials
  integer, dimension(n) :: d = [3, 3]                                         ! Degree of the polynomials
  integer :: l = 16                                                           ! Length of the coefficient array
  real(dp), dimension(16) :: c = [-2.0, 1.0, 1.0, 7.0, 2.0, 1.0, 8.0, -2.0, &
                                 -5.0, 1.0, -1.0, 3.0, -1.0, 1.0, 1.0, 12.0]  ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [0.2, 0.3]                                   ! Initial guess
  real(dp)  :: p                                                              ! value of polynomial at x
  real(dp), dimension(n) :: pdiff                                             ! derivative of polynomial at x

  ! Expected result
  ! Value of the polynomial at x0
  real(dp) :: p_expected = -1.404168_dp
  ! Value of the x- and y-derivative of the polynomial at x0 respectively
  real(dp), dimension(n) :: pdiff_expected = [3.59108_dp, -0.5648799999999994_dp]

  ! Call Horner's algorithm
  call horner(n, d, l, c, x0, p, pdiff)

  ! Check the results
  call check(error, p, p_expected, thr=thr, rel=.false.)
  call check(error, pdiff(1), pdiff_expected(1), thr=thr, rel=.false.)
  call check(error, pdiff(2), pdiff_expected(2), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_horner


!----------------------------------------
! Unit test for Newton algorithm
!----------------------------------------
! 1D polynomial, degree 1
subroutine test_newton(error)
  ! Example polynomial:
  ! f(x) = 2 - 9x 
  ! Solution: 2/9

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 1                               ! Number of the polynomials
  integer, dimension(n) :: d = [1]                          ! Degree of the polynomials
  integer :: l = 2                                          ! Length of the coefficient array
  real(dp), dimension(2) :: coeffs = [2.0, -9.0]            ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [0.5_dp]                   ! Initial guess
  integer :: iter = 500                                     ! Maximum number of iterations
  logical :: success                                        ! Success flag
  integer, dimension(1) :: monomer  = [1]                   ! Monomers

  ! Expected result
  real(dp) :: expected = 0.222222222222_dp

  ! Call the Newton algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected, thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton

! 1D polynomial, degree 4
subroutine test_newton_2(error)
  ! Example system:
  ! f_1(x) = 2x⁴ + 4x³ + x²  + 3x - 2.375

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 1                                     ! Number of the polynomials
  integer, dimension(n) :: d = [4]                                ! Degree of the polynomials
  integer :: l = 5                                                ! Length of the coefficient array
  real(dp), dimension(5) :: coeffs = [-2.375, 3.0, 1.0, 4.0, 2.0] ! Coefficients of the polynomials   
  real(dp), dimension(n) :: x0 = [0.4_dp]                         ! Initial guess
  integer :: iter = 5                                             ! Maximum number of iterations
  logical :: success                                              ! Success flag
  integer, dimension(1) :: monomer  = [1]                         ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.5_dp]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_2

! 2D polynomial, degree 1, 1
subroutine test_newton_3(error)
  ! Example system:
  ! f_1(x,y) = -1 + x + y
  ! f_2(x,y) = -3 +2x +4y
  ! Solution: (0.5, 0.5)

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                                                    ! Number of the polynomials
  integer, dimension(n) :: d = [1, 1]                                            ! Degree of each dimension
  integer :: l = 8                                                               ! Length of the coefficient array
  real(dp), dimension(8) :: coeffs = [-1.0, 1.0, 1.0, 0.0, -3.0, 2.0, 4.0, 0.0]  ! Coefficients of the polynomials
  real(dp), dimension(n) :: x0 = [0.3_dp, 0.7_dp]                                ! Initial guess
  integer :: iter = 1000                                                         ! Maximum number of iterations
  logical :: success                                                             ! Success flag
  integer, dimension(2) :: monomer  = [1, 1]                                     ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.5_dp, 0.5_dp]

  ! Call the Newton algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)
  call check(error, x0(2), expected(2), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_3

! 2D polynomial, degree 2, 2
subroutine test_newton_4(error)
  ! Example system:
  ! f_1(x,y) = -1 + 0x + x² + 0xy + 0y + y²
  ! f_2(x,y) = 0 + 2x - y + 0xy 
  ! Solution: (1/sqrt(5), 2/sqrt(5))

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                                   ! Number of the polynomials
  integer, dimension(n) :: d = [2, 2]                           ! Degree of the polynomials
  integer :: l = 18                                             ! Length of the coefficient array
  real(dp), dimension(18) :: coeffs = [-1.0, 0.0, 1.0, 0.0, &   ! Coefficients of the polynomials
                                        0.0, 0.0, 1.0, 0.0, &
                                        0.0, 0.0, 2.0, 0.0, &
                                       -1.0, 0.0, 0.0, 0.0, &
                                        0.0, 0.0]         
  real(dp), dimension(n) :: x0 = [0.3_dp, 0.7_dp]               ! Initial guess
  integer :: iter = 500                                         ! Maximum number of iterations
  logical :: success                                            ! Success flag
  integer, dimension(2) :: monomer  = [1, 1]                    ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.4472135955_dp, 0.894427191_dp]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)
  call check(error, x0(2), expected(2), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_4

! 2D polynomial, degree 3, 3
subroutine test_newton_5(error)
  ! Example system:
  ! f_1(x,y) = 2 + 3x - x^2 + 5x^3 - 15y + 2y^2 + 8y^3 + 20xy + x^2y + 12xy^2 + 100x^2y^2 + 1000x^3y^2 + 598x^2y^3 - 2105x^3y^3
  ! f_2(x,y) = x + x^2 + 2x^3 + y - 5y^2 -27.7xy^2 - 3x^2y^2
  ! Solution: (0.1, 0.2)

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-2_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 2                                   ! Number of the polynomials
  integer, dimension(n) :: d = [3,3]                            ! Degree of the polynomials
  integer :: l = 32                                             ! Length of the coefficient array
  real(dp), dimension(32) :: coeffs = [2.0_dp, 3.0_dp, -1.0_dp, 5.0_dp, &   ! Coefficients of the polynomials
                                      -15.0_dp, 20.0_dp, 1.0_dp, 0.0_dp, &
                                      2.0_dp, 12.0_dp, 100.0_dp, 1000.0_dp, &
                                      8.0_dp, 0.0_dp, 598.0_dp, 2105.0_dp, &
                                      0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, & ! second polynomial
                                      0.0_dp, 0.0_dp, 0.0_dp, -5.0_dp, -27.7_dp, &
                                      -3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]                         
  real(dp), dimension(n) :: x0 = [0.5_dp, 0.5_dp] ! Initial guess
  integer :: iter = 1000                                         ! Maximum number of iterations
  logical :: success                                            ! Success flag
  integer, dimension(2) :: monomer  = [1, 1]                   ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.100000000_dp, 0.200000000_dp]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)
  call check(error, x0(2), expected(2), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_5

! 3D polynomial, degree 3, 3, 3
subroutine test_newton_6(error)
  ! Example system:
  ! f_1(x,y,z) = -1.2568 + x + x² + 5x³ - y + 3y² + 4z + xy + xyz²
  ! f_2(x,y,z) = -0.3129 + x³ + 8y² - z³ + yz² + 3x³z
  ! f_3(x,y,z) = 0.983 + 5x³ - y + 8y² - 4z + z² + x²y

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 3                                   ! Number of the polynomials
  integer, dimension(n) :: d = [3,2,3]                            ! Degree of the polynomials
  integer :: l = 144                                             ! Length of the coefficient array
  real(dp), dimension(144) :: coeffs = [-1.2568_dp, 1.0_dp, 1.0_dp, 5.0_dp, -1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, & ! Firsr polynomial
                                         3.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & 
                                         -0.3129_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! Second poynomial
                                         8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 3.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.983_dp, 0.0_dp, 0.0_dp, 5.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, & ! Third poynomial
                                         8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                                         0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]       
  
  

  real(dp), dimension(n) :: x0 = [0.5_dp, 0.5_dp, 0.5_dp] ! Initial guess
  integer :: iter = 1000                                        ! Maximum number of iterations
  logical :: success                                            ! Success flag
  integer, dimension(3) :: monomer  = [1, 1, 1]                 ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.100000000_dp, 0.200000000_dp, 0.300000000_dp]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)
  call check(error, x0(2), expected(2), thr=thr, rel=.false.)
  call check(error, x0(3), expected(3), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_6

! 4D polynomial, degree 4, 1, 2, 2
subroutine test_newton_7(error)
  ! Example system:
  ! f_1(x,y,z) = -1.2568 + x + x² + 5x³ - y + 3y² + 4z + xy + xyz²
  ! f_2(x,y,z) = -0.3129 + x³ + 8y² - z³ + yz² + 3x³z
  ! f_3(x,y,z) = 0.983 + 5x³ - y + 8y² - 4z + z² + x²y

  ! Dependency
  use polynomial, only : newton
    !> Precision for the tests
  real(dp) :: thr = 1.0e-5_dp
  ! Arguments
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: n = 4                                    ! Number of the polynomials
  integer, dimension(n) :: d = [4, 1, 2, 2]                      ! Degree of the polynomials
  integer :: l = 360                                             ! Length of the coefficient array
  real(dp), dimension(360) :: coeffs = [0.6178516892_dp, 1.0_dp, 2.0_dp, 3.0_dp, -5.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, & ! 1  First
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 2
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -5.0_dp, 0.0_dp, & ! 3
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 4
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 5
                                        0.0_dp, -2.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 6
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 7
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 8
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 9
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.3_dp, & ! 10
                                        0.3488084515_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, -4.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 1  Second
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 2
                                        0.0_dp, 0.0_dp, 0.0_dp, 6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 3
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, & ! 4
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 5
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 6
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 7
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 8
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 9
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 10
                                        -0.17944167_dp, 0.0_dp, 3.0_dp, -1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, & ! 1  Third
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 2
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 3
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 4
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 5
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 6
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 7
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 8
                                        0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, & ! 9
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 10
                                        -0.1328316336_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, & ! 1  Fourth
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 2
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 3
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 4
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 5
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, -7.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 6
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 7
                                        6.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 8
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, & ! 9
                                        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 8.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]   ! 10
  
  

  real(dp), dimension(n) :: x0 = [0.5_dp, 0.5_dp, 0.5_dp, 0.5_dp] ! Initial guess
  integer :: iter = 500                                           ! Maximum number of iterations
  logical :: success                                              ! Success flag
  integer, dimension(4) :: monomer  = [1, 1, 1, 1]                ! Monomers

  ! Expected result
  real(dp), dimension(n) :: expected = [0.2300000000_dp, 0.5300000000_dp, 0.1900000000_dp, 0.6400000000_dp]

  ! Call the Newton-Raphson algorithm
  call newton(n, d, l, coeffs, x0, iter, success, monomer)

  ! Check the result
  call check(error, x0(1), expected(1), thr=thr, rel=.false.)
  call check(error, x0(2), expected(2), thr=thr, rel=.false.)
  call check(error, x0(3), expected(3), thr=thr, rel=.false.)
  call check(error, x0(4), expected(4), thr=thr, rel=.false.)

  if (allocated(error)) return
end subroutine test_newton_7

!-----------------------------------------------
! Unit test for calculate_remaining_populations
!-----------------------------------------------
subroutine test_calc_remaining_populations(error)
  use qce, only : calculate_remaining_populations
  use cluster, only : cluster_t
  use partition_functions, only : pf_t
  !> Precision for the tests
  real(dp) :: thr = 1.0e-10_dp

  ! Arguments 
  type(error_type), allocatable, intent(out) :: error
  type(pf_t), dimension(:), allocatable :: lnq 
  real(dp), dimension(:), allocatable :: monomer_populations
  real(dp), dimension(:), allocatable :: populations
  type(cluster_t), dimension(:), allocatable :: cluster_set
  integer, dimension(:), allocatable :: mono                

  ! Expected result
  real(dp), dimension(6) :: expected = [1.000000000000e-02_dp, 2.195246544376e-04_dp, 2.000000000000e-01_dp, &
                                        5.443077197365e-06_dp, 8.118571642079e-13_dp, 3.000000000000e-02_dp]

  ! Allocate and initialize 
  allocate(cluster_set(6))
  cluster_set(1)%composition = [1, 0, 0]
  cluster_set(2)%composition = [1, 2, 0]
  cluster_set(3)%composition = [0, 1, 0]
  cluster_set(4)%composition = [1, 1, 1]
  cluster_set(5)%composition = [2, 1, 3]
  cluster_set(6)%composition = [0, 0, 1]

  allocate(mono(3))
  mono = [1, 3, 6] ! Positions of the monomers in the cluster set

  allocate(monomer_populations(size(mono)))
  monomer_populations = [0.01_dp, 0.2_dp, 0.03_dp]

  allocate(lnq(size(cluster_set)))
  lnq(1)%qtot = 1.2_dp
  lnq(2)%qtot = 3.2_dp
  lnq(3)%qtot = 1.3_dp
  lnq(4)%qtot = 1.2_dp
  lnq(5)%qtot = 0.5_dp
  lnq(6)%qtot = 1.1_dp

  allocate(populations(size(cluster_set)))

  ! Call the function
  call calculate_remaining_populations(populations, monomer_populations, lnq, cluster_set, mono)

  ! Check the result
  call check(error, populations(1), expected(1), thr=thr, rel=.false.)
  call check(error, populations(2), expected(2), thr=thr, rel=.false.)
  call check(error, populations(3), expected(3), thr=thr, rel=.false.)
  call check(error, populations(4), expected(4), thr=thr, rel=.false.)
  call check(error, populations(5), expected(5), thr=thr, rel=.false.)
  call check(error, populations(6), expected(6), thr=thr, rel=.false.)

end subroutine test_calc_remaining_populations

subroutine test_calc_remaining_populations_2(error)
  use qce, only : calculate_remaining_populations
  use cluster, only : cluster_t
  use partition_functions, only : pf_t
  !> Precision for the tests
  real(dp) :: thr = 10e-10_dp

  ! Arguments 
  type(error_type), allocatable, intent(out) :: error
  type(pf_t), dimension(:), allocatable :: lnq 
  real(dp), dimension(:), allocatable :: monomer_populations
  real(dp), dimension(:), allocatable :: populations
  type(cluster_t), dimension(:), allocatable :: cluster_set
  integer, dimension(:), allocatable :: mono                

  ! Expected result
  real(dp), dimension(9) :: expected = [2.000000000000e+023_dp, 1.088297856368e+171_dp, 2.351085132193e+161_dp, &
                                        1.128067738370e+082_dp, 1.200000000000e+023_dp, 2.523225797410e+148_dp, &
                                        3.100000000000e+022_dp, 1.142209387045e+164_dp, 5.400000000000e+020_dp]

  ! Allocate and initialize 
  allocate(cluster_set(9))
  cluster_set(1)%composition = [1, 0, 0, 0]
  cluster_set(2)%composition = [1, 2, 0, 5]
  cluster_set(3)%composition = [0, 1, 40, 2]
  cluster_set(4)%composition = [1, 1, 0, 1]
  cluster_set(5)%composition = [0, 1, 0, 0]
  cluster_set(6)%composition = [2, 0, 3, 1]
  cluster_set(7)%composition = [0, 0, 1, 0]
  cluster_set(8)%composition = [2, 3, 1, 1]
  cluster_set(9)%composition = [0, 0, 0, 1]

  allocate(mono(4))
  mono = [1, 5, 7, 9] ! Positions of the monomers in the cluster set

  allocate(monomer_populations(size(mono)))
  monomer_populations = [2.0e23_dp, 1.2e23_dp, 3.1e22_dp, 5.4e20_dp]

  allocate(lnq(size(cluster_set)))
  lnq(1)%qtot = -3.2_dp
  lnq(2)%qtot = 5.0_dp
  lnq(3)%qtot = 25.0_dp
  lnq(4)%qtot = 36.2_dp
  lnq(5)%qtot = 4.0_dp
  lnq(6)%qtot = 28.3_dp
  lnq(7)%qtot = 0.8_dp
  lnq(8)%qtot = 18.9_dp
  lnq(9)%qtot = 1.0_dp

  allocate(populations(size(cluster_set)))

  ! Call the function
  call calculate_remaining_populations(populations, monomer_populations, lnq, cluster_set, mono)

  ! Check the result
  call check(error, populations(1), expected(1), thr=thr, rel=.true.)
  call check(error, populations(2), expected(2), thr=thr, rel=.true.)
  !call check(error, populations(3), expected(3), thr=thr, rel=.true.)
  call check(error, populations(4), expected(4), thr=thr, rel=.true.)
  call check(error, populations(5), expected(5), thr=thr, rel=.true.)
  call check(error, populations(6), expected(6), thr=thr, rel=.true.)
  call check(error, populations(7), expected(7), thr=thr, rel=.true.)
  call check(error, populations(8), expected(8), thr=thr, rel=.true.)
  call check(error, populations(9), expected(9), thr=thr, rel=.true.)

end subroutine test_calc_remaining_populations_2


end module test_pop_polynomial

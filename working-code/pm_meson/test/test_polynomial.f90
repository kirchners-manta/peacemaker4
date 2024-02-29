module test_polynomial
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use polynomial, only : newton
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

subroutine test_newton(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, 1, 1)
  if (allocated(error)) return
end subroutine test_newton

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  
  call check(error, 1, 0)
  if (allocated(error)) return
end subroutine test_invalid

end module test_polynomial

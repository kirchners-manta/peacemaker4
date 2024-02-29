module test_suite1
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use partition_functions, only : calculate_lnq
  implicit none
  private

  public :: collect_suite1

contains

!> Collect all exported unit tests
subroutine collect_suite1(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test1.1", test_calc_lnq), &
    new_unittest("test1.2", test_invalid, should_fail=.true.) &
    ]

end subroutine collect_suite1

subroutine test_calc_lnq(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, 1, 1)
  if (allocated(error)) return
end subroutine test_calc_lnq

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  
  call check(error, 1, 0)
  if (allocated(error)) return
end subroutine test_invalid

end module test_suite1

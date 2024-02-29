module test_suite2
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use shared_data
  implicit none
  private

  public :: collect_suite2

contains

!> Collect all exported unit tests
subroutine collect_suite2(testsuite)
  !> Collection of tests
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("test2.1", test_valid), &
    new_unittest("test2.2", test_invalid, should_fail=.true.) &
    ]

end subroutine collect_suite2

subroutine test_valid(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, 1, 1)
  if (allocated(error)) return
end subroutine test_valid

subroutine test_invalid(error)
  type(error_type), allocatable, intent(out) :: error
  
  call check(error, 1, 0)
  if (allocated(error)) return
end subroutine test_invalid

end module test_suite2

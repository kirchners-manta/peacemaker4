module test_partition_functions
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds, only : dp
    implicit none
    private
  
    public :: collect_partition_functions
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_partition_functions(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_calc_lnqtrans", test_calc_lnqtrans) &
      !new_unittest("test1.2", test_invalid, should_fail=.true.) &
      ]
  
  end subroutine collect_partition_functions

!----------------------------------------
! Unit test for lnqtrans
!----------------------------------------  
  ! The subroutine calculate_lnqtrans must NOT contain any global variables, 
  ! since it can not be unit tested in that case.
  ! Thus the clusterset is explicitly passed as an argument.
  ! Changes can be seen in partition_functions.f90 - calculate_lnqtrans.
  subroutine test_calc_lnqtrans(error)
    use partition_functions, only : calculate_lnqtrans
    use cluster, only : cluster_t
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2) :: lnq 
    real(dp) :: vol = 2.0
    real(dp) :: temp = 298.0
    real(dp) :: bxv = 0.8
    type(cluster_t), dimension(:), allocatable :: cluster_set
    real(dp) :: v_excl = 0.1

    ! Expected result
    real(dp), dimension(2) :: expected = [59.3344, 62.7883]

    allocate(cluster_set(2))
    cluster_set(1)%mass = 0.001
    cluster_set(2)%mass = 0.01

    call calculate_lnqtrans(lnq, bxv, temp, vol, cluster_set, v_excl)
    write(*,*) lnq
  
    call check(error, lnq(1), expected(1), thr=0.001_dp, rel=.false.)
    call check(error, lnq(2), expected(2), thr=0.001_dp, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqtrans
  
  subroutine test_invalid(error)
    type(error_type), allocatable, intent(out) :: error
    
    call check(error, 1, 0)
    if (allocated(error)) return
  end subroutine test_invalid
  
end module test_partition_functions
  
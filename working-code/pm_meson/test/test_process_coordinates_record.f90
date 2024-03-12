module test_process_coords
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds, only : dp
    implicit none
    private
  
    public :: collect_process_coords
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_process_coords(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_diagonalize_3x3", test_diagonalize_3x3) &
      ]
  
  end subroutine collect_process_coords

!---------------------------------------------
! Unit test for diagonalization of 3x3 matrix
!--------------------------------------------- 
  ! Diagonalization of 3x3 matrix
  subroutine test_diagonalize_3x3(error)
    use cluster, only : diagonalize_3x3
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(3,3) :: a
    real(dp), dimension(3) :: eig
    real(dp), dimension(3) :: expected

    ! Initialize the matrix
    a = reshape([2.0_dp, 2.0_dp, -1.0_dp, 2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp, 4.0_dp], [3,3])

    ! Expected eigenvalues
    expected = [4.64575131_dp, 3.0_dp, -0.64575131_dp]

    call diagonalize_3x3(a, eig)
  
    call check(error, eig(1), expected(1), thr=thr, rel=.false.)
    call check(error, eig(2), expected(2), thr=thr, rel=.false.)
    call check(error, eig(3), expected(3), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_diagonalize_3x3


end module test_process_coords
  
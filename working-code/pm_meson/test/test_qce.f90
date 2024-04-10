module test_qce
    use testdrive, only : new_unittest, unittest_type, error_type, check
    use kinds, only : dp
    !> imports go here
    implicit none
    private
  
    public :: collect_qce
  
  contains
  
  !> Collect all exported unit tests
  subroutine collect_qce(testsuite)
    !> Collection of tests
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
  
    testsuite = [ &
      new_unittest("test_check_convergence", test_check_convergence), &
      new_unittest("test_check_convergence_2", test_check_convergence_2) &
      ]
  
  end subroutine collect_qce
  
  !----------------------------------------------
  ! Unit tests for check_convergence
  !----------------------------------------------

  subroutine test_check_convergence(error)
    use qce, only : check_convergence
    use partition_functions, only : pf_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: gibbs = 0.0_dp
    real(dp) :: temp = 298.15_dp
    real(dp) :: vol = 1.0e-3_dp
    real(dp) :: press = 1.0e5_dp
    real(dp), dimension(4) :: populations = [5.2e23_dp, 1.1e23_dp, 3.8e23_dp, 2.3e23_dp]
    type(pf_t), dimension(4) :: lnq
    logical :: converged    
    real(dp) :: max_dev = 1e-7_dp

    ! Expected result
    real(dp) :: expected = 204069.414303223_dp

    ! Initialize lnq
    lnq%qtot = [5.3_dp, -5.0_dp, 38.2_dp, -1.2_dp]

    call check_convergence(gibbs, temp, vol, press, populations, lnq, converged, max_dev)   
  
    call check(error, gibbs, expected, thr=thr, rel=.false.)
    call check(error, converged, .false.)
    if (allocated(error)) return
  end subroutine test_check_convergence

  subroutine test_check_convergence_2(error)
    use qce, only : check_convergence
    use partition_functions, only : pf_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp) :: gibbs = 4284429.71_dp
    real(dp) :: temp = 512.15_dp
    real(dp) :: vol = 1.0e-5_dp
    real(dp) :: press = 1.0e7_dp
    real(dp), dimension(5) :: populations = [5.8e22_dp, 1.1e20_dp, 3.1e25_dp, 2.3e12_dp, 1.0e15_dp]
    type(pf_t), dimension(5) :: lnq
    logical :: converged    
    real(dp) :: max_dev = 1e-9_dp

    ! Expected result
    real(dp) :: expected = 4284429.70805588_dp

    ! Initialize lnq
    lnq%qtot = [25.3_dp, -51.0_dp, 38.2_dp, -1.2_dp, 111.0_dp]

    call check_convergence(gibbs, temp, vol, press, populations, lnq, converged, max_dev)   
  
    call check(error, gibbs, expected, thr=thr, rel=.false.)
    call check(error, converged, .true.)
    if (allocated(error)) return
  end subroutine test_check_convergence_2
  
  end module test_qce
  
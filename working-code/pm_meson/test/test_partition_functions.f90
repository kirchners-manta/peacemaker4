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
      new_unittest("test_calc_lnqtrans", test_calc_lnqtrans), &
      new_unittest("test_calc_lnqvib", test_calc_lnqvib), &
      new_unittest("test_calc_lnqvib_2", test_calc_lnqvib_2) &
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

    ! Allocate cluster_set
    allocate(cluster_set(2))
    cluster_set(1)%mass = 0.001
    cluster_set(2)%mass = 0.01

    call calculate_lnqtrans(lnq, bxv, temp, vol, cluster_set, v_excl)
  
    call check(error, lnq(1), expected(1), thr=0.001_dp, rel=.false.)
    call check(error, lnq(2), expected(2), thr=0.001_dp, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqtrans

!----------------------------------------
! Unit test for lnqvib
!---------------------------------------- 
  subroutine test_calc_lnqvib(error)
    use partition_functions, only : calculate_lnqvib
    use cluster, only : cluster_t
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2):: lnq
    real(dp) :: temp = 298.0
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Only vibration
    real(dp) :: rotor_cutoff = 0.0

    ! Expected result
    real(dp), dimension(2) :: expected = [0.2584469613283448, 21.647769951871787]

    ! Allocate cluster_set
    allocate(cluster_set(2))
    ! Harmonic oscillator
    cluster_set(1)%frequencies = [100.0, 200.0, 300.0]
    cluster_set(1)%anharmonicity = 0.0
    cluster_set(1)%inertia = [1.0, 1.0, 1.0]

    ! Morse oscillator
    cluster_set(2)%frequencies = [5.0, 648.0, 1000.0, 3555.7] !wavenumbers
    cluster_set(2)%anharmonicity = 0.34
    cluster_set(2)%inertia = [1.0, 1.0, 1.0, 1.0]

    call calculate_lnqvib(lnq, temp, cluster_set, rotor_cutoff)

    call check(error, lnq(1), expected(1), thr=0.00001_dp, rel=.false.)
    call check(error, lnq(2), expected(2), thr=0.00001_dp, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqvib

  subroutine test_calc_lnqvib_2(error)
    use partition_functions, only : calculate_lnqvib
    use cluster, only : cluster_t
    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(3):: lnq
    real(dp) :: temp = 298.0
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Grimme's rotor cutoff
    real(dp) :: rotor_cutoff = 3000.0

    ! Expected result
    real(dp), dimension(3) :: expected = [1.7676272441463634, -3.855446885285766, -3.8012813311299443]

    ! Allocate cluster_set
    allocate(cluster_set(3))
    ! Harmonic oscillator
    cluster_set(1)%frequencies = [100.0, 200.0, 300.0]
    cluster_set(1)%anharmonicity = 0.0
    cluster_set(1)%inertia = [1.0, 1.0, 1.0]
    cluster_set(1)%sigma = 1.0

    ! Harmonic oscillator
    cluster_set(2)%frequencies = [5.0, 648.0, 1000.0, 3555.7] !wavenumbers
    cluster_set(2)%anharmonicity = 0.0
    cluster_set(2)%inertia = [3.2, 5.9, 89.3, 1.0]
    cluster_set(2)%sigma = 1.0

    ! Morse oscillator -> Free rotator approximation is NOT used
    cluster_set(3)%frequencies = [180.0, 270.0, 3550.0]
    cluster_set(3)%anharmonicity = 0.6
    cluster_set(3)%inertia = [1.0, 3.0, 1.0]
    cluster_set(3)%sigma = 1.0

    call calculate_lnqvib(lnq, temp, cluster_set, rotor_cutoff)

    call check(error, lnq(1), expected(1), thr=0.00001_dp, rel=.false.)
    call check(error, lnq(2), expected(2), thr=0.00001_dp, rel=.false.)
    call check(error, lnq(3), expected(3), thr=0.00001_dp, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqvib_2
  
  subroutine test_invalid(error)
    type(error_type), allocatable, intent(out) :: error
    
    call check(error, 1, 0)
    if (allocated(error)) return
  end subroutine test_invalid
  
end module test_partition_functions
  
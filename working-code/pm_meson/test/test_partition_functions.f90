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
      new_unittest("test_calc_lnqvib_2", test_calc_lnqvib_2), &
      new_unittest("test_calc_lnqrot", test_calc_lnqrot), &
      new_unittest("test_calc_lnqelec", test_calc_lnqelec), &
      new_unittest("test_calc_lnqint", test_calc_lnqint), &
      new_unittest("test_calc_dlnqtrans", test_calc_dlnqtrans), &
      new_unittest("test_calc_dlnqtrans_2", test_calc_dlnqtrans_2, should_fail=.true.), &
      new_unittest("test_calc_dlnqvib", test_calc_dlnqvib), &
      new_unittest("test_calc_dlnqvib_2", test_calc_dlnqvib_2), &
      new_unittest("test_calc_dlnqrot", test_calc_dlnqrot) &
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
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2) :: lnq 
    real(dp) :: bxv = 0.8
    real(dp) :: temp = 298.0
    real(dp) :: vol = 2.0
    type(cluster_t), dimension(:), allocatable :: cluster_set
    real(dp) :: v_excl = 0.1

    ! Expected result
    real(dp), dimension(2) :: expected = [59.33445379941631, 62.78833143890738]

    ! Allocate cluster_set
    allocate(cluster_set(2))
    cluster_set(1)%mass = 0.001
    cluster_set(2)%mass = 0.01

    call calculate_lnqtrans(lnq, bxv, temp, vol, cluster_set, v_excl)
  
    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqtrans

!----------------------------------------
! Unit test for lnqvib
!---------------------------------------- 
  subroutine test_calc_lnqvib(error)
    use partition_functions, only : calculate_lnqvib
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2):: lnq
    real(dp) :: temp = 315.0
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Only vibration
    real(dp) :: rotor_cutoff = 0.0

    ! Expected result
    real(dp), dimension(2) :: expected = [-1.3967910870430646, 23.9041034243484]

    ! Allocate cluster_set
    allocate(cluster_set(2))
    ! Harmonic oscillator
    cluster_set(1)%frequencies = [11.4, 25.6, 68.7, 134.9, 3555.9]
    cluster_set(1)%anharmonicity = 0.0
    cluster_set(1)%inertia = [212.4, 212.4]
    cluster_set(1)%sigma = 2

    ! Morse oscillator
    cluster_set(2)%frequencies = [5.0, 648.0, 1000.0, 3555.7] !wavenumbers
    cluster_set(2)%anharmonicity = 0.34
    cluster_set(2)%inertia = [38.6, 43.9, 112.9]
    cluster_set(2)%sigma = 1

    call calculate_lnqvib(lnq, temp, cluster_set, rotor_cutoff)

    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqvib

  subroutine test_calc_lnqvib_2(error)
    use partition_functions, only : calculate_lnqvib
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

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

    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    call check(error, lnq(3), expected(3), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqvib_2

!----------------------------------------
! Unit test for lnqrot
!---------------------------------------- 
  subroutine test_calc_lnqrot(error)
    use partition_functions, only : calculate_lnqrot
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(3) :: lnq 
    real(dp) :: temp = 298.0

    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Expected result
    real(dp), dimension(3) :: expected = [0.0, 6.130167446075413, 9.17613133433584]

    ! Allocate cluster_set
    allocate(cluster_set(3))
    ! Case 1 : Cluster is an atom
    cluster_set(1)%inertia = [0.0]
    cluster_set(1)%sigma = 1.0
    cluster_set(1)%atom = .true.
    cluster_set(1)%linear = .false.

    ! Case 2 : Cluster is linear
    cluster_set(2)%inertia = [37.4, 37.4]
    cluster_set(2)%sigma = 1.0
    cluster_set(2)%atom = .false.
    cluster_set(2)%linear = .true.

    ! Case 3 : Cluster is nonlinear
    cluster_set(3)%inertia = [21.2, 43.2, 70.0]
    cluster_set(3)%sigma = 2.0
    cluster_set(3)%atom = .false.
    cluster_set(3)%linear = .false.


    call calculate_lnqrot(lnq, temp, cluster_set)
  
    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    call check(error, lnq(3), expected(3), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqrot

!----------------------------------------
! Unit test for lnqelec
!---------------------------------------- 

  subroutine test_calc_lnqelec(error)
    use partition_functions, only : calculate_lnqelec
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2) :: lnq 
    real(dp) :: temp = 318.15
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Expected result
    real(dp), dimension(2) :: expected = [25.793440947562555, 10.653072928364542]

    ! Allocate cluster_set
    allocate(cluster_set(2))
    cluster_set(1)%energy = -68.23
    cluster_set(2)%energy = -28.18
    
    call calculate_lnqelec(lnq, temp, cluster_set)
  
    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqelec

!----------------------------------------
! Unit test for lnqelec
!---------------------------------------- 

  subroutine test_calc_lnqint(error)
    use partition_functions, only : calculate_lnqint
    use cluster, only : cluster_t
    use constants, only : avogadro
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(2) :: lnq 
    real(dp) :: amf = 1.0e-48_dp
    real(dp) :: temp = 450.0
    real(dp) :: vol = 3.5e-3
    real(dp), dimension(3) :: ntot = [0.4*avogadro, 0.2*avogadro, 0.4*avogadro]  
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Expected result
    real(dp), dimension(2) :: expected = [0.27694094093862087, 0.2492468468447588]

    ! Allocate cluster_set
    allocate(cluster_set(2))
    cluster_set(1)%composition = [2, 3, 5]
    cluster_set(2)%composition = [1, 6, 2]
    
    call calculate_lnqint(lnq, amf, temp, vol, cluster_set, ntot)
  
    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_lnqint

!----------------------------------------
! Unit test for dlnqtrans
!----------------------------------------  
  ! Temperature derivative of the natural logarithm of the translational partition function
  subroutine test_calc_dlnqtrans(error)
    use partition_functions, only : calculate_dlnqtrans
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(1) :: dlnq 
    real(dp) :: bxv = 0.8
    real(dp) :: bxv_temp = 0.9
    real(dp) :: temp = 431.0
    real(dp) :: vol = 2.5e-3
    real(dp) :: vexcl = 1e-4

    ! Expected result
    real(dp), dimension(1) :: expected = [-0.033709804222354325]

    call calculate_dlnqtrans(dlnq, bxv, bxv_temp, temp, vol, vexcl)

    call check(error, dlnq(1), expected(1), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_dlnqtrans

  ! What happens if the denominator is zero?
  ! This test should fail
  subroutine test_calc_dlnqtrans_2(error)
    use partition_functions, only : calculate_dlnqtrans
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(1) :: dlnq 
    real(dp) :: bxv = 2.0
    real(dp) :: bxv_temp = 0.9
    real(dp) :: temp = 431.0
    real(dp) :: vol = 2
    real(dp) :: vexcl = 1.0

    ! Expected result
    real(dp), dimension(1) :: expected = [0.0]

    call calculate_dlnqtrans(dlnq, bxv, bxv_temp, temp, vol, vexcl)

    call check(error, dlnq(1), expected(1), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_dlnqtrans_2

  !----------------------------------------
  ! Unit test for dlnqvib
  !----------------------------------------
  subroutine test_calc_dlnqvib(error)
    use partition_functions, only : calculate_dlnqvib
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(1):: dlnq
    real(dp) :: temp = 315.0
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Grimme's rotor cutoff
    real(dp) :: rotor_cutoff = 0.0

    ! Expected result
    real(dp), dimension(2) :: expected = [0.03860909908840071, 0.6617367002746789]

    ! Allocate cluster_set
    allocate(cluster_set(2))

    ! Harmonic oscillator, rotor_cutoff = 0.0
    cluster_set(1)%frequencies = [11.4, 25.6, 68.7, 134.9, 3555.9]
    cluster_set(1)%anharmonicity = 0.0

    ! Morse oscillator, rotor_cutoff = 0.0
    cluster_set(2)%frequencies = [11.4, 25.6, 68.7, 134.9, 3555.9, 10000.0]
    cluster_set(2)%anharmonicity = 3.5

    call calculate_dlnqvib(dlnq, temp, cluster_set, rotor_cutoff)
  
    call check(error, dlnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, dlnq(2), expected(2), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_dlnqvib

  subroutine test_calc_dlnqvib_2(error)
    use partition_functions, only : calculate_dlnqvib
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(1):: dlnq
    real(dp) :: temp = 315.0
    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Grimme's rotor cutoff
    real(dp) :: rotor_cutoff = 2398.0

    ! Expected result
    real(dp), dimension(1) :: expected = [0.02798357483860576]

    ! Allocate cluster_set
    allocate(cluster_set(1))
    cluster_set(1)%frequencies = [11.4, 25.6, 68.7, 134.9, 3555.9]
    cluster_set(1)%anharmonicity = 0.0

    call calculate_dlnqvib(dlnq, temp, cluster_set, rotor_cutoff)
  
    call check(error, dlnq(1), expected(1), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_dlnqvib_2
 
!----------------------------------------
! Unit test for dlnqrot
!---------------------------------------- 
  subroutine test_calc_dlnqrot(error)
    use partition_functions, only : calculate_dlnqrot
    use cluster, only : cluster_t
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(3) :: lnq 
    real(dp) :: temp = 512.0

    type(cluster_t), dimension(:), allocatable :: cluster_set

    ! Expected result
    real(dp), dimension(3) :: expected = [0.0, 0.001953125, 0.0029296875]

    ! Allocate cluster_set
    allocate(cluster_set(3))
    ! Case 1 : Cluster is an atom
    cluster_set(1)%atom = .true.
    cluster_set(1)%linear = .false.

    ! Case 2 : Cluster is linear
    cluster_set(2)%atom = .false.
    cluster_set(2)%linear = .true.

    ! Case 3 : Cluster is nonlinear
    cluster_set(3)%atom = .false.
    cluster_set(3)%linear = .false.


    call calculate_dlnqrot(lnq, temp, cluster_set)
  
    call check(error, lnq(1), expected(1), thr=thr, rel=.false.)
    call check(error, lnq(2), expected(2), thr=thr, rel=.false.)
    call check(error, lnq(3), expected(3), thr=thr, rel=.false.)
    if (allocated(error)) return
  end subroutine test_calc_dlnqrot
  
end module test_partition_functions
  
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
      new_unittest("test_check_convergence_2", test_check_convergence_2), &
      new_unittest("test_downhill_simplex_reorder", test_downhill_simplex_reorder), &
      new_unittest("test_downhill_simplex_reorder_2", test_downhill_simplex_reorder_2), &
      new_unittest("test_downhill_simplex_reorder_3", test_downhill_simplex_reorder_3), &
      new_unittest("test_downhill_simplex_reorder_4", test_downhill_simplex_reorder_4), &
      new_unittest("test_downhill_simplex_middle", test_downhill_simplex_middle), &
      new_unittest("test_downhill_simplex_middle_2", test_downhill_simplex_middle_2), &
      new_unittest("test_downhill_simplex_middle_3", test_downhill_simplex_middle_3), &
      new_unittest("test_downhill_simplex_reflection", test_downhill_simplex_reflection) &
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

  !----------------------------------------------
  ! Unit tests for downhill_simplex reorder
  !----------------------------------------------
  ! The reorder subroutine reorders the rows of the simplex based on the function values
  ! which are stored in the last column of the simplex. The row with the lowest function
  ! value is moved to the first row, the row with the second lowest function value is moved
  ! to the second row, and so on. 

  ! One parameter is optimized
  ! n = 2
  subroutine test_downhill_simplex_reorder(error)
    use qce, only : reorder

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), allocatable, dimension(:,:) :: simplex 

    ! Expected result
    real(dp), dimension(2, 2) :: expected = reshape([1.5_dp, 1.0_dp, 2.0_dp, 4.0_dp], [2, 2])
    ! [1.5, 2.0]
    ! [1.0, 4.0]

    ! Allocate and initialize simplex
    allocate(simplex(2, 2))
    simplex = reshape([1.0_dp, 1.5_dp, 4.0_dp, 2.0_dp], [2, 2])
    ! [1.0, 4.0]
    ! [1.5, 2.0]

    ! Call reorder
    call reorder(simplex)

    ! Check each element
    do i = 1, 2
      do j = 1, 2
        call check(error, simplex(i, j), expected(i, j))
      end do
    end do    

  end subroutine test_downhill_simplex_reorder

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_reorder_2(error)
    use qce, only : reorder

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), allocatable, dimension(:,:) :: simplex 

    ! Expected result
    real(dp), dimension(3, 3) :: expected = reshape([2.0_dp, 3.0_dp, 1.0_dp, &
                                                    5.2_dp, 0.8_dp, -2.0_dp, &
                                                    -2.3_dp, 1.1_dp, 91.2_dp], [3, 3])
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]
    ! [1.0, -2.0, 91.2]
    

    ! Allocate and initialize simplex
    allocate(simplex(3, 3))
    simplex = reshape([1.0_dp, 2.0_dp, 3.0_dp, &
                      -2.0_dp, 5.2_dp, 0.8_dp, &
                       91.2_dp, -2.3_dp, 1.1_dp], [3, 3])
    ! [1.0, -2.0, 91.2]
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]

    ! Call reorder
    call reorder(simplex)

    ! Check each element
    do i = 1, 3
      do j = 1, 3
        call check(error, simplex(i, j), expected(i, j))
      end do
    end do    

  end subroutine test_downhill_simplex_reorder_2

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_reorder_3(error)
    use qce, only : reorder

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), allocatable, dimension(:,:) :: simplex 

    ! Expected result
    real(dp), dimension(4, 4) :: expected = reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, &
                                                    1.1_dp, 5.2_dp, 0.8_dp, -2.0_dp, &
                                                    0.8_dp, -2.3_dp, 1.1_dp, 91.2_dp, &
                                                    1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [4, 4])     
    ! [4.0, 1.1, 0.8, 1.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [1.0, -2.0, 91.2, 4.0]                                      
    
    ! Allocate and initialize simplex
    allocate(simplex(4, 4))
    simplex = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                      -2.0_dp, 5.2_dp, 0.8_dp, 1.1_dp, &
                       91.2_dp, -2.3_dp, 1.1_dp, 0.8_dp, &
                       4.0_dp, 3.0_dp, 3.0_dp, 1.0_dp], [4, 4])
    ! [1.0, -2.0, 91.2, 4.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [4.0, 1.1, 0.8, 1.0]

    ! Call reorder
    call reorder(simplex)

    ! Check each element
    do i = 1, 4
      do j = 1, 4
        call check(error, simplex(i, j), expected(i, j))
      end do
    end do    

  end subroutine test_downhill_simplex_reorder_3

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_reorder_4(error)
    use qce, only : reorder

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), allocatable, dimension(:,:) :: simplex 

    ! Expected result
    real(dp), dimension(5, 5) :: expected = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                                                    -2.0_dp, 5.2_dp, 0.8_dp, 1.1_dp, 2.0_dp, &
                                                     91.2_dp, -2.3_dp, 1.1_dp, 0.8_dp, 3.0_dp, &
                                                     4.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, &
                                                     1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [5, 5])                                    
    
    ! Allocate and initialize simplex
    allocate(simplex(5, 5))
    simplex = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                      -2.0_dp, 5.2_dp, 0.8_dp, 1.1_dp, 2.0_dp, &
                       91.2_dp, -2.3_dp, 1.1_dp, 0.8_dp, 3.0_dp, &
                       4.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, &
                       1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [5, 5])
    ! [1.0, -2.0, 91.2, 4.0, 1.0]
    ! [2.0, 5.2, -2.3, 1.0, 2.0]
    ! [3.0, 0.8, 1.1, 3.0, 3.0]
    ! [4.0, 1.1, 0.8, 3.0, 4.0]
    ! [5.0, 2.0, 3.0, 4.0, 5.0]

    ! Call reorder
    call reorder(simplex)

    ! Check each element
    do i = 1, 5
      do j = 1, 5
        call check(error, simplex(i, j), expected(i, j))
      end do
    end do    

  end subroutine test_downhill_simplex_reorder_4

  !----------------------------------------------
  ! Unit tests for downhill_simplex middle
  !----------------------------------------------
  ! Calculates midpoint m of n-1 best points x1, ..., xn-1

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_middle(error)
    use qce, only : middle

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(3) :: m
    logical :: skip_qce = .true.
    real(dp), dimension(3, 3) :: simplex = reshape([2.0_dp, 3.0_dp, 1.0_dp, &
                                                    5.2_dp, 0.8_dp, -2.0_dp, &
                                                    -2.3_dp, 1.1_dp, 91.2_dp], [3, 3])
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]
    ! [1.0, -2.0, 91.2]                                                    

    ! Expected result
    real(dp), dimension(3) :: expected = [2.5_dp, 3.0_dp, 0.0_dp]
    ! The third entry is the value of the point (2.5, 3.0). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call middle
    call middle(simplex, m, skip_qce)

    ! Check each element
    do i = 1, 2
      call check(error, m(i), expected(i))
    end do    

  end subroutine test_downhill_simplex_middle

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_middle_2(error)
    use qce, only : middle
    real(dp) :: thr = 1.0e-3_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(4) :: m
    logical :: skip_qce = .true.
    real(dp), dimension(4, 4) :: simplex = reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, &
                                                    1.1_dp, 5.2_dp, 0.8_dp, -2.0_dp, &
                                                    0.8_dp, -2.3_dp, 1.1_dp, 91.2_dp, &
                                                    1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [4, 4])
    ! [4.0, 1.1, 0.8, 1.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [1.0, -2.0, 91.2, 4.0]                      

    ! Expected result
    real(dp), dimension(4) :: expected = [3.0_dp, 2.3666_dp, -0.1333_dp, 0.0_dp]
    ! The third entry is the value of the point (3.0, 2.3666, -0.1333). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call middle
    call middle(simplex, m, skip_qce)

    ! Check each element
    do i = 1, 3
      call check(error, m(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_middle_2

  ! Four parameter are optimized
  ! n = 5
  subroutine test_downhill_simplex_middle_3(error)
    use qce, only : middle
    real(dp) :: thr = 1.0e-1_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(5) :: m
    logical :: skip_qce = .true.
    real(dp), dimension(5, 5) :: simplex = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                                                    -2.0_dp, 5.2_dp, 0.8_dp, 1.1_dp, 2.0_dp, &
                                                     91.2_dp, -2.3_dp, 1.1_dp, 0.8_dp, 3.0_dp, &
                                                     4.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, &
                                                     1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [5, 5])
    ! [1.0, -2.0, 91.2, 4.0, 1.0]
    ! [2.0, 5.2, -2.3, 1.0, 2.0]
    ! [3.0, 0.8, 1.1, 3.0, 3.0]
    ! [4.0, 1.1, 0.8, 3.0, 4.0]
    ! [5.0, 2.0, 3.0, 4.0, 5.0]             

    ! Expected result
    real(dp), dimension(5) :: expected = [2.5_dp, 1.275_dp, 22.7_dp, 2.75_dp, 0.0_dp]
    ! The third entry is the value of the point (2.5, 1.275, 22.7, 2.75). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call middle
    call middle(simplex, m, skip_qce)

    ! Check each element
    do i = 1, 3
      call check(error, m(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_middle_3

  !---------------------------------------------
  ! Unit tests for downhill_simplex reflection
  !---------------------------------------------
  ! Reflects the worst point xn at midpoint m to form point r

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_reflection(error)
    use qce, only : reflection

    real(dp) :: thr = 1.0e-1_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(3) :: r
    logical :: skip_qce = .true.
    real(dp), dimension(3) :: m = [2.5_dp, 3.0_dp, 0.0_dp]
    real(dp), dimension(3, 3) :: simplex = reshape([2.0_dp, 3.0_dp, 1.0_dp, &
                                                    5.2_dp, 0.8_dp, -2.0_dp, &
                                                    -2.3_dp, 1.1_dp, 91.2_dp], [3, 3])
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]
    ! [1.0, -2.0, 91.2]   
    

    ! Expected result
    real(dp), dimension(3) :: expected = [4.0_dp, 8.0_dp, 0.0_dp]
    ! The third entry is the value of the point (4, 8). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call reflect
    call reflection(simplex, m, r, skip_qce)

    ! Check each element
    do i = 1, 3
      call check(error, r(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_reflection

  end module test_qce
  
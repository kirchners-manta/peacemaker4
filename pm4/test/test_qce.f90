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
      new_unittest("test_downhill_simplex_reflection", test_downhill_simplex_reflection), &
      new_unittest("test_downhill_simplex_reflection_2", test_downhill_simplex_reflection_2), &
      new_unittest("test_downhill_simplex_reflection_3", test_downhill_simplex_reflection_3), &
      new_unittest("test_downhill_simplex_expansion", test_downhill_simplex_expansion), &
      new_unittest("test_downhill_simplex_expansion_2", test_downhill_simplex_expansion_2), &
      new_unittest("test_downhill_simplex_expansion_3", test_downhill_simplex_expansion_3), &
      new_unittest("test_downhill_simplex_contraction", test_downhill_simplex_contraction), &
      new_unittest("test_downhill_simplex_contraction_2", test_downhill_simplex_contraction_2), &
      new_unittest("test_downhill_simplex_contraction_3", test_downhill_simplex_contraction_3), &
      new_unittest("test_downhill_simplex_compression", test_downhill_simplex_compression), &
      new_unittest("test_downhill_simplex_compression_2", test_downhill_simplex_compression_2), &
      new_unittest("test_downhill_simplex_compression_3", test_downhill_simplex_compression_3) &
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
    call middle(simplex, m)

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
    call middle(simplex, m)

    ! Check each element
    do i = 1, 3
      call check(error, m(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_middle_2

  ! Four parameter are optimized
  ! n = 5
  subroutine test_downhill_simplex_middle_3(error)
    use qce, only : middle
    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(5) :: m
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
    call middle(simplex, m)

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

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(3) :: r
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
    call reflection(simplex, m, r)

    ! Check each element
    do i = 1, 2
      call check(error, r(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_reflection

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_reflection_2(error)
    use qce, only : reflection

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(4) :: r
    real(dp), dimension(4) :: m = [3.0_dp, 2.3666_dp, -0.1333_dp, 0.0_dp]
    real(dp), dimension(4, 4) :: simplex = reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, &
                                                    1.1_dp, 5.2_dp, 0.8_dp, -2.0_dp, &
                                                    0.8_dp, -2.3_dp, 1.1_dp, 91.2_dp, &
                                                    1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [4, 4])
    ! [4.0, 1.1, 0.8, 1.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [1.0, -2.0, 91.2, 4.0]       

    ! Expected result
    real(dp), dimension(4) :: expected = [5.0_dp, 6.7332_dp, -91.4666_dp, 0.0_dp]
    ! The third entry is the value of the point (5.0, 6.7332, 91.4666). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call reflect
    call reflection(simplex, m, r)

    ! Check each element
    do i = 1, 3
      call check(error, r(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_reflection_2

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_reflection_3(error)
    use qce, only : reflection

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(5) :: r
    real(dp), dimension(5) :: m = [2.5_dp, 1.275_dp, 22.7_dp, 2.75_dp, 0.0_dp]
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
    real(dp), dimension(5) :: expected = [0.0_dp, 0.55_dp, 42.4_dp, 1.5_dp, 0.0_dp]
    ! The third entry is the value of the point (0.0, 0.55, 42.4, 1.5). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call reflect
    call reflection(simplex, m, r)

    ! Check each element
    do i = 1, 4
      call check(error, r(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_reflection_3

  !---------------------------------------------
  ! Unit tests for downhill_simplex expansion
  !---------------------------------------------
  ! Expands the worst point xn at midpoint m to form point e

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_expansion(error)
    use qce, only : expansion

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(3) :: e
    real(dp), dimension(3) :: m = [2.5_dp, 3.0_dp, 0.0_dp]
    real(dp), dimension(3, 3) :: simplex = reshape([2.0_dp, 3.0_dp, 1.0_dp, &
                                                    5.2_dp, 0.8_dp, -2.0_dp, &
                                                    -2.3_dp, 1.1_dp, 91.2_dp], [3, 3])
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]
    ! [1.0, -2.0, 91.2]   

    ! Expected result
    real(dp), dimension(3) :: expected = [5.5_dp, 13.0_dp, 0.0_dp]
    ! The third entry is the value of the point (3.0, 3.5). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call expand
    call expansion(simplex, m, e)

    ! Check each element
    do i = 1, 2
      call check(error, e(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_expansion

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_expansion_2(error)
    use qce, only : expansion

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(4) :: e
    real(dp), dimension(4) :: m = [3.0_dp, 2.3666_dp, -0.1333_dp, 0.0_dp]
    real(dp), dimension(4, 4) :: simplex = reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, &
                                                    1.1_dp, 5.2_dp, 0.8_dp, -2.0_dp, &
                                                    0.8_dp, -2.3_dp, 1.1_dp, 91.2_dp, &
                                                    1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [4, 4])
    ! [4.0, 1.1, 0.8, 1.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [1.0, -2.0, 91.2, 4.0]       

    ! Expected result
    real(dp), dimension(4) :: expected = [7.0_dp, 11.0998_dp, -182.7999_dp, 0.0_dp]
    ! The third entry is the value of the point (7.0, 11.0998, -182.799). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call expand
    call expansion(simplex, m, e)

    ! Check each element
    do i = 1, 3
      call check(error, e(i), expected(i), thr=thr, rel=.false.)
    end do    

  end subroutine test_downhill_simplex_expansion_2

  ! Four parameter are optimized
  ! n = 5
  subroutine test_downhill_simplex_expansion_3(error)
    use qce, only : expansion

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(5) :: e
    real(dp), dimension(5) :: m = [2.5_dp, 1.275_dp, 22.7_dp, 2.75_dp, 0.0_dp]
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
    real(dp), dimension(5) :: expected = [-2.5_dp, -0.175_dp, 62.1_dp, 0.25_dp, 0.0_dp]
    ! The third entry is the value of the point (-2.5, -0.175, 62.1, 0.25). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call expand
    call expansion(simplex, m, e)

    ! Check each element
    do i = 1, 4
      call check(error, e(i), expected(i), thr=thr, rel=.false.)
    end do

  end subroutine test_downhill_simplex_expansion_3

  !---------------------------------------------
  ! Unit tests for downhill_simplex contraction
  !---------------------------------------------
  ! Contracts the better point h of r or xn with the midpoint m to form point c

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_contraction(error)
    use qce, only : contraction

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(3) :: c
    logical :: skip_qce = .true.
    real(dp), dimension(3) :: m = [2.5_dp, 3.0_dp, 0.0_dp]
    real(dp), dimension(3) :: h = [1.0_dp, -2.0_dp, 0.0_dp] 

    ! Expected result
    real(dp), dimension(3) :: expected = [1.75_dp, 0.5_dp, 0.0_dp]
    ! The third entry is the value of the point (1.75, 0.5). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call contract
    call contraction(m, h, c)

    ! Check each element
    do i = 1, 2
      call check(error, c(i), expected(i), thr=thr, rel=.false.)
    end do

  end subroutine test_downhill_simplex_contraction

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_contraction_2(error)
    use qce, only : contraction

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(4) :: c
    real(dp), dimension(4) :: m = [3.0_dp, 2.3666_dp, -0.1333_dp, 0.0_dp]
    real(dp), dimension(4) :: h = [4.0_dp, 1.1_dp, 0.8_dp, 0.0_dp] 

    ! Expected result
    real(dp), dimension(4) :: expected = [3.5_dp, 1.7333_dp, 0.33335_dp, 0.0_dp]
    ! The third entry is the value of the point (3.5, 1.7333, 0.33335). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call contract
    call contraction(m, h, c)

    ! Check each element
    do i = 1, 3
      call check(error, c(i), expected(i), thr=thr, rel=.false.)
    end do

  end subroutine test_downhill_simplex_contraction_2

  ! Four parameter are optimized
  ! n = 5
  subroutine test_downhill_simplex_contraction_3(error)
    use qce, only : contraction

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i
    real(dp), dimension(5) :: c
    real(dp), dimension(5) :: m = [2.5_dp, 1.275_dp, 22.7_dp, 2.75_dp, 0.0_dp]
    real(dp), dimension(5) :: h = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp] 

    ! Expected result
    real(dp), dimension(5) :: expected = [1.75_dp, 1.6375_dp, 12.85_dp, 3.375_dp, 0.0_dp]
    ! The third entry is the value of the point (1.74, 1.6375, 12.85, 3.375). This is not calculated
    ! for unit tests, so it is set to 0.0.

    ! Call contract
    call contraction(m, h, c)

    ! Check each element
    do i = 1, 4
      call check(error, c(i), expected(i), thr=thr, rel=.false.)
    end do

  end subroutine test_downhill_simplex_contraction_3

  !---------------------------------------------
  ! Unit tests for downhill_simplex compression
  !---------------------------------------------
  ! Compresses the simplex

  ! Two parameter are optimized
  ! n = 3
  subroutine test_downhill_simplex_compression(error)
    use qce, only : compression

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), dimension(3, 3) :: simplex = reshape([2.0_dp, 3.0_dp, 1.0_dp, &
                                                    5.2_dp, 0.8_dp, -2.0_dp, &
                                                    -2.3_dp, 1.1_dp, 91.2_dp], [3, 3])
    ! [2.0, 5.2, -2.3]
    ! [3.0, 0.8, 1.1]
    ! [1.0, -2.0, 91.2]

    ! Expected result
    real(dp), dimension(3, 3) :: expected = reshape([2.0_dp, 2.5_dp, 1.5_dp, &
                                                     5.2_dp, 3.0_dp, 1.6_dp, &
                                                     -2.3_dp, 1.1_dp, 0.0_dp], [3, 3])

    ! Call compress
    call compression(simplex)

    ! Check each element
    do i = 1, 3
      do j = 1, 2
        call check(error, simplex(i, j), expected(i, j), thr=thr, rel=.false.)
      end do
    end do    

  end subroutine test_downhill_simplex_compression

  ! Three parameter are optimized
  ! n = 4
  subroutine test_downhill_simplex_compression_2(error)
    use qce, only : compression

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), dimension(4, 4) :: simplex = reshape([4.0_dp, 2.0_dp, 3.0_dp, 1.0_dp, &
                                                    1.1_dp, 5.2_dp, 0.8_dp, -2.0_dp, &
                                                    0.8_dp, -2.3_dp, 1.1_dp, 91.2_dp, &
                                                    1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp], [4, 4])
    ! [4.0, 1.1, 0.8, 1.0]
    ! [2.0, 5.2, -2.3, 3.0]
    ! [3.0, 0.8, 1.1, 3.0]
    ! [1.0, -2.0, 91.2, 4.0]    

    ! Expected result
    real(dp), dimension(4, 4) :: expected = reshape([4.0_dp, 3.0_dp, 3.5_dp, 2.5_dp, &
                                                     1.1_dp, 3.15_dp, 0.95_dp, -0.45_dp, &
                                                     0.8_dp, -0.75_dp, 0.95_dp, 46.0_dp, &
                                                     1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [4, 4])

    ! Call compress
    call compression(simplex)

    ! Check each element
    do i = 1, 4
      do j = 1, 3
        call check(error, simplex(i,j), expected(i,j), thr=thr, rel=.false.)
      end do
    end do

  end subroutine test_downhill_simplex_compression_2

  ! Four parameter are optimized
  ! n = 5
  subroutine test_downhill_simplex_compression_3(error)
    use qce, only : compression

    real(dp) :: thr = 1.0e-2_dp

    ! Declare error and the subroutine check if not yet declared
    type(error_type), allocatable, intent(out) :: error
    integer :: i, j
    real(dp), dimension(5, 5) :: simplex = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, &
                                                    -2.0_dp, 5.2_dp, 0.8_dp, 1.1_dp, 2.0_dp, &
                                                     91.2_dp, -2.3_dp, 1.1_dp, 0.8_dp, 3.0_dp, &
                                                     4.0_dp, 1.0_dp, 3.0_dp, 3.0_dp, 4.0_dp, &
                                                     1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], [5, 5])

    ! [1.0, -2.0, 91.2, 4.0, 1.0]
    ! [2.0, 5.2, -2.3, 1.0, 2.0]
    ! [3.0, 0.8, 1.1, 3.0, 3.0]
    ! [4.0, 1.1, 0.8, 3.0, 4.0]

    ! Expected result
    real(dp), dimension(5, 5) :: expected = reshape([1.0_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp, &
                                                     -2.0_dp, 1.6_dp, -0.6_dp, -0.45_dp, 0.0_dp, &
                                                     91.2_dp, 44.45_dp, 46.15_dp, 46.0_dp, 47.1_dp, &
                                                     4.0_dp, 2.5_dp, 3.5_dp, 3.5_dp, 4.0_dp, &
                                                     0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], [5, 5])

    ! Call compress
    call compression(simplex)

    ! Check each element
    do i = 1, 5
      do j = 1, 4
        call check(error, simplex(i,j), expected(i,j), thr=thr, rel=.false.)
      end do
    end do

  end subroutine test_downhill_simplex_compression_3

  end module test_qce
  
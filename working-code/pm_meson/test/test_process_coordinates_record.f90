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
      new_unittest("test_diagonalize_3x3", test_diagonalize_3x3), &
      new_unittest("test_center_of_mass", test_center_of_mass) &
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

!---------------------------------------------
! Unit test for center_of_mass
!--------------------------------------------- 
  ! Diagonalization of 3x3 matrix
  subroutine test_center_of_mass(error)
    use cluster, only : center_of_mass
    !> Precision for the tests
    real(dp) :: thr = 1.0e-5_dp

    ! Arguments
    type(error_type), allocatable, intent(out) :: error
    real(dp), dimension(3) :: expected_com
    real(dp), dimension(20,3) :: expected_xyz
    integer :: nr_atoms
    real(dp), dimension(3) :: com
    real(dp), dimension(20) :: mass
    real(dp), dimension(20,3) :: xyz
    integer :: i

    ! Initialize
    nr_atoms = 20
    mass = [12.01, 35.45, 35.45, 1.008, 35.45, 12.01, 1.008, 16.0, 1.008, 1.008, &
              1.008, 16.0, 1.008, 1.008, 16.0, 1.008, 1.008, 16.0, 1.008, 1.008]
    xyz = reshape([-1.16898944033827_dp, -2.34784356765353_dp, -0.65822727470503_dp, &   
                   -0.29281244953045_dp, -1.88998652011726_dp,  2.76478224676894_dp, &   
                    3.04017750865408_dp,  3.01513896067189_dp,  1.71213083048118_dp, &   
                    3.37938957253550_dp,  2.47344440334455_dp,  1.38749152659684_dp, &   
                    1.75428367438544_dp,  1.39549826682973_dp,  1.34486902387060_dp, &   
                    1.71731567397097_dp,  0.47187877834489_dp,  2.14004917691883_dp, &   
                    2.62134584522111_dp,  1.40266719753673_dp, -0.14805428493990_dp, &   
                   -1.28928508067893_dp, -0.69129854923607_dp, -0.09624032526792_dp, &
                    1.46443581796416_dp, -2.07459581216681_dp, -2.64738864366591_dp, &
                   -0.70150365130769_dp, -2.26384057006596_dp, -2.40571568513168_dp, &
                   -0.37364082813805_dp,  0.47508052337209_dp,  0.63637236685189_dp, &
                    1.34971325441160_dp,  2.68669110647591_dp,  2.30907836219124_dp, &
                    3.01687474618809_dp,  1.26284234041995_dp,  0.56485487037962_dp, &
                    0.82485404829678_dp,  0.32745185542377_dp,  0.99930307590600_dp, &
                   -1.27452373497366_dp,  0.98689811267938_dp,  0.20551396958876_dp, &
                   -0.19885824933856_dp,  0.69127857857294_dp, -0.02078584488278_dp, &
                   -0.43213510816834_dp, -1.03494010845601_dp,  0.72736968185670_dp, &
                    1.86216790536171_dp,  2.73283384502543_dp,  1.39534117990936_dp, &
                    0.41214803245629_dp, -0.41677149677762_dp,  0.18984991489869_dp, &
                   -1.72480591317545_dp, -1.23969806978062_dp,  -2.15680622354843_dp], [20,3])
                   

    ! Expected moments of inertia
    expected_com = [-0.04072005, 0.07544383, 0.04374451]
    expected_xyz = reshape([-1.12826939_dp, -2.30712351_dp, -0.61750722_dp, -0.25209239_dp, -1.84926647_dp, &
                             2.80550230_dp,  3.08089756_dp,  3.05585902_dp,  1.75285089_dp,  3.42010963_dp, &
                             2.51416446_dp,  1.42821158_dp,  1.79500373_dp,  1.43621832_dp,  1.38558908_dp, &
                             1.75803573_dp,  0.51259883_dp,  2.18076923_dp,  2.66206590_dp,  1.44338725_dp, &
                            -0.22349812_dp, -1.36472891_dp, -0.76674238_dp, -0.17168416_dp,  1.38899199_dp, &
                            -2.15003964_dp, -2.72283247_dp, -0.77694748_dp, -2.33928440_dp, -2.48115952_dp, &
                            -0.44908466_dp,  0.39963669_dp,  0.56092854_dp,  1.27426942_dp,  2.61124728_dp, &
                             2.23363453_dp,  2.94143092_dp,  1.18739851_dp,  0.48941104_dp,  0.74941022_dp, &
                             0.28370734_dp,  0.95555856_dp, -1.31826825_dp,  0.94315360_dp,  0.16176946_dp, &
                            -0.24260276_dp,  0.64753407_dp, -0.06453036_dp, -0.47587962_dp, -1.07868462_dp, &
                             0.68362517_dp,  1.81842339_dp,  2.68908933_dp,  1.35159667_dp,  0.36840352_dp, &
                            -0.46051601_dp,  0.14610540_dp, -1.76855043_dp, -1.28344258_dp, -2.20055074_dp], [20,3])

    call center_of_mass(nr_atoms, com, mass, xyz)
  
    call check(error, com(1), expected_com(1), thr=thr, rel=.false.)
    call check(error, com(2), expected_com(2), thr=thr, rel=.false.)
    call check(error, com(3), expected_com(3), thr=thr, rel=.false.)

    do i=1,20
      call check(error, xyz(i,1), expected_xyz(i,1), thr=thr, rel=.false.)
      call check(error, xyz(i,2), expected_xyz(i,2), thr=thr, rel=.false.)
      call check(error, xyz(i,3), expected_xyz(i,3), thr=thr, rel=.false.)
    end do  
    

    if (allocated(error)) return
  end subroutine test_center_of_mass


end module test_process_coords
  
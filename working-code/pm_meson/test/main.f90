program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  !> > > > > >   BLOCK 1   > > > > > > > > > > >
  use test_pop_polynomial, only : collect_pop_polynomial
  use test_vol_polynomial, only : collect_vol_polynomial
  use test_partition_functions, only : collect_partition_functions
  use test_process_coords, only : collect_process_coords
  use test_qce, only : collect_qce
  use test_thermo, only : collect_thermo
  use test_auxiliary, only : collect_auxiliary
  !> ^- Hier use  test_<name>, only : collect_<name>  einfügen

  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    !> > > > > >   BLOCK 2   > > > > > > > > > >
    new_testsuite("pop_polynomial", collect_pop_polynomial), &
    new_testsuite("vol_polynomial", collect_vol_polynomial), &
    new_testsuite("partition_functions", collect_partition_functions), &
    new_testsuite("process_coordinates_record", collect_process_coords), &
    new_testsuite("qce", collect_qce), &
    new_testsuite("thermo", collect_thermo), &
    new_testsuite("auxiliary", collect_auxiliary) &
    !> ^- Hier  new_testsuite("<name>", collect_<name>), &  einfügen
    ]

  do is = 1, size(testsuites)
    write(error_unit, fmt) "Testing:", testsuites(is)%name
    call run_testsuite(testsuites(is)%collect, error_unit, stat)
  end do

  if (stat > 0) then
    write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
    error stop
  end if

end program tester

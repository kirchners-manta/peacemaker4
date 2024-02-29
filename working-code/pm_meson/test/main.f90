program tester
  use, intrinsic :: iso_fortran_env, only : error_unit
  use testdrive, only : run_testsuite, new_testsuite, testsuite_type
  !> > > > > >   BLOCK 1   > > > > > > > > > > >
  use test_polynomial, only : collect_polynomial
  use test_suite2, only : collect_suite2
  !> ^- Hier use  test_<name>, only : collect_<name>  einfügen

  implicit none
  integer :: stat, is
  type(testsuite_type), allocatable :: testsuites(:)
  character(len=*), parameter :: fmt = '("#", *(1x, a))'

  stat = 0

  testsuites = [ &
    !> > > > > >   BLOCK 2   > > > > > > > > > >
    new_testsuite("polynomial", collect_polynomial), &
    new_testsuite("suite2", collect_suite2) &
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

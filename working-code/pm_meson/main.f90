!=========================================================================================
! Peacemaker -- A Quantum Cluster Equilibrium Code.
!
! Copyright 2004-2006 Barbara Kirchner, University of Bonn
! Copyright 2007-2012 Barbara Kirchner, University of Leipzig
! Copyright 2013-2022 Barbara Kirchner, University of Bonn
!
! The Peacemaker team:
!
!   * Michael von Domaros
!   * Eva Perlt
!   * Johannes Ingenmey
!   * Marc Bruessel
!   * Christian Spickermann
!   * Sebastian Lehmann
!
!   * Barbara Kirchner
!
! This file is part of Peacemaker.
!
! Peacemaker is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Peacemaker is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Peacemaker.  If not, see <http://www.gnu.org/licenses/>
!=========================================================================================
! This is Peacemaker!
program main
    use kinds
    use lengths
    use info
    use atomic_data
    use input
    use cluster
    use qce
    use tomlf, only: toml_table
    implicit none
    !=====================================================================================
    ! Input and clusterset file names and configurations.
    character(path_len):: input_file
    character(path_len):: clusterset_file
    !=====================================================================================
    ! Time measurement
    integer:: pmk_start
    integer:: pmk_stop
    integer:: count_rate
    !=====================================================================================
    ! Check command line argument count.
    if (command_argument_count() /= 2) then
        call print_usage_info()
        stop
    end if
    !=====================================================================================
    ! Hello, user!
    call print_welcome_info()
    !=====================================================================================
    ! Initialize periodic table.
    call periodic_table%init()
    !=====================================================================================
    ! Parse input and clusterset file.
    call get_command_argument(1, input_file)
    call get_command_argument(2, clusterset_file)
    !=====================================================================================
    ! Process input, perform sanity checks and print the processed input.
    call process_input(table, input_file)
    call check_input()
    call print_input()
    !=====================================================================================
    ! Setup clusterset, perform sanity check and print the processed input.
    call process_clusterset(cluster_table, clusterset_file)
    call check_clusterset()
    call print_clusterset()
    !=====================================================================================
    ! Warn user about unread input/clusterset entries.
    !call input_cfg%check()
    !call clusterset_cfg%check()
    !=====================================================================================
    ! Start time measurement.
    call system_clock(pmk_start)
    !=====================================================================================
    ! Perform QCE calculations.
    call qce_prepare()

    ! Stop time measurement and report.
    call system_clock(pmk_stop, count_rate)
    write(*,'(4X,A,1X,G0.2,1X,A)') &
        "Elapsed time - prepare:", real(pmk_stop-pmk_start) / real(count_rate), "seconds"
    write(*,*)
    call system_clock(pmk_start)

    call qce_start()

    ! Stop time measurement and report.
    call system_clock(pmk_stop, count_rate)
    write(*,'(4X,A,1X,G0.2,1X,A)') &
        "Elapsed time - start:", real(pmk_stop-pmk_start) / real(count_rate), "seconds"
    write(*,*)
    call system_clock(pmk_start)

    call qce_finalize()
    
    ! Stop time measurement and report.
    call system_clock(pmk_stop, count_rate)
    write(*,'(4X,A,1X,G0.2,1X,A)') &
        "Elapsed time - finalize:", real(pmk_stop-pmk_start) / real(count_rate), "seconds"
    write(*,*)

    !=====================================================================================
    ! Print citations.
    call print_citation_info()
    !=====================================================================================
end program main
!=========================================================================================

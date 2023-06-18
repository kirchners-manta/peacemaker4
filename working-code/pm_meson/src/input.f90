!=========================================================================================
! Peacemaker -- A Quantum Cluster Equilibrium Code.
! 
! Copyright 2004-2006 Barbara Kirchner, University of Bonn
! Copyright 2007-2012 Barbara Kirchner, University of Leipzig
! Copyright 2013-2022 Barbara Kirchner, University of Bonn
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
! This module provides the input_t data type and its only instantiation pmk_input, which
! contains the processed input from the input configuration file. The module also provides
! the subroutines process_input(), check_input(), and print_input(), which process the
! input, perform sanity checks and print the processed input, respectively.
!
! How to add a new keyword to the input and get data into peacemaker:
!    1. Add a variable corresponding to the keyword to the input_t data type.
!    2. Set a default value inside process_input(). If the variable is allocatable,
!       allocate it to a default size. Make sure that default values are consistent with
!       other defaults.
!    3. Proccess the keyword by adding a new subroutine in analogy to the other
!       subroutines.
!    4. Perform sanity checks on the processed data inside check_input()
!    5. Print the processed input inside print_input().
!=========================================================================================
module input
    use kinds
    use lengths
    use iso_varying_string
    use config
    use error
    use auxiliary
    use constants

    ! for readin toml
    !use reader, only : read_data, convert_helpers
    use tomlf, only : toml_table, toml_parse, toml_error, toml_array, get_value, len, toml_key

    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: pmk_input
    public :: process_input
    public :: check_input
    public :: print_input
    public :: input_data
    public :: table
    public :: toml_parse
    !=====================================================================================
    ! The input_t data type.
    ! Input type for toml input
    type :: input_data
        ! Input read fromt the [system] section.
        integer :: components

        ! Input read from the [ensemble] section.
        type(range_t) :: temperature
        real(dp) :: pressure
        real(dp), dimension(:), allocatable :: monomer_amounts

        ! Input read from the [qce] section.
        type(range_t) :: amf, bxv
        type(range_t) :: amf_temp, bxv_temp
        real(dp) :: max_deviation, volume_damping_factor, rotor_cutoff
        integer :: qce_iterations, newton_iterations, grid_iterations, optimizer

        ! Input read from the [reference] section.
        logical :: compare, compare_density, compare_isobar, compare_phase_transition
        real(dp) :: ref_phase_transition, ref_density, ref_density_temperature
        real(dp) :: ref_phase_transition_weight, ref_density_weight, ref_isobar_weight
        real(dp), dimension(:), allocatable :: ref_isobar_temperature, ref_isobar_volume
        character(len=:), allocatable :: ref_isobar_file_helper
        type(varying_string) :: ref_isobar_file

        ! Input read from the [output] section.
        logical :: contrib, helmholtz_contrib, internal_contrib, entropy_contrib, &
            cv_contrib, imode
        logical :: progress_bar
    end type input_data

      
    type(toml_table), allocatable :: table
    type(input_data) :: pmk_input

    !=====================================================================================
    contains
        !=================================================================================
        ! Processes the input from the input configuration file.
        subroutine process_input(table, input_file)
            type(toml_table), allocatable, intent(inout) :: table
            character(len=*), intent(in) :: input_file
            integer :: open_unit, ios
            type(toml_error), allocatable :: error

            ! Set default values. Make sure that these defaults are consistent.

            ! Defaults for the [system] section.
            pmk_input%components = 1

            ! Defaults for the [ensemble] section.
            call set_range(pmk_input%temperature, 298.15_dp, 298.15_dp, 1) ! in K
            pmk_input%pressure = 1.013250_dp ! in bar
            allocate(pmk_input%monomer_amounts(pmk_input%components))
            pmk_input%monomer_amounts = 1.0_dp/real(pmk_input%components, dp) ! in mol

            ! Defaults for the [qce] section.
            call set_range(pmk_input%amf, 0.0_dp, 0.0_dp, 1) ! in Jm^3/mol^2
            call set_range(pmk_input%bxv, 1.0_dp, 1.0_dp, 1)
            call set_range(pmk_input%amf_temp, 0.0_dp, 0.0_dp, 1) ! in Jm^3/(K mol^2)
            call set_range(pmk_input%bxv_temp, 0.0_dp, 0.0_dp, 1) ! in 1/K
            pmk_input%max_deviation = 1.0e-9_dp
            pmk_input%volume_damping_factor = 0.01_dp
            pmk_input%rotor_cutoff = 0.00_dp
            pmk_input%qce_iterations = 100
            pmk_input%newton_iterations = 500
            pmk_input%grid_iterations = 1
            pmk_input%optimizer = 0
            pmk_input%imode = .false.

            ! Defaults for the [reference] section.
            ! Changed to .true. for testing purposes.
            pmk_input%compare = .true.
            pmk_input%compare_isobar = .false.
            pmk_input%compare_density = .true.
            pmk_input%compare_phase_transition = .true.
            pmk_input%ref_isobar_weight = 1.0_dp
            pmk_input%ref_density_weight = 1.0_dp
            pmk_input%ref_phase_transition_weight = 1.0_dp

            ! Defaults for the [output] section.
            pmk_input%contrib = .false.
            pmk_input%helmholtz_contrib = .false.
            pmk_input%internal_contrib = .false.
            pmk_input%entropy_contrib = .false.
            pmk_input%cv_contrib = .false.
            pmk_input%progress_bar = .true.

            ! Open the QCE input file.
            block
                open(newunit=open_unit, file=input_file, status="old", iostat=ios)
                if (ios /= 0) call pmk_error("could not open '" // trim(input_file) // "'")   
                call toml_parse(table, open_unit, error)
                close(ios)
                if (allocated(error)) then
                  print '(a)', "Error: "//error%message
                  stop 1
                end if
            end block

            ! Overwrite the default values with user-specified values.
            call read_data(table, pmk_input)
            call convert_helpers(pmk_input)
        end subroutine process_input

        !=================================================================================

        !------------------------------------------------------------------------
        !> Read input data from TOML file
        !------------------------------------------------------------------------
        subroutine read_data(table, input)

            type(toml_table), intent(inout) :: table
            type(input_data), intent(out) :: input
            type(toml_table), pointer :: child
            type(toml_array), pointer :: array
            logical :: reverse
            integer :: ival
        
            ! The default values are set at the end of each get_value call
            !------------------------------------------------------------------------
            !> Read [system] section
            !------------------------------------------------------------------------
            call get_value(table, "system", child)
        
            !> components
            call get_value(child, "components", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%components)
                else 
                    call pmk_argument_count_error("components", "system")
                end if
            else 
                input%components = 1
            end if
        
            !------------------------------------------------------------------------
            !> Read [qce] section
            !------------------------------------------------------------------------
            call get_value(table, "qce", child)
            if (.not.associated(child)) then
              write(*,*) "Error: No [qce] section found in input file"
              stop
            end if
        
            !> amf
            call get_value(child, "amf", array, requested=.false.)
            if (associated(array)) then
                if (len(array) == 1) then 
                    call get_value(array, 1, input%amf%first)
                else if (len(array) == 3) then
                    call get_value(array, 1, input%amf%first)
                    call get_value(array, 2, input%amf%last)
                    call get_value(array, 3, input%amf%num)
                else
                    call pmk_argument_error("amf", "qce")
                end if  
            else 
                call set_range(pmk_input%amf, 0.0_dp, 0.0_dp, 1) ! in Jm^3/mol^2
            end if
        
            !> amf_temp
            call get_value(child, "amf_temp", array, requested=.false.)
            if (associated(array)) then
                if (len(array) == 1) then 
                    call get_value(array, 1, input%amf_temp%first)
                else if (len(array) == 3) then
                    call get_value(array, 1, input%amf_temp%first)
                    call get_value(array, 2, input%amf_temp%last)
                    call get_value(array, 3, input%amf_temp%num)
                else
                    call pmk_argument_error("amf_temp", "qce")
                end if  
            else 
                call set_range(pmk_input%amf_temp, 0.0_dp, 0.0_dp, 1) ! in Jm^3/(K mol^2)
            end if

            !> bxv
            call get_value(child, "bxv", array, requested=.false.)
            if (associated(array)) then
                if (len(array) == 1) then 
                    call get_value(array, 1, input%bxv%first)
                else if (len(array) == 3) then
                    call get_value(array, 1, input%bxv%first)
                    call get_value(array, 2, input%bxv%last)
                    call get_value(array, 3, input%bxv%num)
                else
                    call pmk_argument_error("bxv", "qce")
                end if  
            else 
                call set_range(pmk_input%bxv, 1.0_dp, 1.0_dp, 1)
            end if

          
            !> bxv_temp
            call get_value(child, "bxv", array, requested=.false.)
            if (associated(array)) then
                if (len(array) == 1) then 
                    call get_value(array, 1, input%bxv_temp%first)
                else if (len(array) == 3) then
                    call get_value(array, 1, input%bxv_temp%first)
                    call get_value(array, 2, input%bxv_temp%last)
                    call get_value(array, 3, input%bxv_temp%num)
                else
                    call pmk_argument_error("bxv_temp", "qce")
                end if  
            else 
                call set_range(pmk_input%bxv_temp, 0.0_dp, 0.0_dp, 1) ! in 1/K
            end if
          
            !> qce iterations
            call get_value(child, "iterations", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%qce_iterations)
                else 
                    call pmk_argument_count_error("qce_iterations", "qce")
                end if
            else 
                input%qce_iterations = 100
            end if
          
            !> newton iterations
            call get_value(child, "newton_iterations", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%newton_iterations)
                else 
                    call pmk_argument_count_error("newton_iterations", "qce")
                end if
            else 
                input%newton_iterations = 500
            end if
          
            !> grid iterations
            call get_value(child, "grid_iterations", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%grid_iterations)
                else 
                    call pmk_argument_count_error("grid_iterations", "qce")
                end if
            else 
                input%grid_iterations = 1
            end if
            
          
            !> optimizer
            call get_value(child, "optimizer", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%optimizer)
                else 
                    call pmk_argument_count_error("optimizer", "qce")
                end if
            else 
                input%optimizer = 0
            end if
          
            !> interface mode
            call get_value(child, "imode", array, requested=.false.)
            if (associated(array)) then
                input%imode = .true.
            else 
                input%imode = .false.
            end if
        
            !> maximum relative deviation
            call get_value(child, "max_deviation", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%max_deviation)
                else 
                    call pmk_argument_count_error("max_deviation", "qce")
                end if
            else 
                input%max_deviation = 1.0e-9_dp
            end if
        
            !> volume damping factor
            call get_value(child, "volume_damping_factor", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%volume_damping_factor)
                else 
                    call pmk_argument_count_error("volume_damping_factor", "qce")
                end if
            else 
                input%volume_damping_factor = 0.01_dp
            end if
        
            !> Read value from entry "rotor_cutoff"
            call get_value(child, "rotor_cutoff", array, requested=.false.)
            if (associated(array)) then
                if (len(array)==1) then
                    call get_value(array, 1, input%rotor_cutoff)
                else 
                    call pmk_argument_count_error("rotor_cutoff", "qce")
                end if
            else 
                input%rotor_cutoff = 0.00_dp
            end if
        
            !------------------------------------------------------------------------
            !> Read [ensemble] section
            !------------------------------------------------------------------------
            call get_value(table, "ensemble", child)
        
            !> temperature
            call get_value(child, "temperature", array, requested=.false.)
            if (associated(array)) then
                if (len(array) == 1) then 
                    call get_value(array, 1, input%temperature%first)
                else if (len(array) == 3) then
                    call get_value(array, 1, input%temperature%first)
                    call get_value(array, 2, input%temperature%last)
                    call get_value(array, 3, input%temperature%num)
                else
                    call pmk_argument_error("temperature", "ensemble")
                end if  
            else 
                call set_range(pmk_input%temperature, 298.15_dp, 298.15_dp, 1) ! in K
            end if

            !> pressure
            call get_value(child, "pressure", input%pressure, 1.013250_dp)
        
            ! monomer_amounts
            call get_value(child, "monomer_amounts", array)
            if (associated(array)) then
                allocate(input%monomer_amounts(len(array)))
                do ival = 1, size(input%monomer_amounts)
                  call get_value(array, ival, input%monomer_amounts(ival))
                end do
                call get_value(child, "reverse", reverse, .false.)
                if (reverse) input%monomer_amounts(:) = input%monomer_amounts(size(input%monomer_amounts):1:-1)
            else
                pmk_input%monomer_amounts = 1.0_dp/real(pmk_input%components, dp) ! in mol
            end if
          
            !------------------------------------------------------------------------
            !> Read [reference] section
            !------------------------------------------------------------------------
            ! The reference section is optional
            call get_value(table, "reference", child, requested=.false.)
            if (associated(child)) then
          
                !> reference temperature of phase transition
                !  if only one value is given, it's stored in ref_phase_transition
                !  if two values are given, the first is stored in ref_phase_transition, the second in 
                !  ref_phase_transition_weight
                call get_value(child, "phase_transition", array, requested=.false.)
                if (associated(array)) then
                    input%compare = .true.
                    input%compare_phase_transition = .true.
                    if (len(array) == 1) then
                        call get_value(array, 1, input%ref_phase_transition)
                        input%ref_phase_transition_weight = 1.0_dp
                    else if (len(array) == 2) then
                        call get_value(array, 1, input%ref_phase_transition)
                        call get_value(array, 2, input%ref_phase_transition_weight)
                    else
                        call pmk_argument_count_error("phase_transition", "reference")
                    end if
                else
                    input%compare = .false.
                    input%compare_phase_transition = .false.
                    input%ref_phase_transition_weight = 1.0_dp
                end if
            
                !> reference density
                !  if two values are given, the first is stored in ref_density_temperature, the second in ref_density
                !  if three values are given, the first is stored in ref_density_temperature, the second in ref_density 
                !  and the third in ref_density_weight
                call get_value(child, "density", array, requested=.false.)
                if (associated(array)) then
                    print *, "density associated"
                    input%compare = .true.
                    input%compare_density = .true.
                    if (len(array) == 2) then 
                        call get_value(array, 1, input%ref_density_temperature)
                        call get_value(array, 2, input%ref_density)
                        input%ref_density_weight = 1.0_dp
                    else if (len(array) == 3) then
                        call get_value(array, 1, input%ref_density_temperature)
                        call get_value(array, 2, input%ref_density)
                        call get_value(array, 3, input%ref_density_weight)
                    else
                        call pmk_argument_count_error("density", "reference")
                    end if  
                else 
                    input%compare_density = .false.
                    input%ref_density_weight = 1.0_dp
                end if
            
                !> reference isobar 
                ! if reference isobar is given, compare_isobar is set to true
                call get_value(child, "isobar", array, requested=.false.)
                if (associated(array)) then
                    input%compare = .true.
                    input%compare_isobar = .true.
                    if (len(array) == 1) then 
                        call get_value(array, 1, input%ref_isobar_file_helper)
                    else if (len(array) == 2) then
                        call get_value(array, 1, input%ref_isobar_file_helper)
                        call get_value(array, 2, input%ref_isobar_weight)
                    else
                        call pmk_argument_count_error("isobar", "reference")
                    end if
                else
                    input%compare_isobar = .false.
                    input%ref_isobar_weight = 1.0_dp
                end if

            else 
            ! Defaults for the [reference] section.
                input%compare = .false.
                input%compare_isobar = .false.
                input%compare_density = .false.
                input%compare_phase_transition = .false.
                input%ref_isobar_weight = 1.0_dp
                input%ref_density_weight = 1.0_dp
                input%ref_phase_transition_weight = 1.0_dp
            end if             
              
            !------------------------------------------------------------------------
            !> Read [output] section
            !------------------------------------------------------------------------
            call get_value(table, "output", child)
              
            !> cortribuion
            call get_value(child, "contributions", input%contrib, .false.)
              
            !> helmholtz contribution
            call get_value(child, "helmholtz_contributions", input%helmholtz_contrib, .false.)
              
            !> internal contributions
            call get_value(child, "internal_contributions", input%internal_contrib, .false.)
              
            !> entropy contributions
            call get_value(child, "entropy_contributions", input%entropy_contrib, .false.)
              
            !> cv contributions
            call get_value(child, "cv_contributions", input%cv_contrib, .false.)
              
            ! progress bar
            call get_value(child, "progress_bar", input%progress_bar, .false.)
              
        end subroutine read_data
        
        !> Subroutine to convert helpers in range_t
        subroutine convert_helpers(input)
            type(input_data), intent(inout) :: input
            integer :: i

            if (allocated(input%ref_isobar_file_helper)) then
                do i = 1, len(input%ref_isobar_file_helper)
                    if (input%ref_isobar_file_helper(i:i) == ' ') then
                      input%ref_isobar_file_helper(i:i) = '_'
                    end if
                end do
            end if
            input%ref_isobar_file = trim(input%ref_isobar_file_helper)

        end subroutine convert_helpers


        !=================================================================================
        ! Performs sanity checks on the input.
        subroutine check_input()
            write(*,*) "Checking input..."
            ! Check amf.
            if (pmk_input%amf%first < 0.0_dp) &
                call pmk_unphysical_argument_error("amf", "qce")
            if (.not. check_range(pmk_input%amf)) &
                call pmk_illegal_range_error("amf", "qce")

            ! Check bxv.
            if (pmk_input%bxv%first <= 0.0_dp) &
                call pmk_unphysical_argument_error("bxv", "qce")
            if (.not. check_range(pmk_input%bxv)) &
                call pmk_illegal_range_error("bxv", "qce")

            ! Check qce_iterations.
            if (pmk_input%qce_iterations <= 0_dp) &
                call pmk_argument_error("qce_iterations", "qce")

            ! Check newton_iterations.
            if (pmk_input%newton_iterations <= 0_dp) &
                call pmk_argument_error("newton_iterations", "qce")

            ! Check grid_iterations.
            if (pmk_input%grid_iterations <= 0_dp) &
                call pmk_argument_error("grid_iterations", "qce")

            ! Check temperature.
            if (pmk_input%temperature%first < 0.0_dp) &
                call pmk_unphysical_argument_error("temperature", "ensemble")
            if (.not. check_range(pmk_input%temperature)) &
                call pmk_illegal_range_error("temperature", "ensemble")

            ! Check pressure.
            if (pmk_input%pressure < 0.0_dp) &
                call pmk_unphysical_argument_error("pressure", "ensemble")

            ! Check components.
            if (pmk_input%components <= 0) &
                call pmk_unphysical_argument_error("components", "system")

            ! Check monomer amounts.
            if (any(pmk_input%monomer_amounts < 0.0_dp)) &
                call pmk_unphysical_argument_error("monomer_amounts", "ensemble")
            if (any(pmk_input%monomer_amounts < global_eps)) &
                call pmk_error("monomer_amounts must not be zero (indeterminate system); treat as pure substance instead")

            ! Check maximum relative relative deviation.
            if (pmk_input%max_deviation <= 0.0_dp) &
                call pmk_argument_error("max_deviation", "qce")

            ! Check volume damping factor.
            if (pmk_input%volume_damping_factor <= 0.0_dp .or. &
                pmk_input%volume_damping_factor >= 1.0_dp) &
                call pmk_argument_error("volume_damping_factor", "qce")

            ! Check free rotator frequency threshold.
            if (pmk_input%rotor_cutoff < 0.0_dp) &
                call pmk_argument_error("rotor_cutoff", "qce")

            ! Check reference section.
            if (pmk_input%amf%num*pmk_input%bxv%num > 1 .and. .not. pmk_input%compare) &
                call pmk_error("amf and/or bxv interval specified, " // &
                "but no reference section")

            ! Check density.
            if (pmk_input%compare_density) then
                if (pmk_input%ref_density <= 0.0_dp) &
                    call pmk_unphysical_argument_error("density", "reference")
                if (pmk_input%ref_density_weight < 0.0_dp) &
                    call pmk_unphysical_argument_error("density", "reference")
                if (pmk_input%ref_density_temperature <= 0.0_dp) &
                    call pmk_unphysical_argument_error("density", "reference")
                if (pmk_input%ref_density_temperature <= pmk_input%temperature%first .or. &
                    pmk_input%ref_density_temperature >= pmk_input%temperature%last) &
                    call pmk_error("density reference temperature must " // &
                    "be within the investigated temperature range")
            end if

            ! Check phase transition.
            if (pmk_input%compare_phase_transition) then
                if (pmk_input%ref_phase_transition < 0.0_dp) &
                    call pmk_unphysical_argument_error("phase_transition", "reference")
                if (pmk_input%ref_phase_transition_weight < 0.0_dp) &
                    call pmk_unphysical_argument_error("phase_transition", "reference")
                if (pmk_input%ref_phase_transition <= pmk_input%temperature%first .or. &
                    pmk_input%ref_phase_transition >= pmk_input%temperature%last) &
                    call pmk_error("reference temperature of phase transition must " // &
                    "be within the investigated temperature range")
                if (pmk_input%temperature%num < 2) &
                    call pmk_error("need two temperature points for phase " // &
                    "transition determination")
            end if

            ! Check isobar.
            if (pmk_input%compare_isobar) then
                if (pmk_input%ref_isobar_weight < 0.0_dp) &
                    call pmk_unphysical_argument_error("isobar", "reference")
                if (any(pmk_input%ref_isobar_temperature < 0.0_dp)) &
                    call pmk_unphysical_argument_error("isobar", "reference")
                if (any(pmk_input%ref_isobar_volume <= 0.0_dp)) &
                    call pmk_unphysical_argument_error("isobar", "reference")
                if (any(pmk_input%ref_isobar_temperature < &
                    pmk_input%temperature%first) .or. &
                    any(pmk_input%ref_isobar_temperature > pmk_input%temperature%last)) &
                    call pmk_error("reference isobar temperatures must be within " // &
                    "the investigated temperature range")

            end if
        end subroutine check_input
        !=================================================================================
        ! Prints the processed input.
        subroutine print_input()
            character(fmt_len) :: monomer_amounts_fmt

            ! Write format strings.
            write(monomer_amounts_fmt, '(A, G0, A)') &
                "(12X,A,", pmk_input%components, "(1X,G0.6),1X,A)"

            ! Print input.
            write(*, '(4X,A)') 'Using the following input:'
            write(*, *) 

            ! Print [system] section.
            write(*, '(8X,A)') "[system]"
            write(*, '(12X,A,1X,G0)') "components:", pmk_input%components
            write(*, *)

            ! Print [ensemble] section.
            write(*, '(8X,A)') "[ensemble]"
            write(*, '(12X,A,1X,G0.6,1X,A)') "pressure:", pmk_input%pressure, "[bar]"
            write(*, '(12X,A,1X)', advance = "no") "temperature:"
                call write_range(pmk_input%temperature)
                write(*, '(1X,A)') "[K]"
            write(*, monomer_amounts_fmt) &
                "monomer amounts:", pmk_input%monomer_amounts, "[mol]"
            write(*, *)

            ! Print [qce] section.
            write(*, '(8X,A)') "[qce]"

            write(*, '(12X,A,1X)', advance = "no") "amf:"
            call write_range(pmk_input%amf)
            write(*, '(1X,A)') "[J*m^3/mol^2]"
            write(*, '(12X,A,1X)', advance = "no") "bxv:"
            call write_range(pmk_input%bxv)
            write(*, *)
            write(*, '(12X,A,1X,G0)') "maximum number of QCE iterations:", &
                pmk_input%qce_iterations
            write(*, '(12X,A,1X,G0)') "maximum number of Newton-Raphson iterations:", &
                pmk_input%newton_iterations
            write(*, '(12X,A,1X,G0.6)') "maximum relative deviation:", &
                pmk_input%max_deviation
            write(*, '(12X,A,1X,G0)') "number of grid iterations:", &
                pmk_input%grid_iterations
            if (pmk_input%optimizer > 0) then
               write(*, '(12X,A)') "optimizer: Downhill-Simplex"
            end if
            write(*, '(12X,A,1X,G0.6)') "volume damping factor:", &
                pmk_input%volume_damping_factor
            write(*, '(12X,A,1X,G0.6)') "free rotator correction threshold frequency:", &
                pmk_input%rotor_cutoff
            write(*, *)

            ! Print [reference] section.
            if (pmk_input%compare) then
                write(*,'(8X,A)') "[reference]"
                if (pmk_input%compare_isobar) then
                    write(*,'(12X,A,1X,G0.6,A,A)') "isobar (weight =", &
                        pmk_input%ref_isobar_weight, "): " // &
                        char(pmk_input%ref_isobar_file)
                end if
                if (pmk_input%compare_density) then
                    write(*,'(12X,A,1X,G0.6,A,1X,G0.6,1X,G0.6,1X,A)') &
                        "density (weight =", pmk_input%ref_density_weight, "):", &
                        pmk_input%ref_density_temperature, &
                        pmk_input%ref_density , "[K; g/cm^3]"
                end if
                if (pmk_input%compare_phase_transition) then
                    write(*,'(12X,A,1X,G0.6,A,1X,G0.6,1X,A)') &
                        "phase transition (weight =", &
                        pmk_input%ref_phase_transition_weight, "):", &
                        pmk_input%ref_phase_transition, "[K]"
                end if
                write(*,*)
            end if
        end subroutine print_input
        !=================================================================================
        ! Reads the [system] section of the input.
        subroutine read_system_section(cfg)
            type(config_t), intent(inout) :: cfg
            type(record_t), pointer :: p

            ! Set components.
            p => cfg%get_record("system", "components")
            if (associated(p)) call process_components_record(p)
        end subroutine read_system_section
        !=================================================================================
        ! Reads the [ensemble] section of the input.
        subroutine read_ensemble_section(cfg)
            type(config_t), intent(inout) :: cfg
            type(record_t), pointer :: p

            ! Set temperature.
            p => cfg%get_record("ensemble", "temperature")
            if (associated(p)) call process_temperature_record(p)

            ! Set pressure.
            p => cfg%get_record("ensemble", "pressure")
            if (associated(p)) call process_pressure_record(p)

            ! Set monomers amounts.
            p => cfg%get_record("ensemble", "monomer_amounts")
            if (associated(p)) call process_monomer_amounts_record(p)
        end subroutine read_ensemble_section
        !=================================================================================
        ! Reads the [qce] section of the input.
        subroutine read_qce_section(cfg)
            type(config_t), intent(inout) :: cfg
            type(record_t), pointer :: p

            ! Set amf.
            p => cfg%get_record("qce", "amf")
            if (associated(p)) call process_amf_record(p)

            ! Set amf_temp.
            p => cfg%get_record("qce", "amf_temp")
            if (associated(p)) call process_amf_temp_record(p)

            ! Set bxv.
            p => cfg%get_record("qce", "bxv")
            if (associated(p)) call process_bxv_record(p)

            ! Set bxv_temp.
            call set_range(pmk_input%bxv_temp, 0.0_dp, 0.0_dp, 1)
            p => cfg%get_record("qce", "bxv_temp")
            if (associated(p)) call process_bxv_temp_record(p)

            ! Set qce_iterations.
            p => cfg%get_record("qce", "qce_iterations")
            if (associated(p)) call process_qce_iterations_record(p)

            ! Set newton_iterations.
            p => cfg%get_record("qce", "newton_iterations")
            if (associated(p)) call process_newton_iterations_record(p)

            ! Set grid_iterations.
            p => cfg%get_record("qce", "grid_iterations")
            if (associated(p)) call process_grid_iterations_record(p)

            ! Set optimizer.
            p => cfg%get_record("qce", "optimizer")
            if (associated(p)) call process_optimizer_record(p)

            ! Get interface mode
            p => cfg%get_record("qce", "interface_mode")
            if (associated(p)) then
                pmk_input%imode = .true.
            end if

            ! Set maximum relative deviation
            p => cfg%get_record("qce", "max_deviation")
            if (associated(p)) call process_max_deviation_record(p)

            ! Set volume damping factor
            p => cfg%get_record("qce", "volume_damping_factor")
            if (associated(p)) call process_volume_damping_factor_record(p)

            ! Set volume damping factor
            p => cfg%get_record("qce", "rotor_cutoff")
            if (associated(p)) call process_rotor_cutoff_record(p)
        end subroutine read_qce_section
        !=================================================================================
        ! Reads the [reference] section of the input.
        subroutine read_reference_section(cfg)
            type(config_t), intent(inout) :: cfg
            type(record_t), pointer :: p

            ! Set reference temperature of phase transition and weight.
            p=> cfg%get_record("reference", "phase_transition")
            if (associated(p)) then
                pmk_input%compare = .true.
                pmk_input%compare_phase_transition = .true.
                call process_phase_transition_record(p)
            end if

            ! Set reference density and weight.
            p=> cfg%get_record("reference", "density")
            if (associated(p)) then
                pmk_input%compare = .true.
                pmk_input%compare_density = .true.
                call process_density_record(p)
            end if

            ! Set reference isobar and weight.
            p=> cfg%get_record("reference", "isobar")
            if (associated(p)) then
                pmk_input%compare = .true.
                pmk_input%compare_isobar = .true.
                call process_isobar_record(p)
            end if
        end subroutine read_reference_section
        !=================================================================================
        ! Read the [output] section of the input.
        subroutine read_output_section(cfg)
            type(config_t), intent(inout) :: cfg
            type(record_t), pointer :: p, p2

            ! Get contributions.
            p => cfg%get_record("output", "contributions")
            if (associated(p)) then
                call process_contrib_record(p)
            end if

            ! Get progress.
            p => cfg%get_record("output", "progress")
            p2 => cfg%get_record("output", "noprogress")
            if (associated(p) .and. associated(p2)) then
                call pmk_error("cannot have both 'progress' and no 'progress' in section [output]")
            end if
            if (associated(p)) then
                pmk_input%progress_bar = .true.
            else if (associated(p2)) then
                pmk_input%progress_bar = .false.
            end if
        end subroutine read_output_section
        !=================================================================================
        ! Processes the contributions record.
        subroutine process_contrib_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: iarg

            if (p%nr_args == 0) then
                pmk_input%contrib = .true.
            else
                do iarg = 1, p%nr_args
                    select case (char(p%args(iarg)))
                        case ("helmholtz")
                            pmk_input%helmholtz_contrib =  .true.
                        case ("internal")
                            pmk_input%internal_contrib =  .true.
                        case ("entropy")
                            pmk_input%entropy_contrib =  .true.
                        case ("cv")
                            pmk_input%cv_contrib =  .true.
                    end select
                end do
            end if
        end subroutine process_contrib_record
        !=================================================================================
        ! Process the temperature record.
        subroutine process_temperature_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios
            real(dp) :: first, last
            integer :: num

            if (p%nr_args == 1) then
                ! Single temperature point.
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("temperature", "ensemble")
                call set_range(pmk_input%temperature, first, first, 1)
            else if (p%nr_args == 3) then
                ! Temperature range specified by initial temperature (t0), final
                ! temeprature (t1), and number of temperature points (n).
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("temperature", "ensemble")
                last = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("temperature", "ensemble")
                num = string2int(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("temperature", "ensemble")
                call set_range(pmk_input%temperature, first, last, num)
            else
                call pmk_argument_count_error("temperature", "ensemble")
            end if
        end subroutine process_temperature_record
        !=================================================================================
        ! Process the pressure record.
        subroutine process_pressure_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%pressure = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("pressure", "ensemble")
            else
                call pmk_argument_count_error("pressure", "ensemble")
            end if
        end subroutine process_pressure_record
        !=================================================================================
        ! Process the monomer amounts record.
        subroutine process_monomer_amounts_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: i, ios

            ! Deallocate array, as it is allocated by default.
            deallocate(pmk_input%monomer_amounts)

            if (p%nr_args == pmk_input%components) then
                allocate(pmk_input%monomer_amounts(pmk_input%components))
                do i = 1, pmk_input%components
                    pmk_input%monomer_amounts(i) = string2real(p%args(i), ios)
                    if (ios /= 0) &
                        call pmk_argument_error("amounts", "ensemble")
                end do
            else
                call pmk_argument_count_error("amounts", "ensemble")
            end if
        end subroutine process_monomer_amounts_record
        !=================================================================================
        ! Process the amf record.
        subroutine process_amf_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios
            real(dp) :: first, last
            integer :: num

            if (p%nr_args == 1) then
                ! single amf value
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("amf", "qce")
                call set_range(pmk_input%amf, first, first, 1)
            else if (p%nr_args == 3) then
                ! amf range
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("amf", "qce")
                last = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("amf", "qce")
                num = string2int(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("amf", "qce")
                call set_range(pmk_input%amf, first, last, num)
            else
                call pmk_argument_count_error("amf", "qce")
            end if
        end subroutine process_amf_record
        !=================================================================================
        ! Process the amf_temp record.
        subroutine process_amf_temp_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios
            real(dp) :: first, last
            integer :: num

            if (p%nr_args == 1) then
                ! single amf_temp value
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("amf_temp", "qce")
                call set_range(pmk_input%amf_temp, first, first, 1)
            else if (p%nr_args == 3) then
                ! amf_temp range
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("amf_temp", "qce")
                last = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("amf_temp", "qce")
                num = string2int(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("amf_temp", "qce")
                call set_range(pmk_input%amf_temp, first, last, num)
            else
                call pmk_argument_count_error("amf_temp", "qce")
            end if
        end subroutine process_amf_temp_record
        !=================================================================================
        ! Process the bxv record.
        subroutine process_bxv_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios
            real(dp) :: first, last
            integer :: num

            if (p%nr_args == 1) then
                ! single bxv value
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("bxv", "qce")
                call set_range(pmk_input%bxv, first, first, 1)
            else if (p%nr_args == 3) then
                ! bxv range
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("bxv", "qce")
                last = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("bxv", "qce")
                num = string2int(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("bxv", "qce")
                call set_range(pmk_input%bxv, first, last, num)
            else
                call pmk_argument_count_error("bxv", "qce")
            end if
        end subroutine process_bxv_record
        !=================================================================================
        ! Process the bxv_temp record.
        subroutine process_bxv_temp_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios
            real(dp) :: first, last
            integer :: num

            if (p%nr_args == 1) then
                ! single bxv_temp value
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("bxv_temp", "qce")
                call set_range(pmk_input%bxv_temp, first, first, 1)
            else if (p%nr_args == 3) then
                ! bxv_temp range
                first = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("bxv_temp", "qce")
                last = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("bxv_temp", "qce")
                num = string2int(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("bxv_temp", "qce")
                call set_range(pmk_input%bxv_temp, first, last, num)
            else
                call pmk_argument_count_error("bxv_temp", "qce")
            end if
        end subroutine process_bxv_temp_record
        !=================================================================================
        ! Process the qce_iterations record.
        subroutine process_qce_iterations_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%qce_iterations = string2int(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("qce_iterations", "qce")
            else
                call pmk_argument_count_error("qce_iterations", "qce")
            end if
        end subroutine process_qce_iterations_record
        !=================================================================================
        ! Process the newton_iterations record.
        subroutine process_newton_iterations_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%newton_iterations = string2int(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("newton_iterations", "qce")
            else
                call pmk_argument_count_error("newton_iterations", "qce")
            end if
        end subroutine process_newton_iterations_record
        !=================================================================================
        ! Process the grid_iterations record.
        subroutine process_grid_iterations_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%grid_iterations = string2int(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("grid_iterations", "qce")
            else
                call pmk_argument_count_error("grid_iterations", "qce")
            end if
        end subroutine process_grid_iterations_record
        !=================================================================================
        ! Process the optimizer record.
        subroutine process_optimizer_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%optimizer = string2int(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("optimizer", "qce")
            else
                call pmk_argument_count_error("optimizer", "qce")
            end if
        end subroutine process_optimizer_record
        !=================================================================================
        subroutine process_max_deviation_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%max_deviation = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("max_deviation", "qce")
            else
                call pmk_argument_count_error("max_deviation", "qce")
            end if
        end subroutine process_max_deviation_record
        !=================================================================================
        subroutine process_volume_damping_factor_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%volume_damping_factor = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("volume_damping_factor", "qce")
            else
                call pmk_argument_count_error("volume_damping_factor", "qce")
            end if
        end subroutine process_volume_damping_factor_record
        !=================================================================================
        subroutine process_rotor_cutoff_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%rotor_cutoff = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("rotor_cutoff", "qce")
            else
                call pmk_argument_count_error("rotor_cutoff", "qce")
            end if
        end subroutine process_rotor_cutoff_record
        !=================================================================================
        subroutine process_components_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%components = string2int(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("components", "system")
            else
                call pmk_argument_count_error("components", "system")
            end if
        end subroutine process_components_record
        !=================================================================================
        subroutine process_density_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 2) then
                pmk_input%ref_density_temperature = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("density", "reference")
                pmk_input%ref_density = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("density", "reference")
            else if (p%nr_args == 3) then
                pmk_input%ref_density_temperature = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("density", "reference")
                pmk_input%ref_density = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("density", "reference")
                pmk_input%ref_density_weight = string2real(p%args(3), ios)
                if (ios /= 0) call pmk_argument_error("density", "reference")
            else
                call pmk_argument_count_error("density", "reference")
            end if
        end subroutine process_density_record
        !=================================================================================
        subroutine process_phase_transition_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%ref_phase_transition = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("phase_transition", "reference")
            else if (p%nr_args == 2) then
                pmk_input%ref_phase_transition = string2real(p%args(1), ios)
                if (ios /= 0) call pmk_argument_error("phase_transition", "reference")
                pmk_input%ref_phase_transition_weight = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("phase_transition", "reference")
            else
                call pmk_argument_count_error("phase_transition", "reference")
            end if
        end subroutine process_phase_transition_record
        !=================================================================================
        subroutine process_isobar_record(p)
            type(record_t), pointer, intent(in) :: p
            integer :: ios

            if (p%nr_args == 1) then
                pmk_input%ref_isobar_file = p%args(1)
            else if (p%nr_args == 2) then
                pmk_input%ref_isobar_file = p%args(1)
                pmk_input%ref_isobar_weight = string2real(p%args(2), ios)
                if (ios /= 0) call pmk_argument_error("isobar", "reference")
            else
                call pmk_argument_count_error("isobar", "reference")
            end if

            ! Read isobar file.
            call read_isobar_file()
        end subroutine process_isobar_record
        !=================================================================================
        subroutine read_isobar_file()
            integer :: my_unit, ios, n, i
            character(:), allocatable :: fn

            ! Open unit.
            allocate(fn, source = char(pmk_input%ref_isobar_file))
            open(newunit = my_unit, file = fn, action = 'read', status = 'old', &
                iostat = ios)
            if (ios /= 0) call pmk_error("could not open '" // fn // "'")

            ! Count lines.
            n = 0
            count_loop: do
                read(my_unit, *, iostat = ios)
                if (ios /= 0) exit count_loop
                n = n + 1
            end do count_loop
            rewind(my_unit)
            if (n == 0) call pmk_error("'" // fn // "' is empty")

            ! Read frequencies.
            allocate(pmk_input%ref_isobar_temperature(n))
            allocate(pmk_input%ref_isobar_volume(n))
            read_loop: do i = 1, n
                read(my_unit, *, iostat = ios) pmk_input%ref_isobar_temperature(i), &
                    pmk_input%ref_isobar_volume(i)
                if (ios > 0) then
                    call pmk_error("could not read '" // fn // "'")
                else if (ios < 0) then
                    exit read_loop
                end if
            end do read_loop

            close(my_unit)
        end subroutine read_isobar_file
        !=================================================================================
end module input
!=========================================================================================

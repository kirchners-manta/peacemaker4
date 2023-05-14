!module reader
!    use tomlf, only : toml_table, toml_array, get_value, len
!    implicit none
!    private
!  
!    public :: read_data
!  
!    contains
!  
!      !------------------------------------------------------------------------
!      !> Read input data from TOML file
!      !------------------------------------------------------------------------
!      subroutine read_data(table, input)
!
!          type(toml_table), intent(inout) :: table
!          type(input_data), intent(out) :: input
!          type(toml_table), pointer :: child
!          type(toml_array), pointer :: array
!          logical :: reverse
!          integer :: ival
!
!          ! The default values are set at the end of each get_value call
!          !------------------------------------------------------------------------
!          !> Read [system] section
!          !------------------------------------------------------------------------
!          call get_value(table, "system", child)
!
!          !> components
!          call get_value(child, "components", input%components)
!
!          !------------------------------------------------------------------------
!          !> Read [qce] section
!          !------------------------------------------------------------------------
!          call get_value(table, "qce", child)
!          if (.not.associated(child)) then
!            write(*,*) "Error: No [qce] section found in input file"
!            stop
!          end if
!
!          !> amf
!          call get_value(child, "amf", array)
!          allocate(input%amf_helper(len(array)))
!          do ival = 1, size(input%amf_helper)
!            call get_value(array, ival, input%amf_helper(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%amf_helper(:) = input%amf_helper(size(input%amf_helper):1:-1)
!
!          !> amf_temp
!          call get_value(child, "amf_temp", array)
!          allocate(input%amf_temp_helper(len(array)))
!          do ival = 1, size(input%amf_temp_helper)
!            call get_value(array, ival, input%amf_temp_helper(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%amf_temp_helper(:) = input%amf_temp_helper(size(input%amf_temp_helper):1:-1)
!
!          !> bxv
!          call get_value(child, "bxv", array)
!          allocate(input%bxv_helper(len(array)))
!          do ival = 1, size(input%bxv_helper)
!            call get_value(array, ival, input%bxv_helper(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%bxv_helper(:) = input%bxv_helper(size(input%bxv_helper):1:-1)
!
!          !> bxv_temp
!          call get_value(child, "bxv_temp", array)
!          allocate(input%bxv_temp_helper(len(array)))
!          do ival = 1, size(input%bxv_temp_helper)
!            call get_value(array, ival, input%bxv_temp_helper(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%bxv_temp_helper(:) = input%bxv_temp_helper(size(input%bxv_temp_helper):1:-1)
!
!          !> qce iterations
!          call get_value(child, "iterations", input%qce_iterations)
!
!          !> newton iterations
!          call get_value(child, "newton_iterations", input%newton_iterations)
!
!          !> grid iterations
!          call get_value(child, "grid_iterations", input%grid_iterations)
!
!          !> optimizer
!          call get_value(child, "optimizer", input%optimizer)
!
!          !> interface mode (imode, logical)
!          !  if imode = true, then the interface is used
!          call get_value(child, "imode", input%imode)
!
!          !> maximum relative deviation
!          call get_value(child, "max_deviation", input%max_deviation)
!
!          !> volume damping factor
!          call get_value(child, "volume_damping_factor", input%volume_damping_factor)
!
!          !> Read value from entry "rotor_cutoff"
!          call get_value(child, "rotor_cutoff", input%rotor_cutoff)
!
!          !------------------------------------------------------------------------
!          !> Read [ensemble] section
!          !------------------------------------------------------------------------
!          call get_value(table, "ensemble", child)
!
!          !> temperature
!          call get_value(child, "temperature", array)
!          allocate(input%temperature(len(array)))
!          do ival = 1, size(input%temperature)
!            call get_value(array, ival, input%temperature(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%temperature(:) = input%temperature(size(input%temperature):1:-1)
!
!          !> pressure
!          call get_value(child, "pressure", input%pressure)
!
!          ! monomer_amounts
!          call get_value(child, "monomer_amounts", array)
!          allocate(input%monomer_amounts(len(array)))
!          do ival = 1, size(input%monomer_amounts)
!            call get_value(array, ival, input%monomer_amounts(ival))
!          end do
!          call get_value(child, "reverse", reverse, .false.)
!          if (reverse) input%monomer_amounts(:) = input%monomer_amounts(size(input%monomer_amounts):1:-1)
!
!          !------------------------------------------------------------------------
!          !> Read [reference] section
!          !------------------------------------------------------------------------
!          call get_value(table, "reference", child)
!
!          !> reference temperature of phase transition
!          !  if only one value is given, it's stored in ref_phase_transition
!          !  if two values are given, the first is stored in ref_phase_transition, the second in 
!          !  ref_phase_transition_weight
!          call get_value(child, "phase_transition", array)
!          if (associated(array)) then
!            do ival = 1, len(array)
!              if (ival == 1) then
!                call get_value(array, ival, input%ref_phase_transition)
!              else if (ival == 2) then
!                call get_value(array, ival, input%ref_phase_transition_weight)
!              end if
!            end do
!          else 
!            call get_value(child, "phase_transition", input%ref_phase_transition)      
!          end if
!
!          !> reference density
!          !  if two values are given, the first is stored in ref_density_temperature, the second in ref_density
!          !  if three values are given, the first is stored in ref_density_temperature, the second in ref_density 
!          !  and the third in ref_density_weight
!          call get_value(child, "density", array)
!          if (associated(array)) then
!            do ival = 1, len(array)
!              if (ival == 1) then
!                call get_value(array, ival, input%ref_density_temperature)
!              else if (ival == 2) then
!                call get_value(array, ival, input%ref_density)
!              else if (ival == 3) then
!                call get_value(array, ival, input%ref_density_weight)
!              end if
!            end do
!          else 
!            ! error catching     
!          end if
!
!          !> reference isobar 
!          call get_value(child, "isobar_file", input%ref_isobar_file_helper)
!          call get_value(child, "isobar_weight", input%ref_isobar_weight)
!
!
!          !------------------------------------------------------------------------
!          !> Read [output] section
!          !------------------------------------------------------------------------
!          call get_value(table, "output", child)
!
!          !> cortribuion
!          call get_value(child, "contributions", input%contrib)
!
!          !> helmholtz contribution
!          call get_value(child, "helmholtz_contributions", input%helmholtz_contrib)
!
!          !> internal contributions
!          call get_value(child, "internal_contributions", input%internal_contrib)
!
!          !> entropy contributions
!          call get_value(child, "entropy_contributions", input%entropy_contrib)
!
!          !> cv contributions
!          call get_value(child, "cv_contributions", input%cv_contrib)
!
!          ! progress bar
!          call get_value(child, "progress_bar", input%progress_bar)
!          
!      end subroutine read_data
!
!        ! Subroutine to convert helpers in range_t
!      subroutine convert_helpers(input)
!          type(input_data), intent(inout) :: input
!          integer :: i
!
!          input%amf%first = input%amf_helper(1)
!          input%amf%last = input%amf_helper(2)
!          input%amf%delta = input%amf_helper(3)
!
!          input%bxv%first = input%bxv_helper(1)
!          input%bxv%last = input%bxv_helper(2)
!          input%bxv%delta = input%bxv_helper(3)
!
!          input%amf_temp%first = input%amf_temp_helper(1)
!          input%amf_temp%last = input%amf_temp_helper(2)
!          input%amf_temp%delta = input%amf_temp_helper(3)
!
!          input%bxv_temp%first = input%bxv_temp_helper(1)
!          input%bxv_temp%last = input%bxv_temp_helper(2)
!          input%bxv_temp%delta = input%bxv_temp_helper(3)
!
!          do i = 1, len(input%ref_isobar_file_helper)
!            if (input%ref_isobar_file_helper(i:i) == ' ') then
!              input%ref_isobar_file_helper(i:i) = '_'
!            end if
!          end do
!          input%ref_isobar_file = trim(input%ref_isobar_file_helper)
!      end subroutine convert_helpers
!
!    
!  end module reader
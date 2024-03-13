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
! This module provides the cluster_t data type, which represents a cluster, and associated
! procedures. The cluster_t data type contains static information about the clusters,
! that are shared across all QCE calculations. By a QCE calculation we mean a complete
! cycle of iterations for a given thermodynamic state (N, p, T) and parameter set 
! (amf, bxv). Peacemaker usually performs many of these calculations, as it samples
! temperatures, parameters, etc. Thus quantities such as moments of inertia and
! frequencies go in here, but temperature or pressure don't, as these may differ across
! multiple QCE calculations.
module cluster
    use kinds
    use iso_varying_string
    use error
    use input
    use auxiliary
    use constants

    !> TOML
    use tomlf, only : toml_table, toml_parse, toml_error, toml_array, get_value, len, toml_key

    implicit none
    private
    !=====================================================================================
    ! Public entities.
    public :: clusterset
    public :: monomer
    public :: process_clusterset
    public :: print_clusterset
    public :: check_clusterset
    public :: cluster_table
    public :: toml_parse
    public :: cluster_t
    public :: process_coordinates_record, diagonalize_3x3, center_of_mass
    !=====================================================================================
    ! The cluster_t data type, which represents a cluster. Reasonable defaults should be
    ! set here, if there are any.
    type :: cluster_t
        ! Total mass of the cluster in amu.
        real(dp) :: mass
        ! Rotational symmetry number sigma.
        integer :: sigma = 1
        ! Vibrational frequencies in 1/cm.
        real(dp), dimension(:), allocatable :: frequencies
        ! Principal moments of inertia in amu*Angstrom^2.
        real(dp), dimension(:), allocatable :: inertia
        ! Composition array (number of monomers).
        integer, dimension(:), allocatable :: composition
        ! Adiabatic interaction energy in kJ/mol.
        real(dp) :: energy
        ! Cluster volume in Angstrom^3
        real(dp) :: volume
        ! Frequency scaling factor.
        real(dp) :: fscale = 1.0_dp
        ! Anharmonicity constant
        real(dp) :: anharmonicity = 0.0_dp
        ! Flags classifying the cluster.
        logical :: linear = .false.
        logical :: atom = .false.
        logical :: monomer = .false.
        ! Cluster label.
        type(varying_string) :: label
    end type cluster_t

    type(toml_table), allocatable :: cluster_table

    !=====================================================================================
    ! The clusterset array and monomer array. The clusterset array is the central quantity
    ! in Peacemaker. It contains all the static information about all the clusters. The
    ! monomer array contains the indices of monomers in the cluster set array.
    type(cluster_t), dimension(:), allocatable, target :: clusterset
    integer, dimension(:), allocatable :: monomer

    !=====================================================================================
    contains
        !=================================================================================
        !> Subroutine reading in the clusterset data from the TOML configuration file.
        subroutine process_clusterset(cluster_table, clusterset_file)

            type(toml_table), allocatable, intent(inout) :: cluster_table
            character(len=*) :: clusterset_file
            integer :: open_unit, ios
            type(toml_error), allocatable :: error

            !type(cluster_t), dimension(:), allocatable, target :: clusterset
            !integer, dimension(:), allocatable :: monomer
            type(toml_table), pointer :: child
            type(toml_array), pointer :: array
            logical :: reverse
            integer :: ival, i, j
            integer :: nr_clusters

            !> A pointer to the current cluster.
            type(cluster_t), pointer :: c

            !> Path helper.
            character(len=:), allocatable :: path_string
            type(varying_string), dimension(:), allocatable :: path

            !> Character variables for the cluster label.
            character(len=256) :: current_line
            character(len=255) :: cluster_label
            character(len=255), dimension(:), allocatable :: cluster_labels

            !> Open the clusterset file and parse it to the cluster_table.
            block
                open(newunit=open_unit, file=clusterset_file, status="old", iostat=ios)
                if (ios /= 0) call pmk_error("could not open '" // trim(clusterset_file) // "'")   
                call toml_parse(cluster_table, open_unit, error)
                close(ios)
                if (allocated(error)) then
                  print '(a)', "Error: "//error%message
                  stop 1
                end if
            end block

            !------------------------------------------------------------------------
            !> Read clusterset data from TOML file
            !------------------------------------------------------------------------

            !> The clusterset input file is opened again, because the names of the clusters
            !  are unknown at this point.

            !>  If the line starts with a '[' and ends with a ']', it is a cluster label.

            !> Count the number of clusters.
            open(newunit=open_unit, file=clusterset_file, status="old", iostat=ios, action="read")
            nr_clusters = 0
            do
                read(open_unit, '(a)', iostat=ios) current_line
                if (ios /= 0) exit 
                
                if (current_line(1:1) == '[') then
                    if (ios /= 0) exit
                    nr_clusters = nr_clusters + 1
                end if
            end do

            if (nr_clusters == 0) call pmk_error("empty clusterset")

            !> Allocate the clusterset array.
            allocate(clusterset(nr_clusters))
            !> Allocare cluster_labels array.
            allocate(cluster_labels(nr_clusters))

            !> Rewind the file.
            rewind open_unit

            nr_clusters = 0
            do
                read(open_unit, '(a)', iostat=ios) current_line
                if (ios /= 0) exit 
                
                if (current_line(1:1) == '[') then
                    if (ios /= 0) exit
                    nr_clusters = nr_clusters + 1
                    if (current_line(1:1) == '[' .and. &
                        current_line(len_trim(current_line):len_trim(current_line)) == ']') then
                        cluster_label = current_line(2:len_trim(current_line)-1)

                        cluster_labels(nr_clusters) = cluster_label
                    end if
                end if 
            end do 
            close(open_unit)

            !> Check for duplicate cluster labels.
            do i = 1, nr_clusters
                do j = i+1, nr_clusters
                    if (cluster_labels(i) == cluster_labels(j)) then
                        call pmk_error("duplicate cluster label '" // trim(cluster_labels(i)) // "'")
                        stop 1
                    end if
                end do
            end do

            !> Read clusterset data from TOML file.
            do i=1, nr_clusters

                        !> Store the cluster in the cluster_t data type.
                        call get_value(cluster_table, cluster_labels(i), child)
                        
                        !> Convert cluster_label to iso_varying_string
                        clusterset(i)%label = cluster_labels(i)

                        !> Assign the c pointer.
                        c => clusterset(i)

                        !> Check for the monomer flag.
                        call get_value(child, "isMonomer", clusterset(i)%monomer, .false.)
                        c%monomer = clusterset(i)%monomer
                        !> If the cluster is a monomer, read the volume.
                        if (clusterset(i)%monomer) then
                            call get_value(child, "volume", clusterset(i)%volume, 0.0_dp)
                            c%volume = clusterset(i)%volume
                        end if

                        !> Get moments of inertia and mass from xyz file.
                        call get_value(child, "coordinates", path_string)
                        if (allocated(path_string)) then
                            write(*,*) len(path_string)
                            allocate(path(1))
                            path(1) = path_string
                            call process_coordinates_record(c, 1, path(1))
                            deallocate(path)
                        else
                            call pmk_missing_key_error("coordinates", clusterset(i)%label)
                        end if

                        !> Get the cluster composition.
                        call get_value(child, "composition", array)
                        if (associated(array)) then
                            allocate(clusterset(i)%composition(len(array)))
                            do j = 1, len(array)
                                call get_value(array, j, clusterset(i)%composition(j))
                                c%composition(j) = clusterset(i)%composition(j)
                            end do
                        else
                            call pmk_missing_key_error("composition", clusterset(i)%label)
                        end if

                        !> Get the frequency scaling factor.
                        call get_value(child, "frequency_scale", clusterset(i)%fscale, 1.0_dp)
                        c%fscale = clusterset(i)%fscale

                        !> Get the frequencies.
                        call get_value(child, "frequencies", path_string)
                        if (allocated(path_string)) then
                            allocate(path(1))
                            path(1) = path_string
                            call process_frequencies_record(c, 1, path(1))
                            deallocate(path)
                        else
                            call pmk_missing_key_error("frequencies", clusterset(i)%label)
                        end if

                        !> Get the anharmonicity constant.
                        call get_value(child, "anharmonicity", clusterset(i)%anharmonicity, 0.0_dp)
                        c%anharmonicity = clusterset(i)%anharmonicity

                        !> Get the interaction energy. 
                        call get_value(child, "energy", clusterset(i)%energy, 0.0_dp)
                        c%energy = clusterset(i)%energy
                        call get_value(child, "energy", array, requested=.false.)
                        if (associated(array)) then
                            if (len(array)==1) then
                                call get_value(array, 1, clusterset(i)%energy)
                                c%energy = clusterset(i)%energy
                            else 
                                call pmk_argument_count_error("energy", c%label)
                            end if
                        end if

                        !> Get the rotational symmetry number sigma.
                        call get_value(child, "sigma", clusterset(i)%sigma, 1)
                        c%sigma = clusterset(i)%sigma
                        call get_value(child, "sigma", array, requested=.false.)
                        if (associated(array)) then
                            if (len(array)==1) then
                                call get_value(array, 1, clusterset(i)%sigma)
                                c%sigma = clusterset(i)%sigma
                            else 
                                call pmk_argument_count_error("sigma", c%label)
                            end if
                        end if

            end do
                                    
            !> Set up the monomer array, assuming that everything is alright. We will check
            !  for errors later.
            allocate(monomer(pmk_input%components))
            monomer = 0
            do i = 1, nr_clusters
                if (clusterset(i)%monomer) then
                    do j = 1, pmk_input%components
                        if (clusterset(i)%composition(j) == 1) monomer(j) = i
                    end do
                end if
            end do

            !> Calculate cluster volumes, assuming that everything is alright. We will
            !  check for errors later.
            do i = 1, nr_clusters
                if (clusterset(i)%monomer) cycle
                clusterset(i)%volume = 0.0_dp
                do j = 1, pmk_input%components
                    clusterset(i)%volume = clusterset(i)%volume + &
                        real(clusterset(i)%composition(j), dp) * &
                        clusterset(monomer(j))%volume
                end do
            end do

        end subroutine process_clusterset

        !=================================================================================
        ! Prints information about the processed cluster set to the screen.
        subroutine print_clusterset()
            ! A pointer to the current cluster for convenience.
            type(cluster_t), pointer :: c
            integer:: i
    
            write(*,'(4X,A)') "Using the following clusterset:"
            write(*,*)
            do i = 1, size(clusterset)
                c => clusterset(i)
    
                ! Print cluster label.
                if (c%monomer) then
                    write(*,'(8X,3A)') "[", trim(char(c%label)), "] (monomer)"
                else
                    write(*,'(8X,3A)') "[", trim(char(c%label)), "]"
                end if
    
                ! Print composition.
                write(*, '(12X,A,1X)', advance = "no") "composition:"
                call array_sample(c%composition)
                write(*, *)
    
                ! Print rotational symmetry number, energy, volume, mass.
                write(*, '(12X,A,1X,G0)') "sigma:", c%sigma
                write(*, '(12X,A,1X,G0.6,1X,A)') "energy:", c%energy, "[kJ/mol]"
                write(*, '(12X,A,1X,G0.6,1X,A)') "volume:", c%volume, "[A^3]"
                write(*, '(12X,A,1X,G0.6,1X,A)') "mass:", c%mass, "[amu]"
    
                ! Print intertia array.
                write(*, '(12X,A,1X,3(G0.6,1X))', advance = "no") "inertia:", &
                    c%inertia(:)
                write(*, '(A)') "[amu*Angstrom^2]"
    
                ! Print anhamronicity constant.
                if (c%anharmonicity > 0.0_dp) then
                    write(*, '(12X,A,1X,G0.6)') "anharmonicity constant", c%anharmonicity
                end if
    
                ! Print frequencies and scaling factor.
                if (abs(c%fscale-1.0_dp) > global_eps) then
                    write(*, '(12X,A,G0.6,A)', advance = "no") &
                        "frequencies (scaled by ", c%fscale, "): "
                else
                    write(*, '(12X,A,1X)', advance = "no") "frequencies:"
                end if
                call array_sample(c%frequencies)
                write(*, '(1X,A)') "[1/cm]"
                write(*, *)
            end do
        end subroutine print_clusterset

        !=================================================================================
        ! Processes a coordinates record.
        ! The coordinates record is a file name. The file contains the number of atoms in
        ! the first line, followed by an empty/comment line, followed by one line per atom.
        ! Each line contains the atom label and the x, y, and z coordinates in Angstrom.
        ! The total mass of the cluster is calculated and the origin is shifted to the
        ! center of mass. The inertia tensor is calculated and diagonalized. The cluster
        ! is classified as atom, linear, or non-linear. The moments of inertia are assigned
        ! accordingly.
        subroutine process_coordinates_record(c, nr_args, args)
            use atomic_data, only: periodic_table
            use constants
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: i
            integer:: ios
            integer:: my_unit
            integer:: nr_atoms
    
            character(2), dimension(:), allocatable :: label
            real(dp), dimension(:), allocatable :: mass
            real(dp), dimension(:, :), allocatable :: xyz
            real(dp), dimension(3) :: com
            real(dp), dimension(3, 3) :: inertia
            real(dp), dimension(3) :: eig
            integer:: n
    
            if (nr_args == 1) then
                ! Open unit.
                open(newunit = my_unit, file = char(args(1)), action = 'read', &
                    status = 'old', iostat = ios)
                if (ios /= 0) call pmk_error("could not open '" // char(args(1)) // "'")
    
                ! Read number of atoms.
                read(my_unit, *, iostat = ios) nr_atoms
                if (ios /= 0) call pmk_error("illegal file format in '" // &
                    char(args(1)) // "'")
    
                ! Comment line.
                read(my_unit, *, iostat = ios)
    
                ! Read coordinates.
                allocate(label(nr_atoms))
                allocate(mass(nr_atoms))
                allocate(xyz(nr_atoms, 3))
                do i = 1, nr_atoms
                    read(my_unit, *, iostat = ios) &
                        label(i), xyz(i, 1), xyz(i, 2), xyz(i, 3)
                    if (ios /= 0) call pmk_error("unexpected end of file in '" // &
                        char(args(1)) // "'")
                end do
    
                ! Assign masses and calculate total mass.
                do i = 1, nr_atoms
                    mass(i) = periodic_table%mass(label(i))
                end do
                c%mass = sum(mass)
    
                ! Calculate center of mass and shift to origin.
                call center_of_mass(nr_atoms, com, mass, xyz)
    
                ! Calculate inertia tensor.
                inertia = 0.0_dp
                do i = 1, nr_atoms
                    inertia(1,1) = inertia(1,1) + mass(i)*(xyz(i,2)**2 + xyz(i,3)**2)
                    inertia(2,2) = inertia(2,2) + mass(i)*(xyz(i,1)**2 + xyz(i,3)**2)
                    inertia(3,3) = inertia(3,3) + mass(i)*(xyz(i,1)**2 + xyz(i,2)**2)
                    inertia(1,2) = inertia(1,2) - mass(i)*xyz(i,1)*xyz(i,2)
                    inertia(1,3) = inertia(1,3) - mass(i)*xyz(i,1)*xyz(i,3)
                    inertia(2,3) = inertia(2,3) - mass(i)*xyz(i,2)*xyz(i,3)
                end do
                inertia(2, 1) = inertia(1, 2)
                inertia(3, 1) = inertia(1, 3)
                inertia(3, 2) = inertia(2, 3)
    
                ! Diagonalize inertia tensor.
                call diagonalize_3x3(inertia, eig)
    
                ! Check for atoms, linear molecules, and assign moments of inertia.
                n = 3 - count(eig <= global_eps)
                if (n == 3) then
                    allocate(c%inertia(3))
                    c%inertia = eig
                else if (n == 2) then
                    c%linear = .true.
                    allocate(c%inertia(1))
                    c%inertia(1) = eig(1)
                else if (n == 0) then
                    c%atom = .true.
                end if
    
                ! Clean up.
                close(my_unit)
                deallocate(label)
                deallocate(mass)
                deallocate(xyz)
            else
                call pmk_argument_count_error("coordinates", c%label)
            end if
        end subroutine process_coordinates_record

        !=================================================================================
        ! Calculation of the center of mass of a cluster.
        ! Shifts the origin to the center of mass.
        subroutine center_of_mass(nr_atoms, com, mass, xyz)
            !> Number of atoms.
            integer, intent(in) :: nr_atoms
            !> Center of mass.
            real(dp), dimension(3), intent(out) :: com
            !> Masses of the atoms.
            real(dp), dimension(:), intent(in) :: mass
            !> Coordinates of the atoms.
            real(dp), dimension(:, :), intent(inout) :: xyz

            integer:: i

            com = 0.0_dp
            ! Calculate center of mass.
            do i = 1, nr_atoms
                com(:) = com(:) + mass(i)*xyz(i, :)
            end do
            com(:) = com(:) / sum(mass)

            ! Shift to origin.
            do i = 1, nr_atoms
                xyz(i, :) = xyz(i, :) - com(:)
            end do
        end subroutine center_of_mass

        !=================================================================================
        ! Direct method for the diagonalization of a 3x3 symmetric matrix.
        subroutine diagonalize_3x3(a, eig)
            !> Symmetric 3x3 matrix.
            real(dp), dimension(3, 3), intent(in) :: a
            !> Eigenvalues.
            real(dp), dimension(3), intent(out) :: eig

            real(dp):: p1, p2, p, q, r, phi, tmp
            real(dp), dimension(3, 3) :: identity
            real(dp), dimension(3, 3) :: B

            !> Definition of 3x3 identity matrix.
            identity = 0.0_dp
            identity(1,1) = 1.0_dp
            identity(2,2) = 1.0_dp
            identity(3,3) = 1.0_dp

            p1 = a(1,2)**2 + a(1,3)**2 + a(2,3)**2
            if (p1 <= global_eps) then
                ! Matrix is already diagonal.
                eig(1) = a(1,1)
                eig(2) = a(2,2)
                eig(3) = a(3,3)
            else
                q = (a(1,1) + a(2,2) + a(3,3))/3.0_dp
                p2 = (a(1,1)-q)**2 + (a(2,2)-q)**2 + &
                    (a(3,3)-q)**2 + 2.0_dp*p1
                p = sqrt(p2/6.0_dp)
                B = (a(:, :) - q*identity(:, :))/p
                r = B(1,1)*B(2,2)*B(3,3) + B(1,2)*B(2,3)*B(3,1) + &
                    B(1,3)*B(2,1)*B(3,2) - B(1,3)*B(2,2)*B(3,1) - &
                    B(1,2)*B(2,1)*B(3,3) - B(1,1)*B(2,3)*B(3,2)
                r = 0.5_dp * r

                if (r <= -1.0_dp) then
                    phi = pi/3.0_dp
                else if (r >= 1.0_dp) then
                    phi = 0.0_dp
                else
                    phi = acos(r)/3.0_dp
                end if

                eig(1) = q + 2.0_dp*p*cos(phi)
                eig(3) = q + 2.0_dp*p*cos(phi + (2.0_dp*pi/3.0_dp))
                eig(2) = 3.0_dp*q - eig(1) - eig(3)
            end if

            ! Sort the eigenvalues.
            if (eig(1) < eig(2)) then
                tmp = eig(2)
                eig(2) = eig(1)
                eig(1) = tmp
            end if
            if (eig(1) < eig(3)) then
                tmp = eig(3)
                eig(3) = eig(1)
                eig(1) = tmp
            end if
            if (eig(2) < eig(3)) then
                tmp = eig(3)
                eig(3) = eig(2)
                eig(2) = tmp
            end if

        end subroutine diagonalize_3x3


        !=================================================================================
        ! Reads frequencies from a frequency file. A frequency file is similar to an
        ! xyz file. It contains the number of vibrational frequencies in the first line,
        ! followed by an empty/comment line, followed by one line per frequency.
        ! Sorts the frequencies.
        subroutine process_frequencies_record(c, nr_args, args)
            type(cluster_t), pointer, intent(inout) :: c
            integer, intent(in) :: nr_args
            type(varying_string), dimension(nr_args), intent(in) :: args
    
            integer:: nr_frequencies
            integer:: nr_zeros
            integer:: my_unit
            integer:: ios
            integer:: i
            integer:: j
            real(dp):: tmp
    
            if (nr_args == 1) then
                ! Open unit.
                open(newunit = my_unit, file = char(args(1)), action = 'read', &
                    status = 'old', iostat = ios)
                if (ios /= 0) call pmk_error("could not open '" // char(args(1)) // "'")
    
                ! Read number of frequencies.
                read(my_unit, *, iostat = ios) nr_frequencies
                if (ios /= 0) call pmk_error("illegal file format in '" // &
                    char(args(1)) // "'")
    
                ! Comment line.
                read(my_unit, *, iostat = ios)
    
                ! Read number of zeros.
                nr_zeros = 0
                do i = 1, nr_frequencies
                    read(my_unit, *, iostat = ios) tmp
                    if (ios /= 0) call pmk_error("unexpected end of file in '" // &
                        char(args(1)) // "'")
                    if (abs(tmp) <= global_eps) then
                        nr_zeros = nr_zeros + 1
                    end if
                end do
                rewind(my_unit)
                allocate(c%frequencies(nr_frequencies - nr_zeros))
    
                ! Read frequencies.
                nr_zeros = 0
                read(my_unit, *, iostat = ios)
                read(my_unit, *, iostat = ios)
                do i = 1, nr_frequencies
                    read(my_unit, *, iostat = ios) tmp
                    if (abs(tmp) <= global_eps) then
                        nr_zeros = nr_zeros + 1
                    else
                        c%frequencies(i - nr_zeros) = tmp
                    end if
                end do
    
                ! Sort frequencies.
                do i = 1, size(c%frequencies) - 1
                    do j = i + 1, size(c%frequencies)
                        if (c%frequencies(j) < c%frequencies(i)) then
                            tmp = c%frequencies(i)
                            c%frequencies(i) = c%frequencies(j)
                            c%frequencies(j) = tmp
                        end if
                    end do
                end do
    
                c%frequencies = c%fscale*c%frequencies
                close(my_unit)
            else
                call pmk_argument_count_error("frequencies", c%label)
            end if
        end subroutine process_frequencies_record

        !=================================================================================
        ! Performs sanity checks on the clusterset.
        subroutine check_clusterset()
            ! A pointer to the current cluster.
            type(cluster_t), pointer :: c
            ! Loop and status stuff.
            integer:: i
            integer:: j
    
            ! Check monomer count
            if (count(clusterset%monomer) /= pmk_input%components) &
                call pmk_error("invalid number of monomers")
    
            do i = 1, size(clusterset)
                c => clusterset(i)
    
                ! Check compositions
                if (any(c%composition < 0)) &
                    call pmk_unphysical_argument_error("composition", c%label)
    
                ! Check sigma
                if (c%sigma < 0) &
                    call pmk_unphysical_argument_error("sigma", c%label)
    
                ! Check the anharmoncity constant.
                if (c%anharmonicity < 0.0_dp) &
                    call pmk_unphysical_argument_error("anharmonicity", c%label)
    
                ! Check the frequency scaling factor.
                if (c%fscale <= 0.0_dp) &
                    call pmk_unphysical_argument_error("frequency_scale", c%label)
    
                ! Check frequencies
                do j = 1, size(c%frequencies)
                    if (c%frequencies(j) < 0.0_dp) &
                        call pmk_unphysical_argument_error("frequencies", c%label)
                end do
    
                ! Check volume
                if (c%volume <= 0.0_dp) &
                    call pmk_unphysical_argument_error("volume", c%label)
    
                ! Check monomers.
                if (c%monomer) then
    
                    ! Check sum of composition.
                    if (sum(c%composition) /= 1) call pmk_error(&
                        "monomer/composition mismatch in cluster '[" // &
                        char(c%label) // "]'")
    
                    ! Check that this is the only monomer for the given component, by
                    ! comparing it to the monomer array. Assume that there are multiple
                    ! monomers for a particular component. Then because of the way we set
                    ! up the monomer array in setup_clusterset(), the monomer array will
                    ! point to the last 'monomer' for this component. We can use this fact
                    ! to check for multiply defined monomers.
                    do j = 1, pmk_input%components
                        ! Check for missing monomer first
                        if (monomer(j) == 0) call pmk_error("missing monomers")
                        if (c%composition(j) == 1) then
                            if (monomer(j) /= i) call pmk_error( &
                                "more than one monomer per component")
                        end if
                    end do
                end if
            end do
        end subroutine check_clusterset
        !=================================================================================
end module cluster
!=========================================================================================

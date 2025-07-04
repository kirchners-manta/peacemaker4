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
! This module implements the partition functions used by Peacemaker.
module partition_functions
    use kinds
    use cluster
    use constants
    use shared_data
    implicit none
    private
    !=====================================================================================
    ! Public entities
    public :: calculate_lnq
    public :: calculate_lnqtrans, calculate_lnqvib, calculate_lnqrot, calculate_lnqelec, &
              calculate_lnqint
    public :: update_lnq
    public :: calculate_dlnq
    public :: calculate_dlnqtrans, calculate_dlnqvib, calculate_dlnqrot, calculate_dlnqelec, &
              calculate_dlnqint
    public :: calculate_ddlnq
    public :: calculate_ddlnqtrans, calculate_ddlnqvib, calculate_ddlnqrot, calculate_ddlnqelec, &
              calculate_ddlnqint
    public :: pf_t
    !=====================================================================================
    ! Data type representing a cluster partition function.
    type :: pf_t
        real(dp) :: qtrans, qvib, qrot, qelec, qint, qtot
    end type pf_t
    !=====================================================================================
    contains
        !=================================================================================
        ! Calculates all cluster partition functions.
        subroutine calculate_lnq(lnq, amf, bxv, temp, vol)
            type(pf_t), dimension(size(clusterset)), intent(out) :: lnq
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: temp
    
            call calculate_lnqtrans(lnq(:)%qtrans, bxv, temp, vol, clusterset, global_data%vexcl)
            call calculate_lnqvib(lnq(:)%qvib, temp, clusterset, global_data%rotor_cutoff)
            call calculate_lnqrot(lnq(:)%qrot, temp, clusterset)
            call calculate_lnqelec(lnq(:)%qelec, temp, clusterset)
            call calculate_lnqint(lnq(:)%qint, amf, temp, vol, clusterset, global_data%ntot)
            lnq(:)%qtot = lnq(:)%qtrans + lnq(:)%qvib + lnq(:)%qrot + lnq(:)%qelec + &
                lnq(:)%qint
        end subroutine calculate_lnq
        !=================================================================================
        ! Updates the cluster partition functions that depend on the volume.
        subroutine update_lnq(lnq, amf, bxv, temp, vol)
            type(pf_t), dimension(size(clusterset)), intent(inout) :: lnq
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: vol
    
            call calculate_lnqtrans(lnq(:)%qtrans, bxv, temp, vol, clusterset, global_data%vexcl)
            call calculate_lnqint(lnq(:)%qint, amf, temp, vol, clusterset, global_data%ntot)
            lnq(:)%qtot = lnq(:)%qtrans + lnq(:)%qvib + lnq(:)%qrot + lnq(:)%qelec + &
                lnq(:)%qint
        end subroutine update_lnq
        !=================================================================================
        ! Calculates analytical derivatives of the cluster partition functions.
        subroutine calculate_dlnq(dlnq, amf, bxv, amf_temp, bxv_temp, temp, vol)
            type(pf_t), dimension(size(clusterset)), intent(out) :: dlnq
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: amf_temp
            real(dp), intent(in) :: bxv_temp
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: vol
    
            call calculate_dlnqtrans(dlnq(:)%qtrans, bxv, bxv_temp, temp, vol, global_data%vexcl)
            call calculate_dlnqvib(dlnq(:)%qvib, temp, clusterset, global_data%rotor_cutoff)
            call calculate_dlnqrot(dlnq(:)%qrot, temp, clusterset)
            call calculate_dlnqelec(dlnq(:)%qelec, temp, clusterset)
            call calculate_dlnqint(dlnq(:)%qint, temp, amf, amf_temp, vol, clusterset, global_data%ntot)
            dlnq(:)%qtot = dlnq(:)%qtrans + dlnq(:)%qvib + dlnq(:)%qrot + &
                dlnq(:)%qelec + dlnq(:)%qint
        end subroutine calculate_dlnq
        !=================================================================================
        ! Calculates analytical curvatures of the cluster partition functions.
        subroutine calculate_ddlnq(ddlnq, amf, bxv, amf_temp, bxv_temp, temp, vol)
            type(pf_t), dimension(size(clusterset)), intent(out) :: ddlnq
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: amf_temp
            real(dp), intent(in) :: bxv_temp
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: vol
    
            call calculate_ddlnqtrans(ddlnq(:)%qtrans, bxv, bxv_temp, temp, vol, global_data%vexcl)
            call calculate_ddlnqvib(ddlnq(:)%qvib, temp, clusterset, global_data%rotor_cutoff)
            call calculate_ddlnqrot(ddlnq(:)%qrot, temp, clusterset)
            call calculate_ddlnqelec(ddlnq(:)%qelec, temp, clusterset)
            call calculate_ddlnqint(ddlnq(:)%qint, temp, amf, amf_temp, vol, clusterset, global_data%ntot)
            ddlnq(:)%qtot = ddlnq(:)%qtrans + ddlnq(:)%qvib + ddlnq(:)%qrot + &
                ddlnq(:)%qelec + ddlnq(:)%qint
        end subroutine calculate_ddlnq
        !=================================================================================
        ! Calculates the translational cluster partition function.
        ! q_trans = (2*pi*m*kb*T/h^2)^1.5*V
        subroutine calculate_lnqtrans(lnq, bxv, temp, vol, cluster_set, v_excl)
            real(dp), dimension(:), intent(out) :: lnq
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: temp
            real(dp), intent(in) :: bxv
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), intent(in) :: v_excl
    
            real(dp):: lambda
            real(dp):: mass
            integer:: iclust
    
            do iclust = 1, size(cluster_set)
                ! Calculate the de Broglie wave length
                mass = cluster_set(iclust)%mass*amu
                lambda = planck/sqrt(2.0_dp*pi*mass*kb*temp)
    
                ! Calculate the partition function
                lnq(iclust) = log(vol-bxv*v_excl) - 3.0_dp*log(lambda)
            end do
            
            ! If a value is infinite, warn the user.
            if (any(lnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The translational partition function of one cluster'
                write(*,*) '         is infinite.'
            end if

        end subroutine calculate_lnqtrans
        !=================================================================================
        ! Calculates the vibrational cluster partition function.
        ! edit: q_vib is the product, the sum occurs due to the logarithm
        ! q_vib = product over all nu: exp[-h*nu/(2*kb*T)]/{1-exp[-h*nu/(kb*T)]}
        subroutine calculate_lnqvib(lnq, temp, cluster_set, rotor_cutoff)
            real(dp), dimension(:), intent(out) :: lnq
            real(dp), intent(in) :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), intent(in) :: rotor_cutoff
    
            integer:: iclust
            integer:: ifreq, imoment
            real(dp):: factor
            real(dp):: t_vib
            real(dp):: q_ho, q_fr
            real(dp):: moi, emoi, Bav, t_rot, w
            logical :: free_rotator
            
            free_rotator = .false.
            w = 1.0_dp
            if (rotor_cutoff > 0.0_dp) free_rotator = .true.

            ! TODO: Check atoms.
            factor = planck*100.0_dp*speed_of_light/kb
            do iclust = 1, size(cluster_set)
                associate(f => cluster_set(iclust)%frequencies, q => lnq(iclust), &
                    x => cluster_set(iclust)%anharmonicity, c => cluster_set(iclust))
                    q = 0.0_dp
                    q_ho = 0.0_dp
                    q_fr = 0.0_dp
                    
                    if (c%atom) cycle
                    
                    ! Average intertia moment of cluster for free rotator approximation
                    Bav = sum(c%inertia)/real(size(c%inertia,1),dp)*amu*1.0e-20_dp
                    
                    do ifreq = 1, size(f)
                        t_vib = factor*f(ifreq)
                        if (x > 0.0_dp) then
                            ! Morse oscillator.
                            q = q - log(2.0_dp*sinh(0.5_dp*t_vib/temp)) + x*t_vib/temp * &
                                (0.25_dp + 0.5_dp/sinh(0.5_dp*t_vib/temp)**2)
                        else
                            ! Harmonic oscillator.
                            ! TODO: Use log1p and exp1m family of functions.
                            q_ho = -0.5_dp*t_vib/temp - log(1.0_dp - exp(-t_vib/temp))

                            if (free_rotator) then
                                ! Free rotator approximation by Grimme.
                                moi = planck/(8 * pi**2 * 100.0_dp * speed_of_light * f(ifreq))
                                emoi = moi * Bav/(moi + Bav)
                                t_rot = hbar**2/(2.0_dp*emoi*kb)
                                q_fr = log(sqrt(pi*temp/t_rot)/real(c%sigma, dp))
                                w = 1.0_dp/(1.0_dp + (rotor_cutoff/f(ifreq))**4)
                            end if
                            
                            q = q + w * q_ho + (1.0_dp - w) * q_fr
                        end if
                    end do
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(lnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The vibrational partition function of one cluster'
                write(*,*) '         is infinite.'
            end if

        end subroutine calculate_lnqvib
        !=================================================================================
        ! Calculates the rotational cluster partition function.
        ! q_rot = (1/sigma)*(T/t_rot)                               --- linear
        ! q_rot = sqrt(pi)/sigma*(T^3/(t_rot1*t_rot2*t_rot3))^(1/2) --- nonlinear
        ! q_rot = 1                                                 --- atom
        subroutine calculate_lnqrot(lnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: lnq
            real(dp), intent(in) :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
    
            integer:: iclust
            integer:: imoment
            real(dp):: t_rot
    
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust))
                    if (c%atom) then
                        lnq(iclust) = 0.0_dp
                    else
                        lnq(iclust) = 1.0_dp
                        do imoment = 1, size(c%inertia)
                            t_rot = hbar**2/(2.0_dp*c%inertia(imoment)*amu*1.0e-20_dp*kb)
                            lnq(iclust) = lnq(iclust) * temp/t_rot
                        end do
                        if (c%linear) then
                            ! linear cluster:
                            ! q_rot = T/(sigma*t_rot)
                            ! Note: For a linear molecule, the inertia array contains
                            ! two moments of inertia of equal magnitude, which were
                            ! multiplied together above.
                            lnq(iclust) = sqrt(lnq(iclust))/real(c%sigma, dp)
                        else
                            ! polyatomic, nonlinear cluster:
                            ! q_rot = sqrt(pi)/sigma*(T^3/(t_rot1*t_rot2*t_rot3))
                            lnq(iclust) = sqrt(pi*lnq(iclust))/real(c%sigma, dp)
                        end if
                        ! Let's not forget the logarithm.
                        lnq(iclust) = log(lnq(iclust))
                    end if
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(lnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The rotational partition function of one cluster'
                write(*,*) '         is infinite.'
            end if

        end subroutine calculate_lnqrot
        !=================================================================================
        ! Calculates the electronic cluster partition function.
        ! q_elec = exp(-dE/(kb*T))
        subroutine calculate_lnqelec(lnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: lnq
            real(dp), intent(in) :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
    
            ! The adiabatic interaction energy is in units of kJ/mol.
            lnq(:) = (-1000.0_dp/avogadro)*cluster_set(:)%energy/(kb*temp)

            ! If a value is infinite, warn the user.
            if (any(lnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The electronic partition function of one cluster'
                write(*,*) '         is infinite.'
            end if

        end subroutine calculate_lnqelec
        !=================================================================================
        ! Calculates the mean field partition function.
        subroutine calculate_lnqint(lnq, amf, temp, vol, cluster_set, ntot)
            real(dp), dimension(:), intent(out) :: lnq
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: amf
            real(dp), intent(in) :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), dimension(:), intent(in) :: ntot
    
            integer:: iclust
            real(dp):: emf
            
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust))
                    emf = -amf*real(sum(c%composition), dp)*sum(ntot)/vol
                    lnq(iclust) = -emf / (kb*temp)
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(lnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The mean field partition function of one cluster'
                write(*,*) '         is infinite.'
            end if

        end subroutine calculate_lnqint
        !=================================================================================
        ! Calculates the temperature derivative of the translational partition function.
        subroutine calculate_dlnqtrans(dlnq, bxv, bxv_temp, temp, vol, vexcl)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: bxv_temp
            real(dp), intent(in)  :: temp
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: vexcl
    
            ! Fails if -vol+bxv*vexcl is zero.
            dlnq(:) = 1.5_dp/temp + (vexcl * bxv_temp)/(-vol+bxv*vexcl)

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The temperature derivative of the translational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if


        end subroutine calculate_dlnqtrans
        !=================================================================================
        ! Calculates the temperature derivative of the vibrational partition function.
        subroutine calculate_dlnqvib(dlnq, temp, cluster_set, rotor_cutoff)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), intent(in) :: rotor_cutoff
    
            integer :: iclust
            integer:: ifreq
            real(dp):: factor
            real(dp):: t_vib
            real(dp):: fsinh
            real(dp):: fcosh
            real(dp):: q_ho, q_fr
            real(dp):: w
    
            ! TODO: Check atoms.
            factor = planck*100.0_dp*speed_of_light/kb
            do iclust = 1, size(cluster_set)
                associate(f => cluster_set(iclust)%frequencies, q => dlnq(iclust), &
                        x => cluster_set(iclust)%anharmonicity)
                    q = 0.0_dp
                    q_ho = 0.0_dp
                    q_fr = 0.0_dp
                    
                    do ifreq = 1, size(f)
                        t_vib = factor*f(ifreq)
                        if (x > 0.0_dp) then
                            ! Morse oscillator.
                            fsinh = sinh(0.5_dp*t_vib/temp)
                            fcosh = cosh(0.5_dp*t_vib/temp)
                            q = q + 0.5_dp*t_vib*fcosh/(temp**2*fsinh) - &
                                x*t_vib/temp**2*(0.25_dp + 0.5_dp/fsinh**2) + &
                                0.5_dp*x*t_vib**2*fcosh/(temp**3*fsinh**3)
                        else
                            ! Harmonic oscillator.
                            ! TODO: Use log1p and exp1m family of functions.
                            q_ho = 0.5_dp*t_vib/temp**2 + &
                                t_vib/temp**2/(exp(t_vib/temp)-1.0_dp)

                            ! Free rotator approximation by Grimme.
                            q_fr = 0.5_dp/temp
                            
                            w = 1.0_dp/(1.0_dp + (rotor_cutoff/f(ifreq))**4)
                            q = q + w * q_ho + (1.0_dp - w) * q_fr
                        end if
                    end do
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The temperature derivative of the vibrational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_dlnqvib
        !=================================================================================
        ! Calculates the temperature derivative of the rotational partition function.
        subroutine calculate_dlnqrot(dlnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set

            integer :: iclust
    
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust), q => dlnq(iclust))
                    if (c%atom) then
                        q = 0.0_dp
                    else if (c%linear) then
                        ! linear cluster:
                        ! q_rot = T/(sigma*t_rot) -> d/dt(lnq)=1/T
                        q = 1.0_dp/temp
                    else
                        ! polyatomic, nonlinear cluster: d/dT(lnq)=3/(2T)
                        q = 1.5_dp/temp
                    end if
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The temperature derivative of the rotational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_dlnqrot
        !=================================================================================
        ! Calculates the temperature derivative of the electronic partition function.
        subroutine calculate_dlnqelec(dlnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
    
            ! The adiabatic interaction energy is in units of kJ/mol.
            dlnq(:) = (1000.0_dp/avogadro)*cluster_set(:)%energy/(kb*temp**2)

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The temperature derivative of the electronic'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_dlnqelec
        !=================================================================================
        ! Calculates the temperature derivative of the mean field partition function.
        subroutine calculate_dlnqint(dlnq, temp, amf, amf_temp, vol, cluster_set, ntot)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            real(dp), intent(in)  :: amf
            real(dp), intent(in)  :: amf_temp
            real(dp), intent(in)  :: vol
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), dimension(:), intent(in) :: ntot
    
            integer :: iclust
            real(dp):: emf
    
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust))
                    emf = -real(sum(c%composition), dp)*sum(ntot)/vol * (amf - temp * amf_temp)
                    dlnq(iclust) = emf / (kb*temp**2)
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The temperature derivative of the mean field'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_dlnqint
        !=================================================================================
        ! Calculates the second temperature derivative of the translational partition
        ! function.
        subroutine calculate_ddlnqtrans(dlnq, bxv, bxv_temp, temp, vol, vexcl)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in) :: bxv
            real(dp), intent(in) :: bxv_temp
            real(dp), intent(in)  :: temp
            real(dp), intent(in) :: vol
            real(dp), intent(in) :: vexcl

            dlnq(:) = -1.5_dp/temp**2 - &
                     (vexcl * bxv_temp)**2/(-vol+bxv*vexcl)**2

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The second temperature derivative of the translational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_ddlnqtrans
        !=================================================================================
        ! Calculates the second temperature derivative of the vibrational partition
        ! function.
        subroutine calculate_ddlnqvib(dlnq, temp, cluster_set, rotor_cutoff)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), intent(in) :: rotor_cutoff
    
            integer :: iclust
            integer:: ifreq
            real(dp):: factor
            real(dp):: t_vib
            real(dp):: fsinh
            real(dp):: fcosh
            real(dp):: q_ho, q_fr
            real(dp):: w
    
            ! TODO: Check atoms.
            factor = planck*100.0_dp*speed_of_light/kb
            do iclust = 1, size(cluster_set)
                associate(f => cluster_set(iclust)%frequencies, q => dlnq(iclust), &
                        x => cluster_set(iclust)%anharmonicity)
                    q = 0.0_dp
                    q_ho = 0.0_dp
                    q_fr = 0.0_dp
                    
                    do ifreq = 1, size(f)
                        t_vib = factor*f(ifreq)
                        ! Check if anharmonicity constant is set
                        if (x > 0.0_dp) then
                            ! Morse oscillator.
                            fsinh = sinh(0.5_dp*t_vib/temp)
                            fcosh = cosh(0.5_dp*t_vib/temp)
                            q = q - fcosh/fsinh*t_vib/(temp**3) - &
                                2.0_dp*x*t_vib**2/(temp**4)*fcosh/fsinh**3 + &
                                (t_vib/(2.0_dp*temp**2*fsinh))**2 + &
                                2.0_dp*x*t_vib/(temp**3)*(0.25_dp + 0.5_dp/fsinh**2) + &
                                x*(t_vib**3)/(4.0_dp*temp**5)*(2.0_dp*fcosh**2+1.0_dp) / &
                                fsinh**4
                        else
                            ! Harmonic oscillator.
                            ! TODO: Use log1p and exp1m family of functions.
                            q_ho = -t_vib/(temp**3) &
                                - 2.0_dp*t_vib/(temp**3*(exp(t_vib/temp)-1.0_dp)) &
                                + t_vib**2*exp(t_vib/temp) / &
                                (temp**4*(exp(t_vib/temp)-1)**2)

                            ! Free rotator approximation by Grimme.
                            q_fr = -0.5_dp/temp**2
                            
                            w = 1.0_dp/(1.0_dp + (rotor_cutoff/f(ifreq))**4)
                            q = q + w * q_ho + (1.0_dp - w) * q_fr
                        end if
                    end do
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The second temperature derivative of the vibrational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_ddlnqvib
        !=================================================================================
        ! Calculates the second temperature derivative of the rotational partition
        ! function.
        subroutine calculate_ddlnqrot(dlnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set

            integer :: iclust
    
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust), q => dlnq(iclust))
                    if (c%atom) then
                        q = 0.0_dp
                    else if (c%linear) then
                        ! linear cluster:
                        ! q_rot = T/(sigma*t_rot) -> d/dt(lnq)=-1/T**2
                        q = -1.0_dp/temp**2
                    else
                        ! polyatomic, nonlinear cluster: d/dT(lnq)=-3/(2T**2)
                        q = -1.5_dp/temp**2
                    end if
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The second temperature derivative of the rotational'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_ddlnqrot
        !=================================================================================
        ! Calculates the second temperature derivative of the electronic partition
        ! function.
        subroutine calculate_ddlnqelec(dlnq, temp, cluster_set)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            type(cluster_t), dimension(:), intent(in) :: cluster_set
    
            ! The adiabatic interaction energy is in units of kJ/mol.
            dlnq(:) = - (2000.0_dp/avogadro)*cluster_set(:)%energy/(kb*temp**3)

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The second temperature derivative of the electronic'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_ddlnqelec
        !=================================================================================
        ! Calculates the second temperature derivative of the the mean field partition
        ! function.
        subroutine calculate_ddlnqint(dlnq, temp, amf, amf_temp, vol, cluster_set, ntot)
            real(dp), dimension(:), intent(out) :: dlnq
            real(dp), intent(in)  :: temp
            real(dp), intent(in)  :: amf
            real(dp), intent(in)  :: amf_temp
            real(dp), intent(in)  :: vol
            type(cluster_t), dimension(:), intent(in) :: cluster_set
            real(dp), dimension(:), intent(in) :: ntot
    
            integer :: iclust
            real(dp):: emf
    
            do iclust = 1, size(cluster_set)
                associate(c => cluster_set(iclust))
                    emf = -real(sum(c%composition), dp)*sum(ntot)/vol * &
                         (amf - temp * amf_temp)
                    dlnq(iclust) = -2.0_dp*emf/(kb*temp**3)
                end associate
            end do

            ! If a value is infinite, warn the user.
            if (any(dlnq >= huge(0.0_dp))) then
                write(*,*) 'Warning: The second temperature derivative of the mean field'
                write(*,*) '         partition function of one cluster is infinite.'
            end if

        end subroutine calculate_ddlnqint
        !=================================================================================
end module partition_functions
!=========================================================================================

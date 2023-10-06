! *****************************************************************************
!
! File:				swpcoag_model.f90
! Project:			Sweep 2
! Author(s):			Matthew Celnik (msc37) & Rob Patterson (riap2)
!
! Copyright (C) 2006  Matthew S Celnik & Robert Patterson
!
! Licence:
!   This file is part of "sweep2".
!
!   sweep2 is free software; you can redistribute it and/or
!   modify it under the terms of the GNU General Public License
!   as published by the Free Software Foundation; either version 2
!   of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!
! Contact:
!   Dr Markus Kraft
!   Dept of Chemical Engineering
!   University of Cambridge
!   New Museums Site
!   Pembroke Street
!   Cambridge
!   CB2 3RA
!   UK
!
!   Email:   mk306@cam.ac.uk
!   Website: como.cheng.cam.ac.uk
!
! Purpose:
!   Defintions of different coagulation kernel constants and function to
!   calculate kernel values.
!
! Functions:
!   CoagKernel  		-	Calculates the coagualation kernel between two
!						stochastic particles for the given temperature
!						and pressure.
!   FreeMolKernel		-	Calculates the free-molecular coagulation
!						kernel between two stochastic particles.
!   SlipFlowKernel		-	Calculates the slip-flow coagulation kernel
!						between two stochastic particles.
!   HalfHarmonicMean		-	Returns half the harmonic mean of two numbers,
!						this is a helper function for combining
!						coagulation kernels.
! -----------------------------------------------------------------------------
!   NucFreeMolK 		-	Calculates the constant in the free-molecular
!						coagulation kernel for two particles of
!						different sizes.  To get the actual kernel value
!						this function must be multiplied by Sprt(T).
!   NucSlipFlowK		-	Calculates the constants in the slip-flow
!						coagulation kernel for two particles of different
!						sizes.  See function description for notes
!						on how to obtain the true kernel value.
! -----------------------------------------------------------------------------
!   DimerFreeMolK		-	Calculates the constant in the free-molecular
!						coagulation kernel for two particles of the
!						same size.  To get the actual kernel value
!						this function needs to be multiplied by Sqrt(T).
!   DimerSlipFlowK		-	Calculates the constants in the slip-flow
!						coagulation kernel for two particles of the
!						same size.  See function description for notes
!						on how to obtain the true kernel value.
! -----------------------------------------------------------------------------
!   FreeMolCondTerms		-	Calculates the condensation terms for a species
!						of size x condensing on any particle for the
!						free-molecular regime.  These terms are only
!						valid if the species is smaller than the particle.
! *****************************************************************************

Module SWPCOAG_MODEL
    Implicit None
    Public

    ! -------------------------------------------------------
    ! COAGULATION KERNEL CONSTANTS.
    ! -------------------------------------------------------

    ! Free-molecular enhancement factor.
    Real, Parameter, Public  :: E_fm      = 2.2E0 
    Real, Parameter, Public  :: E_cond    = 1.3E0  !1.3E0 !2.2E0 
    Real, Parameter, Public  :: E_nuc     = 2.5E0  !2.5E0 !2.2E0 
    Real, Parameter, Public  :: agg_eff   = 1.0E0
    Real, Parameter, Public  :: cond_eff  = 1.0E0  !1.0E0 !0.3E+0
    Real, Parameter, Public  :: nuc_eff   = 1.0E0  !1.0E0 !5.0E-04
    Real, Parameter, Public  :: dimer_eff = 1.5E-02
    ! Free-molecular kernel constant (need to multiply by T^1/2).
    Real, Parameter, Public  :: C_fm = 1.47265680E-08

    ! Free-molecular majorant factor for coagulation kernel .
    Real, Parameter, Public  :: C_fmmaj = 1.4178E0

    ! Slip-flow kernel constant.
    Real, Parameter, Public  :: C_sf = 9.2046667E-17

    ! -------------------------------------------------------
    ! COAGULATION PARAMETERS.
    ! -------------------------------------------------------

    ! Coagulation majorant enumeration.
    Integer, Parameter, Public :: NOMAJ=0, SLIPMAJ=1, FREEMAJ=2

    Contains
    ! -------------------------------------------------------
    ! KERNEL CONSTANTS ROUTINES.
    !
    !   These routines calculate constant values in coagulation
    !   terms for given particle sizes.
    !
    ! -------------------------------------------------------

    Function CoagKernel(sp1, sp2, T, P, maj) Result(k)
        ! DESCRIPTION:
        !   Calculates the coagulation kernel between
        !   two stochastic particles for the given
        !   temperature and pressure.
        ! RETURNS:
        !   Total coagulation kernel for the two particles.

        Use SWPPART
        Use SWPPARAMS, only: COAG_TYPE
        Implicit None

        ! ARGUMENTS.
        Real :: k ! Return value.
        Type(StochParticle), Intent(IN) :: sp1, sp2 ! Particles for which to calculate the kernel.
        Real, Intent(IN)                :: T, P ! Temperature (K) and pressure (bar).
        Integer, Intent(IN)             :: maj  ! Flag for type of majorant kernel.

        ! VARIABLES.
        Real :: fm, sf

        ! EXECUTABLE CODE.
        If (COAG_TYPE==1) then
            k = FreeMolKernel(sp1, sp2, T, P, .False.)
        Else if (COAG_TYPE==2) then
            k = SlipFlowKernel(sp1, sp2, T, P, .False.)
        Else
            Select Case (maj)
                Case (NOMAJ)
                    ! Return true kernel (transition).
                    fm = FreeMolKernel(sp1, sp2, T, P, .False.)
                    sf = SlipFlowKernel(sp1, sp2, T, P, .False.)
                    k  = HalfHarmonicMean(fm, sf)
                Case (SLIPMAJ)
                    k  = SlipFlowKernel(sp1, sp2, T, P, .True.)
                Case (FREEMAJ)
                    k  = FreeMolKernel(sp1, sp2, T, P, .True.)
            End Select
        End if
    End Function

    ! -------------------------------------------------------

    Function FreeMolKernel(sp1, sp2, T, P, maj) Result(k)
        ! DESCRIPTION:
        !   Calculates the free-molecular coagulation
        !   kernel between two stochastic particles for
        !   the given temperature and pressure.
        !
        ! RETURNS:
        !   Free-molecular coagulation kernel for the
        !   two particles.

        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Real :: k ! Return value.
        Type(StochParticle), Intent(IN) :: sp1, sp2 ! Particles for which to calculate the kernel.
        Real, Intent(IN)                :: T, P ! Temperature (K) and pressure (bar).
        Logical, Intent(IN)             :: maj  ! Flag for majorant kernel.  Set to true to
                                                ! get majorant kernel, false to get real
                                                ! kernel.

        ! EXECUTABLE CODE.
        If (maj) Then
            ! Majorant kernel
            k = agg_eff * E_fm * C_fm * Sqrt(T) & 
                     * (sp1%Properties(iM_1_2) + sp2%Properties(iM_1_2)) &
                     * (sp1%Properties(iD2) + sp2%properties(iD2)) &
                     * C_fmmaj
        Else
            ! True kernel
            k = agg_eff * E_fm * C_fm * Sqrt(T) &
                     * Sqrt((1.0E0/sp1%Properties(iM)) + (1.0E0/sp2%Properties(iM))) &
                     * (sp1%Properties(iD) + sp2%properties(iD))**2
        End If
    End Function

    ! -------------------------------------------------------

    Function SlipFlowKernel(sp1, sp2, T, P, maj) Result(k)
        ! DESCRIPTION:
        !   Calculates the slip-flow coagulation
        !   kernel between two stochastic particles for
        !   the given temperature and pressure.
        ! RETURNS:
        !   Slip-flow coagulation kernel for the
        !   two particles.

        Use SWPCHEM_SHARED, only: Visc, Kn1
        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Real :: k
        Type(StochParticle), Intent(IN) :: sp1, sp2 ! Particles for which to calculate the kernel.
        Real, Intent(IN)                :: T, P ! Temperature (K) and pressure (bar).
        Logical, Intent(IN)             :: maj  ! Flag for majorant kernel.  Set to true to
                                                ! get majorant kernel, false to get real
                                                ! kernel.

        ! EXECUTABLE CODE.
        ! Majorant kernel is the same as the true kernel.
        k = 1.257E0 * Kn1(T,P) * (sp1%Properties(iD_2) + sp2%Properties(iD_2))
        k = k + (sp1%Properties(iD_1) + sp2%Properties(iD_1))
        k = agg_eff * C_sf * k * T * (sp1%Properties(iD) + sp2%properties(iD)) / Visc(T)
    End Function

    ! -------------------------------------------------------

    Function SlipFlowKernel_Rogak(sp1, sp2, T, P, maj) Result(k)
        ! DESCRIPTION:
        !   Calculates the slip-flow coagulation
        !   kernel between two stochastic particles for
        !   the given temperature and pressure.
        ! RETURNS:
        !   Slip-flow coagulation kernel for the
        !   two particles.

        Use SWPCHEM_SHARED, only: Visc, Kn1
        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Real :: k, A = 1.257E0, B = 0.400E0, C = 1.110E0
        Type(StochParticle), Intent(IN) :: sp1, sp2 ! Particles for which to calculate the kernel.
        Real, Intent(IN)                :: T, P ! Temperature (K) and pressure (bar).
        Logical, Intent(IN)             :: maj  ! Flag for majorant kernel.  Set to true to
                                                ! get majorant kernel, false to get real
                                                ! kernel.

        ! EXECUTABLE CODE.
        ! Majorant kernel is the same as the true kernel.
        k = B * Kn1(T,P) * ( sp1%Properties(iD_2) * exp(-C*sp1%Properties(iD)/Kn1(T,P) ) + &
                             sp2%Properties(iD_2) * exp(-C*sp2%Properties(iD)/Kn1(T,P) ) )
        k = k + A * Kn1(T,P) * (sp1%Properties(iD_2) + sp2%Properties(iD_2))
        k = k + (sp1%Properties(iD_1) + sp2%Properties(iD_1))
        k = agg_eff * C_sf * k * T * (sp1%Properties(iD) + sp2%properties(iD)) / Visc(T)
    End Function

    ! -------------------------------------------------------

    Real Function HalfHarmonicMean(x,y)
        ! DESCRIPTION:
        !	Calculates the harmonic mean of two
        !	numbers.
        !
        ! RETURNS:
        !	The harmonic mean.
        Implicit None
        Real, Intent(IN) :: x, y
        HalfHarmonicMean = (x * y) / (x + y)
    End Function

    ! -------------------------------------------------------

    Function RogakTransitonalFlowKernel(sp1, sp2, T, P, maj) Result(k)
        ! DESCRIPTION:
        !   Calculates the Rogak transitional coagulation
        !   kernel between two stochastic particles for
        !   the given temperature and pressure.
        ! RETURNS:
        !   Rogak Transitional coagulation kernel for the
        !   two particles.

        Use SWPCHEM_SHARED, only: Visc, Kn1
        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Real :: tmp, A = 1.257E0, B = 0.400E0, C = 1.110E0
        Real :: k ! Return value.
        Type(StochParticle), Intent(IN) :: sp1, sp2 ! Particles for which to calculate the kernel.
        Real, Intent(IN)                :: T, P ! Temperature (K) and pressure (bar).
        Logical, Intent(IN)             :: maj  ! Flag for majorant kernel.  Set to true to
                                                ! get majorant kernel, false to get real
                                                ! kernel.

        ! EXECUTABLE CODE.
        tmp = B * Kn1(T,P) * ( sp1%Properties(iD_2) * exp(-C*sp1%Properties(iD)/Kn1(T,P) ) + &
                               sp2%Properties(iD_2) * exp(-C*sp2%Properties(iD)/Kn1(T,P) ) )
        tmp = tmp + A * Kn1(T,P) * (sp1%Properties(iD_2) + sp2%Properties(iD_2))
        tmp = tmp + (sp1%Properties(iD_1) + sp2%Properties(iD_1))
        k = 2.0E0 * E_fm * C_fm * sqrt(T) * sqrt((1.0E0/sp1%Properties(iM)) + (1.0E0/sp2%Properties(iM)))
        k = C_sf * tmp * T / (sp1%Properties(iD) + sp2%Properties(iD)) / k / Visc(T)
        k = (1.0E0 + k) / (1.0E0 + 2.0E0*k * (1.0E0 + k) ) 
        k = agg_eff * C_sf * tmp * T * (sp1%Properties(iD) + sp2%Properties(iD)) * k / Visc(T)
    End Function

    ! -------------------------------------------------------
    ! NUCLEATION KERNEL ROUTINES.
    !
    !	These routines can be used to give information about
    !	the nucleation kernels for species.
    !
    ! -------------------------------------------------------

   Real Function NucFreeMolK(d1, d2, m1, m2)
        ! DESCRIPTION:
        !   Calculates the constant in a free-molecular
        !   kernel for two particles of different sizes.  To
        !   get the actual kernel value, this function needs to
        !   be multiplied by Sqrt(Temp).
        ! RETURNS:
        !   Free-molecular kernel without temperature
        !   dependency.        
        Implicit None
        Real, Intent(IN) :: d1, d2 ! Particle diameters (cm).
        Real, Intent(IN) :: m1, m2 ! Particle masses (g).
        NucFreeMolK = E_nuc*C_fm * Sqrt((1.0/m1) + (1.0/m2)) * ((d1 + d2)**2) 
    End Function
    
    ! -------------------------------------------------------

    Subroutine NucSlipFlowK(d1, d2, k1, k2)
        ! DESCRIPTION:
        !   Calculates the particle constant in the
        !   slip-flow kernel for two different particles.
        !   This function is used to initialise inception
        !   reaction rates.
        !   In order to get the true kernel value k1 must
        !   be multiplied by T/mu and k2 must be multiplied by
        !   T^2/muP where mu=Viscosity.  Then the two values
        !   must be added together.  T required in K and pressure
        !   in bar.
        ! RETURNS:
        !   Particle-dependent constants in slip-flow kernel.

        Use SWPPARAMS, only: KNUDSEN_K
        Implicit None

        ! ARGUMENTS.
        Real, Intent(IN)  :: d1, d2 ! Particle diameters (cm).
        Real, Intent(OUT) :: k1, k2 ! Output constants.

        ! VARIABLES.
        Real :: invd1, invd2

        ! EXECUTABLE CODE.
        invd1 = 1.0E0 / d1
        invd2 = 1.0E0 / d2
        k1 = C_sf * (d1 + d2)
        k2 = 1.257E0 * KNUDSEN_K * k1 * ((invd1*invd1) + (invd2*invd2))
        k1 = k1 * (invd1 + invd2)
    End Subroutine

    ! -------------------------------------------------------
    ! DIMERIZATION KERNEL ROUTINES.
    !
    !   These routines can be used to give information about
    !   the dimerization kernels for species.
    !
    ! -------------------------------------------------------

    Real Function DimerFreeMolK(d, m)
        ! DESCRIPTION:
        !   Calculates the constant in a free-molecular
        !   kernel for two particles of the same size.  To
        !   get the actual kernel value, this function needs to
        !   be multiplied by Sqrt(Temp).
        !
        ! RETURNS:
        !   Free-molecular kernel without temperature
        !   dependency.
        Use SWPPARAMS
        Implicit None
        Real, Intent(IN) :: d ! Particle diameter (cm).
        Real, Intent(IN) :: m ! Particle mass (g).
        DimerFreeMolK = E_nuc*C_fm * 4.0E0 * ROOT_TWO * d * d / Sqrt(m)
    End Function

    ! -------------------------------------------------------

    Subroutine DimerSlipFlowK(d, k1, k2)
        ! DESCRIPTION:
        !   Calculates the constants in the slip-flow
        !   kernel for two identical particles.  This
        !   function is used to initialise inception
        !   reaction rates.
        !   In order to get the true kernel value k1 must
        !   be multiplied by T/mu and k2 must be multiplied by
        !   T^2/muP where mu=Viscosity.  Then the two values
        !   must be added together.  T required in K and pressure
        !   in bar.
        ! RETURNS:
        !   Constants in slip-flow kernel.
        Use SWPPARAMS, only: KNUDSEN_K
        Implicit None
        Real, Intent(IN)  :: d      ! Particle diameter (cm).
        Real, Intent(OUT) :: k1, k2 ! Output constants.
        k1 = 4.0E0 * C_sf
        k2 = 1.257E0 * k1 * KNUDSEN_K / (d * d)
    End Subroutine

    ! -------------------------------------------------------
    ! CONDENSATION KERNEL ROUTINES.
    !
    !   These routines can be used to give information about
    !   the condensation kernels for species.
    !
    ! -------------------------------------------------------

    Subroutine FreeMolCondTerms(d, m, k1, k2, k3)
        ! DESCRIPTION:
        !   Calculates the free-molecular condensation terms
        !   for a species/particle of size x onto any soot
        !   particle.  There are three terms in the kernel.
        !
        !   Note: these terms are only valid if x is smaller
        !   than the subsequent particle onto which it
        !   condenses.
        !
        ! RETURNS:
        !   Free-molecular condensation terms via k1-k3.

        Implicit None

        ! ARGUMENTS.
	Real, Intent(IN)  :: d          ! Species diameter (cm).
        Real, Intent(IN)  :: m          ! Species mass (g).
	Real, Intent(OUT) :: k1, k2, k3 ! Condensations terms.

        ! EXECUTABLE CODE.
        k3 = cond_eff*E_cond * C_fm / Sqrt(m)
        k2 = d * k3 * 2.0 ! * 1.13 (A4 shape corrections).
        k1 = d * k2 / 2.0 ! * 1.28 = 1.13 * 1.13
    End Subroutine

End Module

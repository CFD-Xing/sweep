! *****************************************************************************
!
! File:				sharedchem.f90
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
!   Shared chemistry functions that are common to all methods of controlling
!   the gas-phase chemistry in Sweep.
!
! Functions:
!   Init		-	Initialises the shared chemistry components of
!					Sweep.
!	---------------------------------------------------------------------------
!   Visc		-	Calculates the viscosity of air for given temperature.
!   Kn			-	Calculates Knudsen number of a sphere with the 
!					properties of one carbon atom in air.
!   RadicalSiteFraction	-	Calculates the fraction of soot surface active sites
!					which are radical sites (i.e. missing a hydrogen
!					atom).
!   Alpha		-	Calculates the active site fraction given the
!					temperature and first size moment using the correlation
!					given by Frenklach et al. (2000).
! *****************************************************************************


Module SWPCHEM_SHARED
    ! -------------------------------------------------------
    ! IMPORT PUBLIC MODULE FUNCTIONS.
    ! -------------------------------------------------------
    Use SWPERR ! Error codes.
    ! -------------------------------------------------------

    Implicit None
    Public

    Type SharedChemData
        Integer :: NCHEM, NCHEMEXT        ! Counts.
        Integer :: iT, iP, iAlpha, iRad,iRho, ixx   ! Useful indices in chemistry arrays.
        Integer :: iH2, iH, iA4, iC2H2, & ! Species required for radical
                   iOH, iO2, iCO, iH2O    ! site fraction calculation.
    End Type

    Contains

    ! -------------------------------------------------------

    Subroutine Init(chem, speciesNames, mech, flag)
        ! DESCRIPTION:
        !   Initialises shared chemistry components.

        Use StrConv
        Use SWPPARAMS
        Use SWPMECH_TYPES
        Implicit None

        ! ARGUMENTS.
        Type(SharedChemData), Intent(OUT) :: chem      ! Chemistry to initialise.
        Character*(*), Intent(IN)   :: speciesNames(:) ! Species in the gas-phase chemistry.	
        Type(Mechanism), Intent(IN) :: mech            ! Mech to use for initialisation.
        Integer, Intent(OUT)        :: flag            ! Error flag.

        ! VARIABLES.
        Integer :: nsp

        ! EXECUTABLE CODE.
        flag = 0

        If (mech%UseHACA) Then
           ! Get indices of species required for HACA.
           chem%iH2   = IndexOf("H2", speciesNames)
           chem%iH    = IndexOf("H", speciesNames)
           chem%iA4   = IndexOf("A4", speciesNames)
           chem%iC2H2 = IndexOf("C2H2", speciesNames)
           chem%iOH   = IndexOf("OH", speciesNames)
           chem%iO2   = IndexOf("O2", speciesNames)
           chem%iCO   = IndexOf("CO", speciesNames)
           chem%iH2O  = IndexOf("H2O", speciesNames)
        Else
           chem%iH2   = 0
           chem%iH    = 0
           chem%iA4   = 0
           chem%iC2H2 = 0
           chem%iOH   = 0
           chem%iO2   = 0
           chem%iCO   = 0
           chem%iH2O  = 0
        End If

        ! Calculate the memory requirement for the
        ! chemistry (species, temperature and pressure).
        nsp           = Size(speciesNames)
        chem%NCHEM    = nsp + 2
        chem%NCHEMEXT = chem%NCHEM
        chem%iT       = nsp + 1
        chem%iP       = nsp + 2
        chem%iRho     = IndexOf("rho", speciesNames)
        chem%ixx      = IndexOf("xx", speciesNames)

        If (mech%UseHACA) Then
            If (mech%ActSurfModel == ACT_SURF_PROFILE) Then
                chem%NCHEM    = chem%NCHEM + 2
                chem%NCHEMEXT = chem%NCHEMEXT + 1
                chem%iAlpha   = nsp + 3
                chem%iRad     = nsp + 4
            Else
                chem%NCHEM  = chem%NCHEM + 1
                chem%iAlpha = 0
                chem%iRad   = nsp + 3
            End If
        Else
            chem%iAlpha = 0
            chem%iRad   = 0
        End If
    End Subroutine

    ! -------------------------------------------------------
    ! CHEMICAL PROPERTY ROUTINES.
    !
    !   These routines return the values of useful chemical
    !   properties.
    !
    ! -------------------------------------------------------

    Real Function Visc(T)
        ! DESCRIPTION:
        !   Calculate the viscosity or air at the given
        !   temperature (K).
        ! RETURNS:
        !   Viscosity of air.
        Implicit None
        Real, Intent(IN) :: T ! Temperature (K).
        Visc = 14.58E-6 * (T**1.5E0) / (T + 110.4E0)
    End Function

    ! -------------------------------------------------------

    Real Function Kn(d, T, P)
        ! DESCRIPTION:
        !   Calculate the Knudsen number of a sphere with
        !   the given diameter at the given temperature (K)
        !   and pressure (bar).
        ! RETURNS:
        !   Knudsen number.
        Use SWPPARAMS, only: KNUDSEN_K
        Implicit None
        Real, Intent(IN) :: d    ! Particle diameter (cm).
        Real, Intent(IN) :: T, P ! Temperature (K) & pressure (bar).
        Kn = (KNUDSEN_K * T) / (P * d)
    End Function

    ! -------------------------------------------------------

    Real Function Kn1(T, P)
        ! DESCRIPTION:
        !   Calculate the Knudsen number of a sphere with
        !   diameter 1 unit at the given temperature (K)
        !   and pressure (bar).
        ! RETURNS:
        !   Knudsen number of unit diameter (1 cm) sphere.
        Use SWPPARAMS, only: KNUDSEN_K
        Implicit None
        Real, Intent(IN) :: T, P ! Temperature (K) & pressure (bar).
        Kn1 = (KNUDSEN_K * T) / P
    End Function

    ! -------------------------------------------------------

    Real Function RadicalSiteFraction(chem, chemdata)
        ! DESCRIPTION:
        !   Returns the fraction of active sites (given by
        !   alpha) which are radical sites (i.e. missing
        !   a hydrogen atom).
        !
        !   This is given by the HACA model and is a
        !   function of the gas-phase chemistry only.  The
        !   fraction of active sites is a function of the
        !   particle properties, in particular age.
        ! RETURNS:
        !   Steady-state radical site fraction for given
        !   chemical conditions.

        Use SWPPARAMS
        Implicit None

        ! ARGUMENTS.
	Real, Intent(IN) :: chem(:) ! Gas-phase species concentrations.
        Type(SharedChemData), Intent(IN) :: chemdata ! Chemistry data to use.

        ! VARIABLES.
        Double Precision :: r1f, r1b, r2f, r2b, r3f, & ! Reaction rates.
                            r4f, r5f, rdenom, RT

        ! EXECUTABLE CODE.
        ! Reaction rates for the six ABF reactions.
        RT  = Dble(RCAL * chem(chemdata%iT))
        r1f = 4.2D+13 * Exp(-13.0/RT)                                * chem(chemdata%iH)
        r1b = 3.9D+12 * Exp(-11.0/RT)                                * chem(chemdata%iH2)
        r2f = 1.0D+10 * Exp(-1.43/RT) * (chem(chemdata%iT) ** 0.734) * chem(chemdata%iOH)
        r2b = 3.68D+8 * Exp(-17.1/RT) * (chem(chemdata%iT) ** 1.139) * chem(chemdata%iH2O)
        r3f = 2.0D+13                                                * chem(chemdata%iH)
        r4f = 8.0D+07 * Exp( -3.8/RT) * (chem(chemdata%iT) ** 1.56 ) * chem(chemdata%iC2H2)
        r5f = 2.2E+12 * Exp( -7.5/RT)                                * chem(chemdata%iO2)
        rdenom = (r1b+r2b+r3f+r4f+r5f)

        ! Return fraction of active radical sites by using a
        ! steady-state assumption on the above reactions.		
        If (rdenom > 0.0D0) Then
           RadicalSiteFraction = Real((r1f+r2f)/rdenom)
        Else
           RadicalSiteFraction = 0.0E0
        End If
    End Function

    ! -------------------------------------------------------

    Real Function ABFAlpha(T, M1)
        ! DESCRIPTION:
        !   Calculates the active site fraction on a soot
        !   particle surface using the correlation for
        !   temperature given by Frenklach.
        ! RETURNS:
        !   Active site fraction (Alpha).

        Use SWPPARAMS
        Implicit None

        ! ARGUMENTS.
        Real, Intent(IN) :: T, M1 ! Temperature and first size moment
                                  ! for which to correlate alpha.

        ! VARIABLES.
        Real :: a, b ! Intermediate values.

        ! EXECUTABLE CODE.
        If (M1 > 0.0E0) Then
            a = 12.65E0 - 5.63E-03 * T
            b = -1.38E0 + 6.80E-04 * T
            ABFAlpha = Tanh((a / Log10(M1)) + b)
            ABFAlpha = Max(0.0D0, ABFAlpha)
        Else
            ABFAlpha = 0.0E0
        End If
    End Function

End Module

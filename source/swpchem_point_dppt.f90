! *****************************************************************************
!
! File:					swpchem_point_dppt.f90
! Project:				Sweep 2
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
!	Gas-phase chemistry controller for Sweep.  Handles storage and manipulation
!	of chemical species data.  This module does NOT define the surface
!	chemistry of soot particles, that is controlled by the SWPMECH module.
!
!	This flavour of the chemistry module is for a single time point
!	solution for use when the soot is acting on the chemistry and the
!	whole chemistry profile is not available.
!
!   Additionally this flavour accepts chemistry in DOUBLE PRECISION and 
!   takes a pointer to the chemistry array, rather than copying the
!   data each time it is set.
!
!	Note this module requires chemistry as MOLAR CONCENTRATIONS.
!
! Functions:
!	---(Public interface)------------------------------------------------------
!	LoadChem			-	Loads gas-phase chemistry into Sweep for a new
!							simulation.
!	UnloadChem			-	Gets the gas-phase chemistry from Sweep.
!	---(Sweep only routines)---------------------------------------------------
!	InitChem			-	Initialises the chemistry module for Sweep, enabling
!							it to be loaded with chemistry data.
!	ChangeChem			-	Changes the value of a chemistry variable by a
!							given amount.
!	---(Get chem info)---------------------------------------------------------
!	GetChem				-	Gets the gas-phase chemical conditions at a given
!							time.
!	GetTemperature		-	Returns the gas temperature at the given time.
!	---(Auxiliary routines)----------------------------------------------------
!	DataSpacing			-	Returns the time interval between chemistry data
!							points which enclose the given time.
! *****************************************************************************


Module SWPCHEM
	! -------------------------------------------------------
	! IMPORT PUBLIC MODULE FUNCTIONS.
	! -------------------------------------------------------
    Use SWPERR
	Use SWPCHEM_SHARED, InitSharedChem => Init ! Shared chemistry components.
    Use SWPPARAMS, only: FIXED_CHEM, VARIABLE_CHEM
	! -------------------------------------------------------

	Implicit None
	Public

	! -------------------------------------------------------
	! REQUIRED CHEMISTRY PARAMETERS.
	! -------------------------------------------------------

	! Type of chemistry data stored in this module.
	Integer, Parameter	::	CHEM_TYPE = VARIABLE_CHEM

	! -------------------------------------------------------
	! CHEMISTRY DATA.
	! -------------------------------------------------------

    Type ChemData
        Type(SharedChemData) :: Shared ! Data used by all Sweep chem' modules.
        Double Precision, Pointer :: Chem(:)     ! Chemistry data.
        Double Precision, Pointer :: InitChem(:) ! Initial chemistry data (if variable).
        Double Precision          :: InitTime    ! Time applicable for initial chemistry data.
    End Type

	Contains

	! -------------------------------------------------------
	! SWEEP INTERFACE ROUTINES.
	!
	!	These functions are called by external code and not
	!	used by Sweep.
	!
	! -------------------------------------------------------

	Subroutine LoadChem(cdata, t, chem)
		! DESCRIPTION:
		!	Loads the gas-phase chemistry into the data
        !   structures.

        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: cdata ! Where to put chemistry.
		! Flow time at which chemistry values apply and the
		! gas-phase chemistry (incl. Temperature and pressure)
		! to load into Sweep.
		Double Precision, Intent(IN) :: t
        Double Precision, Pointer    :: chem(:)

		! EXECUTABLE CODE.
		cdata%InitTime = t
		cdata%InitChem(1:cdata%Shared%NCHEMEXT) = chem

        If (cdata%Shared%iRad > 0) Then
		    cdata%InitChem(cdata%Shared%iRad) = Dble(RadicalSiteFraction(chem, cdata%Shared))
        End If

		cdata%Chem => chem
	End Subroutine

	! -------------------------------------------------------

	Subroutine UnloadChem(cdata, t, chem)
		! DESCRIPTION:
		!	Gets the gas-phase chemistry from Sweep.

        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: cdata ! Where to get chemistry from.
		! Gas-phase chemistry (incl. Temperature and pressure)
		! to get from Sweep.
		Double Precision, Intent(OUT) :: t, chem(:)

		! VARIABLES.

		! EXECUTABLE CODE.
		t = 0.0E0

        If (Size(chem) >= cdata%Shared%NCHEMEXT) Then
            If (CHEM_TYPE == VARIABLE_CHEM) Then
		        chem(1:cdata%Shared%NCHEMEXT) = cdata%Chem(1:cdata%Shared%NCHEMEXT)
            Else
                chem(1:cdata%Shared%NCHEMEXT) = cdata%InitChem(1:cdata%Shared%NCHEMEXT)
            End If
        Else
            chem = 0.0E0
        End If
	End Subroutine

	! -------------------------------------------------------
	! SWEEP INTERNAL ROUTINES.
	!
	!	These functions are for internal use by Sweep only,
	!	and not to be called by other code.
	!
	! -------------------------------------------------------

	Subroutine InitChem(chem, speciesNames, mech, flag)
		! DESCRIPTION:
		!	Initialises the chemistry components of Sweep.
        !   Basic chemistry components are initialised using
        !   InitShared().  Components only applicable to this
        !   chemistry module are initialised here.

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(OUT) :: chem            ! Chemistry to initialise.
		Character*(*), Intent(IN)   :: speciesNames(:) ! Species in gas-phase chemistry.
        Type(Mechanism), Intent(IN) :: mech            ! Mechanism to use.
		Integer, Intent(OUT)	    :: flag            ! Error flag.

		! EXECUTABLE CODE.
		flag = 0

		! Initialise shared chemistry components.
		Call InitSharedChem(chem%Shared, speciesNames, mech, flag)

		! Clear memory currently used for the chemistry.
		Deallocate(chem%InitChem, STAT=flag)

		! Allocate memory for new chemistry.
		Allocate(chem%InitChem(chem%Shared%NCHEM), STAT=flag)

        If (flag /= 0) Then
            flag = CHEM_INIT_ERR
            Return
        End If

		! Clear memory.
		Nullify(chem%Chem)
		chem%InitChem = 0.0E0
		chem%InitTime = 0.0E0
	End Subroutine

	! -------------------------------------------------------

	Subroutine DeleteChem(chem)
		! DESCRIPTION:
		!	Free memory associated with a chemistry data type.
		Implicit None
        Type(ChemData), Intent(INOUT) :: chem ! Chemistry to free.
        Integer :: err = 0
		Deallocate(chem%InitChem, chem%Chem, STAT=err)
	End Subroutine

	! -------------------------------------------------------
	
	Subroutine ChangeChem(chem, t, id, change)
		! DESCRIPTION:
		!	Changes the value of a chemical species
		!	concentration or temperature or pressure.

		Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: chem ! Chemistry to change.
		Real, Intent(IN)	::	t      ! Time at which to apply the change.
		Integer, Intent(IN)	::	id     ! Index/ID of the variable to change.
		Real, Intent(IN)	::	change ! Change amount.

		! EXECUTABLE CODE.
        If (CHEM_TYPE == VARIABLE_CHEM) Then
		    ! Change chemical variable, but don't allow it
		    ! to become negative.
!            If (MSGS .And. (chem%Chem(id)+Dble(change)<=0.0D0)) Then
!                Print *, "Sweep:\> Chemistry reduced to/below zero!", id
!            End If
		    chem%Chem(id) = Max(0.0D0, chem%Chem(id) + Dble(change))
        End If
	End Subroutine

	! -------------------------------------------------------
	
	Subroutine SetChem(chem, t, id, value)
		! DESCRIPTION:
		!	Set the value of a chemical species
		!	concentration or temperature or pressure.

		Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: chem ! Chemistry to change.
		Real, Intent(IN)	::	t      ! Time at which to apply the change.
		Integer, Intent(IN)	::	id     ! Index/ID of the variable to change.
		Real, Intent(IN)	::	value  ! New value.

		! EXECUTABLE CODE.
        If (CHEM_TYPE == VARIABLE_CHEM) Then
		    ! Change chemical variable, but don't allow it
		    ! to become negative.
		    chem%Chem(id) = Max(0.0D0, Dble(value))
        End If
	End Subroutine

	! -------------------------------------------------------

	Function GetChem(chem, t, actualTime, flag)
		! DESCRIPTION:
		!	Gets the chemical conditions at the given time.  If
        !   the time is in the past, then it returns a linear
        !   interpolation between the initial conditions and the
        !   current conditions.
		! RETURNS:
		!	Chemical conditions at time t.

		Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: chem       ! Chemistry from which to get.
		Real, Intent(IN)           :: t          ! Time for which to get chemical conditions.
        Real, Intent(IN)           :: actualTime ! Actual simulation time.
		Integer, Intent(OUT)       :: flag       ! Error flag.
		Real ::	GetChem(chem%Shared%NCHEM)       ! Return value.


		! EXECUTABLE CODE.
		flag = 0

        If (CHEM_TYPE == VARIABLE_CHEM) Then
		    If (t == actualTime) Then
			    ! Just return the current conditions.
			    GetChem(1:chem%Shared%NCHEMEXT) = Real(chem%Chem)
		    Else
			    ! Return a linear interpolation between the initial
			    ! conditions and current conditions.
                If (actualTime == chem%InitTime) Then
                    GetChem = Real(chem%InitChem)
                Else
			        GetChem(1:chem%Shared%NCHEMEXT) = Real(chem%InitChem(1:chem%Shared%NCHEMEXT) + &
                                                           ((chem%Chem(1:chem%Shared%NCHEMEXT) - &
                                                            chem%InitChem(1:chem%Shared%NCHEMEXT)) * &
                                                           (Dble(t) - chem%InitTime) / &
                                                           (Dble(actualTime) - chem%InitTime)))
                End If
		    End If

            If (chem%Shared%iRad > 0) Then
		        ! Calculate the radical site fraction for this chemistry.
		        GetChem(chem%Shared%iRad) = RadicalSiteFraction(GetChem, chem%Shared)
            End If
        Else
            GetChem = Real(chem%InitChem)
        End If
	End Function

	! -------------------------------------------------------

	Real Function GetTemperature(chem, t, actualTime)
		! DESCRIPTION:
		!	Gets the temperature at the given time.
		! RETURNS:
		!	Temperature at time t.

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: chem       ! Chemistry from which to get.
		Real, Intent(IN)           :: t          ! Time for which to get chemical conditions.
        Real, Intent(IN)           :: actualTime ! Actual simulation time.

		! EXECUTABLE CODE.
        If (CHEM_TYPE == VARIABLE_CHEM) Then
		    If (t == actualTime) Then
			    ! Just return the current temperature.
			    GetTemperature = Real(chem%Chem(chem%Shared%iT))
		    Else
			    ! Return a linear interpolation between the initial
			    ! temperature and current temperature.
			    GetTemperature = Real(chem%InitChem(chem%Shared%iT) + ((chem%Chem(chem%Shared%iT) - chem%InitChem(chem%Shared%iT)) * (Dble(t) - chem%InitTime) / &
												                       (Dble(actualTime) - chem%InitTime)))
		    End If
        Else
            GetTemperature = Real(chem%InitChem(chem%Shared%iT))
        End If
	End Function

	! -------------------------------------------------------

	Real Function DataSpacing(chem, t)
		! DESCRIPTION:
		!	Finds the length of the time interval
		!	containing time.  That is the interval between
		!	the data points which bound t.  For dynamic chemisty
        !   this is the difference between the initial time and the
        !   current time.
		! RETURNS:
		!	Time interval (s).
		Implicit None
        Type(ChemData), Intent(IN) :: chem
		Real, Intent(IN) :: t ! Time (s).
		DataSpacing = 1.0E+30
	End Function

	! -------------------------------------------------------

	Subroutine WriteChemistry(chem, fnum)
		! DESCRIPTION:
		!	Writes chemistry to a binary file.
		Implicit None
        Type(ChemData), Intent(IN) :: chem ! Chemistry to output.
        Integer, Intent(IN)        :: fnum ! UNIT number on which to output.
        Write(fnum) chem%Shared, chem%Chem, chem%InitChem, chem%InitTime
	End Subroutine

	! -------------------------------------------------------

	Subroutine ReadChemistry(chem, fnum)
		! DESCRIPTION:
		!	Reads chemistry from a binary file.
		Implicit None

        ! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: chem ! Chemistry to read.
        Integer, Intent(IN)           :: fnum ! UNIT number on which to output.

        ! EXECUTABLE CODE.
        Read(fnum) chem%Shared, chem%Chem, chem%InitChem, chem%InitTime
	End Subroutine

End Module

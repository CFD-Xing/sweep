! *****************************************************************************
!
! File:					swpchem_profile.f90
! Project:				Sweep 2
! Author(s):			Matthew Celnik (msc37) & Rob Patterson (riap2)
!
! Copyright (C) 2006  Matthew S Celnik
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
!	chemistry of soot particles.
!
!	This flavour of the chemistry module is for an array of time points, i.e.
!	for when the entire chemistry profile is known and is not altered by Sweep.
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
!	---(Get chem info)---------------------------------------------------------
!	GetChem				-	Gets the gas-phase chemical conditions at a given
!							time.
!	GetTemperature		-	Returns the gas temperature at the given time.
!	---(Auxiliary routines)----------------------------------------------------
!	FindTime			-	Locates the index of the last chemistry data point
!							before the given time.
!	DataSpacing			-	Returns the time interval between chemistry data
!							points which enclose the given time.
! *****************************************************************************

Module SWPCHEM
	! -------------------------------------------------------
	! IMPORT PUBLIC MODULE FUNCTIONS.
	! -------------------------------------------------------
    Use SWPERR ! Error codes.
	Use SWPCHEM_SHARED, InitShared => Init ! Shared chemistry components.
    Use SWPPARAMS, only: FIXED_CHEM, VARIABLE_CHEM
	! -------------------------------------------------------

	Implicit None
	Public

	! -------------------------------------------------------
	! PROFILE CHEMISTRY PARAMETERS.
	! -------------------------------------------------------

	! Type of chemistry data stored in this module.  The
	! enumeration is stored in the SWPSharedChemistry module.
	Integer, Parameter	::	CHEM_TYPE = FIXED_CHEM

	! -------------------------------------------------------
	! CHEMISTRY DATA.
	! -------------------------------------------------------

    Type ChemData
        Type(SharedChemData) :: Shared    ! Data used by all Sweep chem' modules.
        Integer              :: NPOINTS=0 ! Counts.
        Real, Pointer        :: Times(:)  ! Time profile.
        Real, Pointer        :: Chem(:,:) ! Chemistry profile.
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
		!	Loads gas-phase chemistry data into a ChemData
        !   structure.

        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: cdata ! Where to put chemistry.
		! Time for which the chemistry profile applies and the
		! gas-phase chemistry (incl. Temperature and pressure)
		! to load into Sweep.
		Real, Intent(IN) :: t(:), chem(:,:)

		! VARIABLES.
		Integer :: LT, LC, err

		! EXECUTABLE CODE.

		! Clear memory currently allocated for chemistry.
		Deallocate(cdata%Times, cdata%Chem, STAT=err)

		! Note the number of data points.
		LT = Size(t)
		LC = Size(chem, 2)
		cdata%NPOINTS = Min(LT, LC)

		! Allocate memory to store chemistry profile.
		Allocate(cdata%Times(cdata%NPOINTS), &
                 cdata%Chem(cdata%Shared%NCHEM,cdata%NPOINTS), STAT=err)

		! Copy the chemistry profile to Sweep.
		cdata%Times = t(1:cdata%NPOINTS)
		cdata%Chem(1:cdata%Shared%NCHEMEXT,1:cdata%NPOINTS) = chem(1:cdata%Shared%NCHEMEXT,1:cdata%NPOINTS)

        If (cdata%Shared%iRad > 0) Then
            ! Calculate HACA parameters.
		    Do LT = 1, cdata%NPOINTS
			    cdata%Chem(cdata%Shared%iRad,LT) = RadicalSiteFraction(cdata%Chem(:,LT), cdata%Shared)
		    End Do
        End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine UnloadChem(cdata, t, chem, flag)
		! DESCRIPTION:
		!	Gets the gas-phase chemistry from a data
        !   type.

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: cdata     ! Data type to get from.
		Real, Intent(OUT)	       :: t(:)      ! Time profile.
        Real, Intent(OUT)          :: chem(:,:) ! Chemistry profile.
		Integer, Intent(OUT)       :: flag      ! Error flag.

		! VARIABLES.
		Integer :: L

		! EXECUTABLE CODE.
		flag = 0

		If (Size(chem, 1) >= cdata%Shared%NCHEMEXT) Then
			! Output the time profile.
			L = Min(Size(t), Size(cdata%Times))
			t(1:L) = cdata%Times(1:L)

			! Output chemistry profile.
			L = Min(Size(chem), Size(cdata%Chem))
			chem(1:cdata%Shared%NCHEMEXT,1:L) = cdata%Chem(1:cdata%Shared%NCHEMEXT,1:L)
		Else
			! There is not enough space in the 
			! output array to store the chemistry.
			flag = UNLOAD_CHEM_ERR
		End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine CopyChem(c1, c2)
		! DESCRIPTION:
		!	Copies data from one chemistry set to another.

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(INOUT) :: c1, c2

        ! EXECUTABLE CODE.
        Call LoadChem(c2, c1%Times, c1%Chem)
        c2%Shared = c1%Shared
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
		!	Initialises the chemistry components.

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(OUT) :: chem            ! Chemistry to initialise.
		Character*(*), Intent(IN)	:: speciesNames(:) ! Species in gas-phase chemistry.
        Type(Mechanism), Intent(IN) :: mech            ! Mech to use for initialisation.
		Integer, Intent(OUT)		:: flag            ! Error flag.

		! EXECUTABLE CODE.
		flag = 0

		! Initialise chemistry components.
		Call InitShared(chem%Shared, speciesNames, mech, flag)
        chem%NPOINTS = 0
	End Subroutine

	! -------------------------------------------------------

	Subroutine DeleteChem(chem)
		! DESCRIPTION:
		!	Free memory associated with a chemistry data type.
		Implicit None
        Type(ChemData), Intent(INOUT) :: chem ! Chemistry to free.
        Integer :: err = 0
		Deallocate(chem%Times, chem%Chem, STAT=err)
	End Subroutine

	! -------------------------------------------------------

	Function GetChem(chem, t, actualTime, flag)
		! DESCRIPTION:
		!	Gets the chemical conditions at the given time.
        ! RETURNS:
        !   Chemical conditions at time t.

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: chem       ! Data from which to get chemistry.
		Real, Intent(IN)           :: t          ! Time for which to get chemical conditions.
        Real, Intent(IN)           :: actualTime ! Actual simulation time.
		Integer, Intent(OUT)       :: flag       ! Error flag.
		Real ::	GetChem(chem%Shared%NCHEM)       ! Return value.

		! VARIABLES.
		Integer :: i, j

		! EXECUTABLE CODE.
		flag = 0

		i = FindTime(t, chem%Times, chem%NPOINTS)
		
		If (i < chem%NPOINTS) Then
            ! Perform linear interpolation between the
            ! two bounding data points.
			j = i + 1
			GetChem = chem%Chem(:,i) + ((chem%Chem(:,j) - chem%Chem(:,i)) * (t - chem%Times(i)) / &
										(chem%Times(j) - chem%Times(i)))
		Else
            ! This is the last data point, and we don't do
            ! extrapolation.
			GetChem = chem%Chem(:,i)
		End If
	End Function

	! -------------------------------------------------------

	Real Function GetTemperature(chem, t, actualTime)
		! DESCRIPTION:
		!	Gets the temperature at the given time.
        ! RETURNS:
        !   Temperature at time t.

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: chem       ! Data from which to get chemistry.
		Real, Intent(IN)           :: t          ! Time for which to get chemical conditions.
        Real, Intent(IN)           :: actualTime ! Actual simulation time.

		! VARIABLES.
		Integer :: i, j

		! EXECUTABLE CODE.
		i = FindTime(t, chem%Times, chem%NPOINTS)
		
		If (i < chem%NPOINTS) Then
            ! Perform linear interpolation between the
            ! two bounding data points.
			j = i + 1
			GetTemperature = chem%Chem(chem%Shared%iT,i) + ((chem%Chem(chem%Shared%iT,j) - chem%Chem(chem%Shared%iT,i)) * (t - chem%Times(i)) / &
												            (chem%Times(j) - chem%Times(i)))
		Else
            ! This is the last data point, and we don't do
            ! extrapolation.
			GetTemperature = chem%Chem(chem%Shared%iT,i)
		End If
	End Function

	! -------------------------------------------------------

	Integer Function FindTime(t, times, N)
		! DESCRIPTION:
		!	Locates the index of the last chemistry data point
		!	before the given time.  This routine uses a binary
        !   search.
		! RETURNS:
		!	Index of chemistry data point.

		Implicit None

		! ARGUMENTS.
        Integer, Intent(IN) :: N
		Real, Intent(IN)    :: t, times(N) ! Time (s).

		! VARIABLES.
		Integer :: LB, UB

		! EXECUTABLE CODE.
		LB = 1
		UB = N

		Do While (UB - LB > 1)
			If (times((UB+LB)/2) <= t ) Then
				! Move to right half of interval.
				LB = (UB + LB) / 2
			Else
				! Move to left half of interval.
				UB = (UB + LB) / 2
			End If
		End Do

		If (times(UB) > t) Then
			! Time is within range of data.
			FindTime = LB
		Else
			FindTime = N
		End If
	End Function

	! -------------------------------------------------------

	Real Function DataSpacing(chem, t)
		! DESCRIPTION:
		!	Finds the length of the time interval
		!	containing time.  That is the interval between
		!	the data points between which lies t.
		! RETURNS:
		!	Time interval (s).

		Implicit None

		! ARGUMENTS.
        Type(ChemData), Intent(IN) :: chem ! Data from which to get chemistry.
		Real, Intent(IN)           :: t    ! Time (s).

		! VARIABLES.
		Integer :: i

		! EXECUTABLE CODE.
		! Find the data point before time.
		i = FindTime(t, chem%Times, chem%NPOINTS)

		If (i < chem%NPOINTS) Then
			! We have data at a later time.
			DataSpacing = chem%Times(i+1) - chem%Times(i)
		Else
			DataSpacing = 1.0D+30 ! Almost +inf
		End If
	End Function

	! -------------------------------------------------------

	Subroutine WriteChemistry(chem, fnum)
		! DESCRIPTION:
		!	Writes chemistry to a binary file.
		Implicit None
        Type(ChemData), Intent(IN) :: chem ! Chemistry to output.
        Integer, Intent(IN)        :: fnum ! UNIT number on which to output.
        Write(fnum) chem%Shared, chem%NPOINTS, chem%Chem, chem%Times
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
        Read(fnum) chem%Shared, chem%NPOINTS, chem%Chem, chem%Times
	End Subroutine

End Module

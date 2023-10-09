! *****************************************************************************
!
! File:				swpmech.f90
! Project:			Sweep 2
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
!   Routines which act on a Sweep mechanism type.
! *****************************************************************************

Module SWPMECH
    Use SWPERR
    Use SWPMECH_TYPES
    Implicit None
    Public

    Contains

        ! -------------------------------------------------------

        Subroutine InitMech(mech, mechFile, speciesNames, flag)
                ! DESCRIPTION:
                !	Initialises a soot mechanism.

                Use SWPMECH_READER
                Implicit None

                ! ARGUMENTS.
                Type(Mechanism), Intent(OUT) :: mech            ! Mechanism to initialise.
                Character(LEN=*), Intent(IN) :: mechFile        ! Mechanism file name.
                Character(LEN=*), Intent(IN) :: speciesNames(:) ! Array of chemical species names.
                Integer, Intent(OUT)         :: flag            ! Error flag.

                ! EXECUTABLE CODE.
                flag = 0

                ! Clear current mechanism memory.
                Call DeleteMech(mech)

                ! Read the mechanism from file into temporary pointer arrays,
                ! and get the counts from the file.
                Call ReadMechXML(mechFile, mech, speciesNames, flag)

                If (flag /= 0) Then
                    flag = LOAD_MECH_ERR
                    Return
                End If
        
                ! Complete groups.
                mech%Groups(mech%GroupCount-2) = NCOAGTERMS
                mech%Groups(mech%GroupCount-1) = 1
                mech%Groups(mech%GroupCount)   = 1
                mech%GroupNames(mech%GroupCount-2) = "Coagulation"
                mech%GroupNames(mech%GroupCount-1) = "Inflow"
                mech%GroupNames(mech%GroupCount)   = "Outflow"

                ! Indices into the sweep internal rate array.
                mech%ISR  = mech%InceptionCount + 1       ! Index of first surface reaction (or condensation).
                mech%ICG  = mech%ISR + mech%ReactionCount ! Index of first coagulation term.
                mech%IIN  = mech%ICG + NCOAGTERMS         ! Index of inflow term.
                mech%IOUT = mech%IIN + 1                  ! Index of outflow term.
        End Subroutine

        ! -------------------------------------------------------

        Subroutine DeleteMech(mech)
                ! DESCRIPTION:
                !	Frees memory allocated to a soot mechanism.
                Implicit None
                Type(Mechanism), Intent(OUT) :: mech ! Mechanism to initialise.
                Integer :: err = 0
                Deallocate(mech%Inceptions, mech%Reactions, mech%Groups, &
                           mech%GroupNames, mech%DeferMask, STAT=err)
        End Subroutine

        ! -------------------------------------------------------
        ! MECHANISM INFORMATION FUNCTIONS.
        !
        !	These routines provide information about the
        !	currently loaded soot mechanism.
        !
        ! -------------------------------------------------------

        Integer Function DeferredProcessCount(mech)
                ! DESCRIPTION:
                !	Returns the number of soot processes that
                !	are being deferred.
                ! RETURNS:
                !	Deferred process count.
                Implicit None
                Type(Mechanism), Intent(IN) :: mech
                DeferredProcessCount = Count(mech%DeferMask)
        End Function

        ! -------------------------------------------------------

        Integer Function NonDeferredProcessCount(mech)
                ! DESCRIPTION:
                !	Returns the number of soot processes that
                !	are not being deferred.
                ! RETURNS:
                !	Non-deferred process count.
                Implicit None
                Type(Mechanism), Intent(IN) :: mech
                NonDeferredProcessCount = Count(.Not. mech%DeferMask)
        End Function

        ! -------------------------------------------------------

        Subroutine SetSinteringParams(mech, a, b)
                ! DESCRIPTION:
                !	Sets the sintering parameters.
                Implicit None
                Type(Mechanism), Intent(INOUT) :: mech
                Real, Intent(IN) :: a, b
                mech%SintParams = (/a,b/)
        End Subroutine

        ! -------------------------------------------------------

        Subroutine GetRateTermNames(mech, names)
                ! DESCRIPTION:
                !	Returns the grouped process names in an array.
                Implicit None
                Type(Mechanism), Intent(IN)   :: mech
                Character(LEN=*), Intent(OUT) :: names(mech%GroupCount)
                names = mech%GroupNames
        End Subroutine
End Module

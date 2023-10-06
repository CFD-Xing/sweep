! *****************************************************************************
!
! File:				swperr.f90
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
!   Definition of error codes that can be reported by Sweep.  In all circumstances
!   an error code of 0 means that the code executed successfully.  Error codes
!   less than zero (<0) are proper errors, and codes greater than zero (>0) are
!   merely warnings.
!
! *****************************************************************************

Module SWPERR
    ! -------------------------------------------------------
    ! IMPORT REQUIRED MODULES.
    ! -------------------------------------------------------
    ! None required.
    ! -------------------------------------------------------

    Implicit None
    Public

    ! -------------------------------------------------------
    ! STANDARD ERROR CODES.
    ! -------------------------------------------------------

    Integer, Parameter :: MEMORY_ERR                = -101

    ! -------------------------------------------------------
    ! STOCHASTIC STEPPER ERROR CODES.
    ! -------------------------------------------------------

    Integer, Parameter :: STEP_END_FAIL_ERR         = -20001
    Integer, Parameter :: RUN_FAIL_ERR              = -20002

    ! -------------------------------------------------------
    ! MECHANISM ERROR CODES.
    ! -------------------------------------------------------

    ! Initialisation errors.
    Integer, Parameter :: LOAD_MECH_ERR             = -10001

    ! Process errors.
    Integer, Parameter :: PROCESS_RANGE_ERR         = -11001
    Integer, Parameter :: NEG_RATE_ERR              = -11002

    ! Process warnings/messages.
    Integer, Parameter :: FICTICIOUS_EVENT          = +11001

    ! Particle ensemble errors.
    Integer, Parameter :: ADD_PARTICLE_ERR          = -12001
    Integer, Parameter :: CHOOSE_PARTICLE_ERR       = -12002
    Integer, Parameter :: GET_PARTICLE_ERR          = -12003
    Integer, Parameter :: REMOVE_PARTICLE_ERR       = -12004
    Integer, Parameter :: REPLACE_PARTICLE_ERR      = -12005
    Integer, Parameter :: GET_LIST_ERR              = -12006
    Integer, Parameter :: SET_LIST_ERR              = -12007
    Integer, Parameter :: CLEAR_ENSEMBLE_ERR        = -12008
    Integer, Parameter :: GET_PROP_ERR              = -12009

    ! -------------------------------------------------------
    ! CHEMISTRY MODULE ERROR CODES.
    ! -------------------------------------------------------

    Integer, Parameter  ::  UNLOAD_CHEM_ERR         = -33001
    Integer, Parameter  ::  CHEM_INIT_ERR           = -33002

    Contains

    ! -------------------------------------------------------

    Subroutine TranslateErr(code, str)
        ! DESCRIPTION:
        !   Provides a string description of a given
        !   error code.
        
        Implicit None

        ! ARGUMENTS.
        ! Error code to translate.
        Integer, Intent(IN) :: code
        ! String to return.
        Character*(*), Intent(OUT) :: str

        ! EXECUTABLE CODE.
        Select Case (code)
            Case (0)
                str = "No error occurred."
            Case (MEMORY_ERR)
                str = "Error allocating or deallocating memory."
            Case (LOAD_MECH_ERR)
                str = "Error loading soot mechanism into Sweep."
            Case (PROCESS_RANGE_ERR)
                str = "An invalid index for a soot process was supplied."
            Case (NEG_RATE_ERR)
                str = "A negative process rate was encountered."
            Case (ADD_PARTICLE_ERR)
                str = "Failed to add a particle to the ensemble."
            Case (CHOOSE_PARTICLE_ERR)
                str = "Failed to select a particle from the ensemble."
            Case (GET_PARTICLE_ERR)
                str = "Failed to get a particle from the ensemble."
            Case (REMOVE_PARTICLE_ERR)
                str = "Failed to remove a particle from the ensemble."
            Case (REPLACE_PARTICLE_ERR)
                str = "Failed to replace a particle in the ensemble."
            Case (GET_LIST_ERR)
                str = "Failed to get a list of particle from the ensemble."
            Case (SET_LIST_ERR)
                str = "Failed to set the list of particle from the ensemble."
            Case (CLEAR_ENSEMBLE_ERR)
                str = "Failed to clear the ensemble of all particles."
            Case (GET_PROP_ERR)
                str = "Failed to get particle properties from the ensemble."
            Case (STEP_END_FAIL_ERR)
                str = "Stochastic step failed to reach stop time."
            Case (RUN_FAIL_ERR)
                str = "Stochastic loop failed to complete."
            Case (UNLOAD_CHEM_ERR)
                str = "Failed to get the chemistry from Sweep."
            Case Default
                str = "Unknown error."
        End Select
    End Subroutine

End Module

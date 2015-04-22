! *****************************************************************************
!
! File:					swpsoln.f90
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
!	Definition of types which describe data used by Sweep to solve a population
!   balance.
! *****************************************************************************

Module SWPSOLN
    Use SWPENSEMBLE, only: Ensemble
    Use SWPCHEM, only: ChemData
    Implicit None
    Public

    Type Solution
        Type(Ensemble) :: Ensemble  ! Particle ensemble.
        Type(ChemData) :: Chemistry ! Chemistry data.
        ! Scaling variables.
        Real :: BaseVolume      = -1.0 ! Base sample volume.
        Real :: BaseTemperature = -1.0 ! Temperature at which base volume applies.
        Real :: LastKnownScale  =  1.0 ! Last ensemble scaling before update.
        ! Solver parameters.
        Real :: SplitRatio = 1.0D9 ! Mean number of non-split events before updating all deferred.
        Integer :: MinStepCount = 1 ! Minimum number of split steps.
        ! Counters.
        Integer(4), Pointer :: ProcessCounter(:)    ! Counts of each process performed.
        Integer(4), Pointer :: FicticiousCounter(:) ! Counts of ficticious events.
        Integer(8) :: StepCount = 0         ! Stochastic step count.
        Integer(8) :: LPDACount = 0         ! Number of LPDA integrates done.
        ! Computation timing (s).
        Real :: CT=0.0, PreSplitCT=0.0, DeferCT=0.0, RatesCT=0.0, &
			    ProcessCT=0.0, ChemCT=0.0
        ! System parameters (should go somewhere ewlse eventually).
        Logical :: SolveInOut = .False. ! Solve inflow & outflow?
        Real    :: ResidenceTime = 0.0 ! PSR, required for inflow & outflow.
    End Type

    Contains

	! -------------------------------------------------------

	Subroutine Init(soln, pcount, spnames, mech, flag)
		! DESCRIPTION:
		!	Initialises a solution type.

        Use SWPENSEMBLE, InitEnsemble=>Init
        Use SWPCHEM
        Use SWPMECH_TYPES
		Implicit None
        
        ! ARGUMENTS.
        Type(Solution), Intent(OUT) :: soln       ! Solution to initialise.
        Integer, Intent(IN)         :: pcount     ! Number of particle to accomodate in ensemble.
        Character*(*), Intent(IN)   :: spnames(:) ! Species in gas-phase chemistry.
        Type(Mechanism), Intent(IN) :: mech       ! Mechanism to use.
        Integer, Intent(OUT)        :: flag       ! Error flag.

        ! EXECTUABLE CODE.
        flag = 0
        Call InitEnsemble(soln%Ensemble, pcount, flag)
        If (flag < 0) Return
        Call InitChem(soln%Chemistry, spnames, mech, flag)
        Allocate(soln%ProcessCounter(mech%ProcessCount), STAT=flag)
        soln%ProcessCounter = 0
        Allocate(soln%FicticiousCounter(mech%ProcessCount), STAT=flag)
        soln%FicticiousCounter = 0
	End Subroutine

	! -------------------------------------------------------

	Subroutine Reset(soln)
		! DESCRIPTION:
		!	Resets a solution type.

        Use SWPENSEMBLE, only: ClearEnsemble
		Implicit None
        
        ! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln

        ! EXECTUABLE CODE.
        Call ClearEnsemble(soln%Ensemble)
        soln%BaseVolume        = -1.0 
        soln%BaseTemperature   = -1.0 
        soln%LastKnownScale    =  1.0 
        soln%ProcessCounter    = 0
        soln%FicticiousCounter = 0
        soln%StepCount         = 0        
        soln%LPDACount         = 0        
        soln%CT                = 0.0
        soln%PreSplitCT        = 0.0
        soln%DeferCT           = 0.0
        soln%RatesCT           = 0.0
	    soln%ProcessCT         = 0.0
        soln%ChemCT            = 0.0
!        soln%SolveInOut        = .False.
!        soln%ResidenceTime     = 0.0
	End Subroutine

	! -------------------------------------------------------

	Subroutine CopySolution(s1, s2)
		! DESCRIPTION:
		!	Copies the data from one solution to another.

        Use SWPCHEM
        USe SWPENSEMBLE
		Implicit None
        
        ! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: s1, s2

        ! VARIABLES.
        Integer :: err = 0

        ! EXECTUABLE CODE.
        s2%ProcessCounter = s1%ProcessCounter
        s2%FicticiousCounter = s1%FicticiousCounter
        s2%BaseVolume = s1%BaseVolume
        s2%BaseTemperature = s1%BaseTemperature
        s2%LastKnownScale = s2%LastKnownScale
        s2%SplitRatio = s1%SplitRatio
        s2%MinStepCount = s1%MinStepCount
        s2%StepCount = s2%StepCount
        s2%LPDACount = s1%LPDACount
        s2%CT = s1%CT
        s2%PreSplitCT = s1%PreSplitCT
        s2%DeferCT = s1%DeferCT
        s2%RatesCT = s1%RatesCT
        s2%ProcessCT = s1%ProcessCT
        s2%ChemCT = s1%ChemCT
        s2%SolveInOut = s1%SolveInOut
        s2%ResidenceTime = s1%ResidenceTime
        Call CopyChem(s1%Chemistry, s2%Chemistry)
        Call CopyEnsemble(s1%Ensemble, s2%Ensemble)
	End Subroutine

	! -------------------------------------------------------

	Subroutine DeleteSolution(soln)
		! DESCRIPTION:
		!	Free memory allocated to a solution.

        Use SWPENSEMBLE, InitEnsemble=>Init
        Use SWPCHEM
		Implicit None
        
        ! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln ! Solution to initialise.

        ! VARIABLES.
        Integer :: err = 0

        ! EXECTUABLE CODE.
        Deallocate(soln%ProcessCounter, soln%FicticiousCounter, STAT=err)
        Call DeleteChem(soln%Chemistry)
        Call DeleteEnsemble(soln%Ensemble)
	End Subroutine

	! -------------------------------------------------------

	Subroutine WriteSolution(soln, fnum)
		! DESCRIPTION:
		!	Writes a solution to a binary file.
        Use SWPENSEMBLE, only: WriteEnsemble
        Use SWPCHEM, only: WriteChemistry
		Implicit None

        ! ARGUMENTS.
        Type(Solution), Intent(IN) :: soln ! Solution to output.
        Integer, Intent(IN)        :: fnum ! UNIT number on which to output.

        ! EXECUTABLE CODE.
        Call WriteEnsemble(soln%Ensemble, fnum)
        Call WriteChemistry(soln%Chemistry, fnum)
        Write(fnum) soln%BaseVolume, soln%BaseTemperature, soln%LastKnownScale, &
                    soln%ProcessCounter, soln%FicticiousCounter, soln%StepCount, &
                    soln%LPDACount, soln%CT, soln%PreSplitCT, soln%DeferCT, &
                    soln%RatesCT, soln%ProcessCT, soln%ChemCT, soln%SolveInOut, &
                    soln%ResidenceTime
	End Subroutine

	! -------------------------------------------------------

	Subroutine ReadSolution(soln, fnum)
		! DESCRIPTION:
		!	Reads a solution from a binary file.
        Use SWPENSEMBLE, only: ReadEnsemble
        Use SWPCHEM, only: ReadChemistry
		Implicit None

        ! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln ! Solution to output.
        Integer, Intent(IN)           :: fnum ! UNIT number on which to output.

        ! EXECUTABLE CODE.
        Call ReadEnsemble(soln%Ensemble, fnum)
        Call ReadChemistry(soln%Chemistry, fnum)
        Read(fnum) soln%BaseVolume, soln%BaseTemperature, soln%LastKnownScale, &
                   soln%ProcessCounter, soln%FicticiousCounter, soln%StepCount, &
                   soln%LPDACount, soln%CT, soln%PreSplitCT, soln%DeferCT, &
                   soln%RatesCT, soln%ProcessCT, soln%ChemCT, soln%SolveInOut, &
                   soln%ResidenceTime
	End Subroutine

	! -------------------------------------------------------

	Pure Real Function SampleVolume(soln, T)
		! DESCRIPTION:
		!	Calculates the sample volume for a given
		!	temperature using the base volume and temperature
		!	taken from the start of the simulation.  This function
		!	uses Boyle's law, and is only suitable for ideal gases.
		! RETURNS:
		!	Sample volume (cm-3).
        Use SWPENSEMBLE, only: ScalingFactor
		Implicit None
        Type(Solution), Intent(IN) :: soln
		Real, Intent(IN) :: T	! Temperature (K).
		SampleVolume = ScalingFactor(soln%Ensemble) * soln%BaseVolume
!		SampleVolume = ScalingFactor(soln%Ensemble) * soln%BaseVolume * T / soln%BaseTemperature
	End Function

	! -------------------------------------------------------

	Subroutine SetScaling(soln, maxM0, temp)
		! DESCRIPTION:
		!	Sets the base sample volume and the corresponding
		!	gas temperature (K).
		Implicit None
        Type(Solution), Intent(INOUT) :: soln
		Real, Intent(IN) :: maxM0 ! Maximum expected number density.
        Real, Intent(IN) :: temp  ! Temperature (K) at predicted max. M0.
		soln%BaseVolume	= Real(soln%Ensemble%Capacity) / maxM0
		soln%BaseTemperature = temp
	End Subroutine

	! -------------------------------------------------------

	Subroutine SetM0(soln, m0)
		! DESCRIPTION:
		!	Sets the M0.
		Implicit None
        Type(Solution), Intent(INOUT) :: soln
		Real, Intent(IN) :: m0 ! Maximum expected number density.
        If (m0 > 0.0E0) Then
		    soln%BaseVolume	= Real(soln%Ensemble%FirstSpace-1) / m0
        End If
	End Subroutine

	! -------------------------------------------------------

	Real Function GetM0(soln)
		! DESCRIPTION:
		!	Gets the M0.
		Implicit None
        Type(Solution), Intent(INOUT) :: soln
		GetM0 = Real(soln%Ensemble%FirstSpace-1) / soln%BaseVolume
	End Function

	! -------------------------------------------------------

	Subroutine SetBaseTemperature(soln, temp, scaleVol)
		! DESCRIPTION:
		!	Sets the base gas temperature (K) of the system.
		!	This involve rescaling the sample volume first.

		Implicit None

		! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln
		Real, Intent(IN) :: temp ! Temperature (K) to set.
		! Set to True to rescale the sample volume to
		! the new temperature, False to keep the sample
		! volume constant.
		Logical, Intent(IN)	::	scaleVol

		! EXECUTABLE CODE.
		If (scaleVol) soln%BaseVolume = SampleVolume(soln, temp)
		soln%BaseTemperature = temp
	End Subroutine

	! -------------------------------------------------------

	Subroutine ResetCounters(soln)
		! DESCRIPTION:
		!	Resets the process counters to 0.
		Implicit None
        Type(Solution), Intent(INOUT) :: soln
		soln%ProcessCounter = 0
        soln%FicticiousCounter = 0
        soln%StepCount = 0
        soln%LPDACount = 0
	End Subroutine

	! -------------------------------------------------------

	Subroutine ResetCT(soln)
		! DESCRIPTION:
		!	Resets the process counters to 0.
		Implicit None
        Type(Solution), Intent(INOUT) :: soln
        soln%CT         = 0.0
        soln%PreSplitCT = 0.0
        soln%DeferCT    = 0.0
        soln%RatesCT    = 0.0
		soln%ProcessCT  = 0.0
        soln%ChemCT     = 0.0
	End Subroutine

End Module

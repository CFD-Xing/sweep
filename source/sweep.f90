! *****************************************************************************
!
! File:				sweep.f90
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
! Prerequisite source code and licencing:
!   In addition to the mops source code, the following code libraries are
!   required to successfully compile this program:
!
!       xmlf90 library (revised BSD licence, included here):
!           Developed by Prof. Alberto Garcia et al. of Universidad del
!           Pais Vasco, Bilbao, Spain (http://lcdx00.wm.lc.ehu.es/ag/xml/).
!
!       MersenneTwister random number generator library (GNU Library GPL):
!           Developed by Makoto Matsumoto and Takuji Nishimura.
!
! Functions:
!    InitSweep		-	Initialises Sweep.
! *****************************************************************************


Module Sweep
    ! -------------------------------------------------------
    ! SWEEP PUBLIC INTERFACE.
    ! -------------------------------------------------------
    ! All these functions/modules will be visible to any
    ! code which uses the Sweep solver.  That is, this list
    ! represents the public interface of Sweep2.
    Use SWPERR
    Use SWPSTATS
    Use SWPPARAMS
    Use SWPPART, GetParticleStats=>GetStats
    Use SWPSTEP, only: RunSweep=>Run
    Use SWPCHEM, only: LoadSweepChem=>LoadChem, UnloadSweepChem=>UnloadChem, &
                       SweepChemData=>ChemData, InitSweepChem=>InitChem, &
                       GetSweepChem=>GetChem
    Use SWPENSEMBLE, only: ClearEnsemble, GetEnsembleArray, &
                           LoadEnsembleArray, ParticleCount
    Use SWPSOLN, only: SetScaling, SampleVolume, SetM0, GetM0
    Use SWPSOLN, only: SweepSolution=>Solution, InitSweepSolution=>Init, &
                       DeleteSweepSolution=>DeleteSolution, ResetSweepSolution=>Reset, &
                       ResetSweepCT=>ResetCT, WriteSweepSolution=>WriteSolution, &
                       ReadSweepSolution=>ReadSolution, CopySweepSolution=>CopySolution
    Use SWPMECH, only: SweepMechanism=>Mechanism, InitSweepMechanism=>InitMech, &
                       DeleteSweepMechanism=>DeleteMech
    Use SWPCHEM_SHARED, only: ABFAlpha
    Use SWPPROCESS, only: SweepChemChangeRates=>ChemChangeRates, SweepRateTerms=>GroupedRateTerms
    ! -------------------------------------------------------

    Implicit none
    Public

    Contains

    ! -------------------------------------------------------

    Subroutine InitSweep()
        ! DESCRIPTION:
        !   Initialises Sweep.
        Use SWPRNG
        Implicit None
        Call SeedRnd(123) ! Seed random number generator.
    End Subroutine
End Module

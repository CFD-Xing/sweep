! *****************************************************************************
!
! File:				profile_driver_program.f90
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
!   A program to run sweep using profiled output from premix.
!
! *****************************************************************************

Program Profile_Driver_Program
    Use Sweep
    Use Profile_Driver_Input
    Use Profile_Driver_Output
    Implicit None

    ! PARAMETERS.
    Character(LEN=*), Parameter :: setfile = "sweep.inp" ! Settings file name.

    ! VARIABLES.
    Type(Settings)    :: sim  ! Simulation settings.
    Type(ChemProfile) :: chem ! Chemistry profile.

    ! Sweep variables.
    Type(SweepSolution)  :: soln
    Type(SweepMechanism) :: mech

    ! Counters and time step variables.
    Integer :: err, i, irun, iint
    Real :: t, tstop, dt, start, finish

    ! EXECUTABLE CODE.
    call cpu_time(start)

    ! Load settings file.
    Call CommandLine(sim)
    Call ReadAllSettings(setfile, sim, err)
    If (err /= 0) Then
        Print *, "FAILED to read settings file!  Quitting..."
        Stop
    End If

    ! Read chemistry profile.
    Call ReadChemProfile(sim%ChemFile, chem, err)
    If (err /= 0) Then
        Print *, "FAILED to read chemistry profile!  Quiting..."
        Stop
    End If

    ! Initialise Sweep.
    Call InitSweep()

    Call InitSweepMechanism(mech, sim%MechFile, chem%Header(1:chem%iT-1), err)
    mech%ParticleModel = sim%ParticleModel
    If (sim%SetActSurf) Then
        mech%ActSurfModel  = sim%ActSurfModel
        If (sim%ActSurfModel == ACT_SURF_CONST) mech%ConstActSurf = sim%ConstActSurf
    End If
    If (err /= 0) Then
        Print *, "FAILED to initialise Sweep mechanism!  Quiting..."
        Stop
    End If

    Call InitSweepSolution(soln, sim%PCount, chem%Header(1:chem%iT-1), mech, err)
    soln%SplitRatio = sim%SRatio
    If (err /= 0) Then
        Print *, "FAILED to initialise Sweep solution!  Quiting..."
        Stop
    End If

    ! Put chemistry into Sweep.
    Call LoadSweepChem(soln%Chemistry, chem%Times, chem%Chem)

    ! Solve system.
    If (sim%Solve) Then
        Do irun = 1, sim%Runs
            Print "(X,A,I3,A,I3)", "Run ", irun, " / ", sim%Runs

            ! Reset Sweep solution.
            Call ClearEnsemble(soln%Ensemble)
            ! Call SetScaling(soln, sim%MaxM0, MaxVal(chem%Chem(chem%iT,:)))
            Call SetScaling(soln, sim%MaxM0, MinVal(chem%Chem(chem%irho,:)))

            ! Start output.
            t = sim%Times%Items(1)
            Call BeginOutput(sim%OutputFile, irun, err)
            Call PrintHeader()
            Call WriteOutput(0, t, soln, mech, sim%Range, err)

            Do iint = 1, sim%Steps%Count
                Call PrintHeader()

                t     = Real(sim%Times%Items(iint))
                tstop = Real(sim%Times%Items(iint+1))
                dt    = (tstop - t) / Real(sim%Steps%Items(iint))

                Do i = 1, sim%Steps%Items(iint)
                    ! Run Sweep.
                    tstop = t + dt
                    Call RunSweep(t, tstop, soln, mech, err)
                    Call WriteOutput(i, t, soln, mech, sim%Range, err)
                End Do

                Call WritePSL(sim%OutputFile, irun, iint, t, soln, mech, err)
            End Do

            Call EndOutput()
        End Do
    End If

    ! Process output files.
    If (sim%PostProcess) Then
        Call ProcessOutput(sim%OutputFile, sim%Runs, Sum(sim%Steps%Items), mech, err)
        Call ProcessPSLs(sim%OutputFile, sim%Runs, sim%Steps%Count, soln%Ensemble%Capacity, mech, err)
    End If

    ! Clear memory.
    Call DeleteSweepSolution(soln)
    Call DeleteSweepMechanism(mech)
    ! CPU time
    call cpu_time(finish)
    print '("Time = ",f4.0," seconds.")',finish-start
End Program


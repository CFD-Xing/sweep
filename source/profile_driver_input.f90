! *****************************************************************************
!
! File:					profile_driver_input.f90
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
!   Input routines for the Sweep profiled chemistry driver program.
!
! Functions:
!
! *****************************************************************************

Module Profile_Driver_Input
    Use LIST_T
    Use Sweep
    Implicit None
    Public
    
    Integer, Parameter :: UIN = 102
    Integer, Parameter :: MAX_POINTS=1000, MAX_VARS=200

    Type Settings
        ! Control parameters.
        Logical :: Solve=.True., PostProcess=.True.
        ! Input files.
        Character(LEN=200) :: MechFile=""  ! Particle mechanism file name.
        Character(LEN=200) :: ChemFile=""  ! Chemistry profile file name.
        ! Time stepping variables.
        Type(LIST_DP)  :: Times  ! Times between PSL outputs (first entry is start time).
        Type(LIST_INT) :: Steps  ! Steps between output times.
        Integer        :: Runs=0 ! Number of runs to perform.
        ! Sweep settings.
        Integer :: PCount=0     ! Max. stochastic particle count in Sweep.
        Real    :: MaxM0=0.0    ! Max. expected number density in system.
        Real    :: SRatio=1.0E9 ! Split ratio.
        Integer :: ParticleModel = SPHERICAL_PARTICLE_MODEL
        Logical :: SetActSurf    = .False.
        Integer :: ActSurfModel  = ACT_SURF_CONST
        Real    :: ConstActSurf  = 1.0E0
        ! Output settings.
        Character(LEN=200) :: OutputFile="" ! Output file name.
        Type(StatRange)    :: Range ! Range of output statistics.
    End Type

    Type ChemProfile
        Integer :: NVARS=0, NPOINTS=0
        Character(LEN=16), Pointer :: Header(:)
        Real, Pointer :: Times(:), Chem(:,:)
        Integer :: iT=0, iP=0, iA=0
    End Type

    Contains

	! -------------------------------------------------------

    Subroutine CommandLine(sim)
        ! DESCRIPTION.
        !   Reads program settings from the command line.

#ifdef _WIN32
	    Use DFLIB, Only: NARGS, GETARG
#endif
        Implicit None

        ! ARGUMENTS.
        Type(Settings), Intent(INOUT) :: sim

        ! VARIABLES.
        Integer(2) :: N, i
        Character(LEN=128) :: arg

        ! EXECUTABLE CODE.
!        N = NARGS()
        N = 1
        i = 1

        Do While (i <= N)
            Call GetArg(i, arg)
            Select Case (arg)
                Case ("-po","/po","\po")
                    sim%Solve = .False.
                    sim%PostProcess = .True.
            End Select
            i = i + 1
        End Do
    End Subroutine
    ! -------------------------------------------------------

    Subroutine ReadAllSettings(file, set, flag)
	    ! DESCRIPTION:
	    !	Reads a standard settings file and returns 
        !   all the program settings for the desired
        !   simulation.

        Use SettingsFileReader
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)  :: file ! File to read.
        Type(Settings), Intent(INOUT) :: set  ! Somewhere to put settings.
        Integer, Intent(OUT)          :: flag ! Error flag.

        ! VARIABLES.
        Integer :: i, int, err
        Double Precision :: dp
        Character(LEN=100) :: str
        Type(FileData) :: fdata

        ! EXECUTABLE CODE.
        flag = 0

        ! Load settings file.
        Call LoadSettingsFile(file, UIN, fdata, err)
        If (err /= 0) Return

        ! Input files.
        i = GetParams(fdata, "MECHFILE", set%MechFile, err)
        If (err /= 0) Then
            Print *, "FAILED to read mechanism file name!"
            flag = -1
            Return
        End If
        i = GetParams(fdata, "CHEMFILE", set%ChemFile, err)
        If (err /= 0) Then
            Print *, "FAILED to read chemistry file name!"
            flag = -1
            Return
        End If

        ! Time step parameters.
        i = GetParams(fdata, "STARTT", dp, err)
        set%Times = set%Times + dp

        i = 1
        Do While(err == 0)
            i = GetParams(fdata, "TIME", "RI", dp, int, err, i)
            If (err /= 0) Exit
            set%Times = set%Times + dp
            set%Steps = set%Steps + int
            i = i + 1
        End Do

        i = GetParams(fdata, "RUNS", set%Runs, err)
        If (err /= 0) Then
            Print *, "FAILED to read run count!"
            flag = -1
            Return
        End If

        ! Sweep parameters.

        i = GetParams(fdata, "PCOUNT", set%PCount, err)
        If (err /= 0) Then
            Print *, "FAILED to read stochastic particle count!"
            flag = -1
            Return
        End If

        i = GetParams(fdata, "MAXM0", dp, err)
        If (err /= 0) Then
            Print *, "FAILED to read max. likely M0!"
            flag = -1
            Return
        Else
            set%MaxM0 = Real(dp)
        End If

        i = GetParams(fdata, "SRATIO", dp, err)
        If (err /= 0) Then
            set%SRatio = 1.0D9
            Print *, "DEFAULT split ratio"
        Else
            set%SRatio = Real(dp)
        End If

        i = GetParams(fdata, "PMODEL", str, err)
        If (err /= 0) Then
            Print *, "DEFAULT particle model"
        Else
            Select Case (str)
                Case ("spherical")
                    set%ParticleModel = SPHERICAL_PARTICLE_MODEL
                Case ("surfvol")
                    set%ParticleModel = SURFACE_VOLUME_MODEL
                Case Default
                    Print *, "DEFAULT particle model"
            End Select
        End If

        i = GetParams(fdata, "ACTSURFMODEL", str, err)
        If (err /= 0) Then
            Print *, "DEFAULT active surface model"
        Else
            set%SetActSurf = .True.
            Select Case (str)
                Case ("const")
                    set%ActSurfModel = ACT_SURF_CONST
                Case ("profile")
                    set%ActSurfModel = ACT_SURF_PROFILE
                Case Default
                    Print *, "DEFAULT active surface model"
            End Select
        End If

        i = GetParams(fdata, "ACTSURF", dp, err)
        If (err /= 0) Then
            Print *, "DEFAULT active surface value"
        Else
            set%ConstActSurf = dp
        End If

        ! Output parameters.
        i = GetParams(fdata, "OUTPUTFILE", set%OutputFile, err)
        If (err /= 0) Then
            Print *, "FAILED to read output file name!"
            flag = -1
            Return
        End If

        i = GetParams(fdata, "LOWSTAT", set%Range%Lower, err)
        If (err /= 0) Then
            set%Range%Lower = 0.0D0
            Print *, "DEFAULT lower stats bound."
        End If

        i = GetParams(fdata, "HIGHSTAT", set%Range%Upper, err)
        If (err /= 0) Then
            set%Range%Upper = 1.0D300
            Print *, "DEFAULT upper stats bound."
        End If
        
        Print "(X)"
    End Subroutine
    
    ! -------------------------------------------------------

    Subroutine ReadChemProfile(file, chem, flag)
	    ! DESCRIPTION:
	    !	Reads a comma separated file of a flame profile
        !   generated by PREMIX.  These files used to be called
        !   PSDF_input.dat.  The following data columns should be
        !   present in the input file, in the given order:
        !
        !   1.  Time (s).
        !   2.  C2H2 conc'n (mol/cm3).
        !   3.  H2 conc'n (mol/cm3).
        !   4.  H conc'n (mol/cm3).
        !   5.  O2 conc'n (mol/cm3).
        !   6.  OH conc'n (mol/cm3).
        !   7.  H2O conc'n (mol/cm3).
        !   8.  CO conc'n (mol/cm3).
        !   9.  A4 conc'n (mol/cm3).
        !  10.  Temperature (K).
        !  11.  Pressure (bar).
        !  12.  Alpha (active site fraction).

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)   :: file ! Name of input file.
        Type(ChemProfile), Intent(OUT) :: chem ! Output chemistry profile.
        Integer, Intent(OUT)           :: flag ! Error flag, returns <0 on error.

        ! VARIABLES.
        Character(LEN=16) :: h1
        Integer :: i, reqT, reqP, reqA
        Real, Allocatable :: col(:)

        ! EXECUTABLE CODE.
        flag = 0

        ! Open input file.
        Open (UNIT=UIN, FILE=file, FORM="FORMATTED", STATUS="OLD", ACTION="READ", IOSTAT=flag)
        If (flag /= 0) Then
            flag = -1
            Return
        End If
        
        ! Read the number of points and variables from the 
        ! file.
        Read (UNIT=UIN, FMT=*) chem%NPOINTS, chem%NVARS

        ! Allocate profile space.
        Deallocate(chem%Header, chem%Times, chem%Chem, STAT=flag)
        Allocate(chem%Header(chem%NVARS), chem%Times(chem%NPOINTS), &
                 chem%Chem(chem%NVARS,chem%NPOINTS), STAT=flag)

        ! Read header row.
        Read (UNIT=UIN, FMT=*) h1, chem%Header

        ! Read profile.
        Do i = 1, chem%NPOINTS
            Read (UNIT=UIN, FMT=*) chem%Times(i), chem%Chem(:,i)
        End Do

        ! Close file.
        Close (UNIT=UIN, IOSTAT=flag)

        ! Sort profile.
        Allocate(col(chem%NPOINTS), STAT=flag)
        chem%iT = IndexOf("T", chem%Header)
        chem%iP = IndexOf("P", chem%Header)
        chem%iA = IndexOf("Alpha", chem%Header)

        If (chem%iA == 0) Then
            ! No Alpha in profile.
            reqT = chem%NVARS - 1
            reqP = chem%NVARS
            reqA = 0
        Else
            reqT = chem%NVARS - 2
            reqP = chem%NVARS - 1
            reqA = chem%NVARS
        End If

        ! Put temperature column in the correct place.
        If (chem%iT /= reqT) Then
            col = chem%Chem(reqT,:)
            chem%Chem(reqT,:) = chem%Chem(chem%iT,:)
            chem%Chem(chem%iT,:) = col
            h1 = chem%Header(reqT)
            chem%Header(reqT) = chem%Header(chem%iT)
            chem%Header(chem%iT) = h1

            ! Check to see if we have displaced P
            ! or Alpha columns.
            If (chem%iP == reqT) Then
                chem%iP = chem%iT
            Else
                If (chem%iA == reqT) chem%iA = chem%iT
            End If
            chem%iT = reqT
        End If

        ! Put pressure column in the correct place.
        If (chem%iP /= reqP) Then
            col = chem%Chem(reqP,:)
            chem%Chem(reqP,:) = chem%Chem(chem%iP,:)
            chem%Chem(chem%iP,:) = col
            h1 = chem%Header(reqP)
            chem%Header(reqP) = chem%Header(chem%iP)
            chem%Header(chem%iP) = h1

            ! Check to see if we have displaced
            ! alpha column.
            If (chem%iA == reqP) chem%iA = chem%iP
            chem%iP = reqP
        End If

        ! Put alpha column the correct place.
        If (chem%iA /= reqA) Then
            col = chem%Chem(reqA,:)
            chem%Chem(reqA,:) = chem%Chem(chem%iA,:)
            chem%Chem(chem%iA,:) = col
            h1 = chem%Header(reqA)
            chem%Header(reqA) = chem%Header(chem%iA)
            chem%Header(chem%iA) = h1
            chem%iA = reqA        
        End If

        ! Sorted!
        Deallocate(col, STAT=flag)
    End Subroutine
    
End Module


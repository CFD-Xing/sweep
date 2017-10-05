! *****************************************************************************
!
! File:					profile_driver_output.f90
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
!   Output routines for the Sweep profiled chemistry driver program.
!
! Functions:
!
! *****************************************************************************

Module Profile_Driver_Output
    Implicit None
    Public
    
    Integer, Parameter :: UOUT = 102, UPSL = 103
    Character(LEN=*), Parameter :: AFMT   = "(A,400(',',A,:))"
    Character(LEN=*), Parameter :: ESFMT  = "(ES,400(',',ES,:))"

    ! Confidence interval parameter.
    Double Precision, Parameter :: CONFA = 3.29D0 ! 99.9% confidence interval

    Contains

    ! -------------------------------------------------------

    Subroutine PrintHeader()
        Use ConsoleTable, ph=>PrintHeader
        Implicit None
	    Call ph("Step", (/"Time (s)", "#SP", "M0", "M1", "Surface", "Avg. D"/), 6)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine BeginOutput(file, run, flag)
	    ! DESCRIPTION:
	    !	Opens a file for output, appending the given run
        !   number to the file name.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: file ! File name root (no extension).
        Integer, Intent(IN)          :: run  ! Run number to append to file name.
        Integer, Intent(OUT)         :: flag ! Error flag.

        ! VARIABLES.
        Logical :: o

        ! EXECUTABLE CODE.
        flag = 0
        
        Inquire(UNIT=UOUT, OPENED=o)

        If (o == .False.) Then
		Open(UNIT=UOUT, FILE=Trim(file)//"_"//Trim(CStr(run))//".dat", FORM="UNFORMATTED", &
		     STATUS="REPLACE", ACTION="READWRITE", IOSTAT=flag)
             If (flag /= 0) flag = -1
        Else
            Print *, "Output file already open!"
            flag = -2
        End If
    End Subroutine

    ! -------------------------------------------------------

    Subroutine WriteOutput(N, t, soln, mech, range, flag)
	    ! DESCRIPTION:
	    !	Writes unformatted output data to file and
        !   writes some values to the console.

        Use Sweep
        Use SWPCHEM, only: GetTemperature,GetRho
        Use SWPSTATS_INDICES
        Use ConsoleTable
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)              :: N     ! Step number.
        Real, Intent(IN)                 :: t     ! Current time (s).
        Type(SweepSolution), Intent(IN)  :: soln  ! Solution to output.
        Type(SweepMechanism), Intent(IN) :: mech  ! Mechanism to use.
        Type(StatRange), Intent(IN)      :: range ! Stats output range.
        Integer, Intent(OUT)             :: flag  ! Error flag.

        ! VARIABLES.
        Integer :: i, err
        Double Precision :: stats(NSTATS), rates(mech%GroupCount)

        ! EXECUTABLE CODE.
        flag = 0
        
      !  stats = GetStats(soln, GetTemperature(soln%Chemistry, Real(t), Real(t)), range)
        stats = GetStats(soln, GetRho(soln%Chemistry, Real(t), Real(t)), range)
        rates = Dble(SweepRateTerms(Real(t), GetSweepChem(soln%Chemistry, Real(t), Real(t), err), soln, mech, flag))

        Call PrintScientificD(N, (/Dble(t), stats(iSP), stats(iM0), stats(iM1), stats(iSurf), stats(iDiam)/), 6)
        Write(UOUT) t, stats, rates
    End Subroutine

    ! -------------------------------------------------------

    Subroutine EndOutput()
	    ! DESCRIPTION:
	    !	Closes the current output file.
        Implicit None
        Integer :: err
        Close (UNIT=UOUT, IOSTAT=err)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine ProcessOutput(file, runs, N, mech, flag)
	    ! DESCRIPTION:
	    !	Processes unformatted output files for
        !   multiple runs and returns averaged data in a CSV
        !   file.

        Use Sweep
        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)     :: file ! File name root (no extension).
        Integer, Intent(IN)              :: runs ! Number of runs.
        Integer, Intent(IN)              :: N    ! Number of output steps per file.
        Type(SweepMechanism), Intent(IN) :: mech ! Mechanism to use.
        Integer, Intent(OUT)             :: flag ! Error flag.

        ! VARIABLES.
        Integer :: i, irun
        Double Precision :: t(N), druns
        Double Precision :: stats(NSTATS), as(NSTATS,N), es(NSTATS,N)
        Double Precision :: rates(mech%GroupCount), ar(mech%GroupCount,N), er(mech%GroupCount,N)

        ! EXECUTABLE CODE.
        flag = 0
        
        as = 0.0D0
        ar = 0.0D0
        es = 0.0D0
        er = 0.0D0

        Do irun = 1, runs
		    Open(UNIT=UOUT, FILE=Trim(file)//"_"//Trim(CStr(irun))//".dat", FORM="UNFORMATTED", &
		         STATUS="OLD", ACTION="READ", IOSTAT=flag)

            If (flag /= 0) Then
                flag = -1
                Return
            End If

            Do i = 1, N
                Read(UOUT, IOSTAT=flag) t(i), stats, rates

                ! Calculate sums.
                as(:,i) = as(:,i) + stats
                ar(:,i) = ar(:,i) + rates

                ! Calculate sums of squares.
                If (runs > 1) Then
                    es(:,i) = es(:,i) + (stats * stats)
                    er(:,i) = er(:,i) + (rates * rates)
                End If
            End Do

            Close(UNIT=UOUT, IOSTAT=flag)
        End Do

        ! Calculate averages over all runs.
        druns = Dble(runs)
        as    = as / druns
        ar    = ar / druns

        If (runs > 1) Then
            ! Calculate standard deviations over all runs.
            es = (es / druns) - (as * as)
            er = (er / druns) - (ar * ar)
            ! Calculate confidence intervals.
            es = CONFA * Sqrt(Abs(es) / druns)
            er = CONFA * Sqrt(Abs(er) / druns)
        End If

        ! Write results to CSV file.
		Open(UNIT=UOUT, FILE=Trim(file)//".csv", FORM="FORMATTED", &
		     STATUS="REPLACE", ACTION="READWRITE", IOSTAT=flag)
    
        Write(UOUT, FMT=AFMT) "Step", "Time (s)", &
                              (Trim(StatNames(i)), i=1, NSTATS), &
                              (Trim(mech%GroupNames(i)), i=1, mech%GroupCount), &
							  (Trim(StatNames(i)) // " Err", i=1, NSTATS), &
                              (Trim(mech%GroupNames(i)) // " Err", i=1, mech%GroupCount)
        Do i = 1, N
            Write(UOUT, FMT=ESFMT) Dble(i-1), t(i), as(:,i), ar(:,i), es(:,i), er(:,i)
        End Do

        Close(UOUT, IOSTAT=flag)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine WritePSL(file, run, num, t, soln, mech, flag)
	    ! DESCRIPTION:
	    !	Writes unformatted output data to file and
        !   writes some values to the console.

        Use Sweep
        Use SWPCHEM, only: GetTemperature,GetRho
        Use SWPSTATS_INDICES
        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)     :: file ! File name root (no extension).
        Integer, Intent(IN)              :: run  ! Run number to append to file name.
        Integer, Intent(IN)              :: num  ! PSL number to append to file name.
        Real, Intent(IN)                 :: t    ! Time (s).
        Type(SweepSolution), Intent(IN)  :: soln ! Solution to output.
        Type(SweepMechanism), Intent(IN) :: mech ! Mechanism to use.
        Integer, Intent(OUT)             :: flag ! Error flag.

        ! VARIABLES.
        Type(StatRange)  :: range
        Double Precision :: psl(PSL_VAR_COUNT,soln%Ensemble%Capacity), &
                            pcl(mech%ComponentCount+mech%TrackerCount,soln%Ensemble%Capacity), &
                            stats(NSTATS)

        ! EXECUTABLE CODE.
		Open(UNIT=UPSL, FILE=Trim(file)//"_"//Trim(CStr(run))//"-psl_"//Trim(CStr(num))//".dat", &
             FORM="UNFORMATTED", STATUS="REPLACE", ACTION="READWRITE", IOSTAT=flag)
       ! stats = GetStats(soln, GetTemperature(soln%Chemistry, t, t), range)
        stats = GetStats(soln, GetRho(soln%Chemistry, t, t), range)
        psl = CompilePSL(soln%Ensemble)
        pcl = CompilePCL(soln%Ensemble, mech)
        Write(UPSL) stats(iM0), psl, pcl
        Close(UPSL, IOSTAT=flag)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine ProcessPSLs(file, runs, nums, N, mech, flag)
	    ! DESCRIPTION:
	    !	Writes unformatted output data to file and
        !   writes some values to the console.

        Use Sweep
        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)     :: file ! File name root (no extension).
        Integer, Intent(IN)              :: runs ! Run number to append to file name.
        Integer, Intent(IN)              :: nums ! PSL number to append to file name.
        Integer, Intent(IN)              :: N    ! Ensemble capacity.
        Type(SweepMechanism), Intent(IN) :: mech ! Mechanism to use.
        Integer, Intent(OUT)             :: flag ! Error flag.

        ! VARIABLES.
        Integer :: irun, inum, i, j, pcount
        Double Precision :: M0(runs), psl(PSL_VAR_COUNT,N,runs), wt, &
                            pcl(mech%ComponentCount+mech%TrackerCount,N,runs)
	    Character(LEN=*), Parameter :: M0FMT = "(ES,400(',,,',ES,:))"

        ! EXECUTABLE CODE.
        Do inum = 1, nums
            ! Compile PSL & PCL for all runs.
            M0  = 0.0D0
            psl = 0.0D0
            Do irun = 1, runs
		        Open(UNIT=UPSL, FILE=Trim(file)//"_"//Trim(CStr(irun))//"-psl_"//Trim(CStr(inum))//".dat", &
                     FORM="UNFORMATTED", STATUS="OLD", ACTION="READ", IOSTAT=flag)
                Read(UPSL) M0(irun), psl(:,:,irun), pcl(:,:,irun)
                Close(UPSL, IOSTAT=flag)
            End Do

            ! Write PSL to CSV file.
		    Open(UNIT=UPSL, FILE=Trim(file)//"-psl_"//Trim(CStr(inum))//".csv", &
                 FORM="FORMATTED", STATUS="REPLACE", ACTION="READWRITE", IOSTAT=flag)
!            Write(UPSL, FMT=M0FMT) M0(1:runs)
!            Write(UPSL, FMT=AFMT) ("v", "s", "d", i=1, runs)
!            Do i = 1, N
!                Write(UPSL, FMT=ESFMT) Reshape(psl(:,i,1:runs),(/3*runs/))
!            End Do
            Write(UPSL, FMT=AFMT) "wt", "vol", "surf", "diam", "type", &
                                  (Trim(mech%Components(i)%Symbol), i=1, mech%ComponentCount), &
                                  (Trim(mech%Trackers(i)), i=1, mech%TrackerCount)
            Do i = 1, runs
                pcount = Count(psl(1,:,i) > 0.0D0)
                wt = M0(i) / Dble(pcount * runs)
                Do j = 1, pcount
                    Write(UPSL, FMT=ESFMT) wt, psl(1:4,j,i), pcl(:,j,i)
                End Do
            End Do
            Close(UPSL, IOSTAT=flag)

            ! Write PCL to CSV file.
!		    Open(UNIT=UPSL, FILE=Trim(file)//"-pcl("//Trim(CStr(inum))//").csv", &
!                 FORM="FORMATTED", STATUS="REPLACE", ACTION="READWRITE", IOSTAT=flag)
!            Write(UPSL, FMT=AFMT) 
!            Write(UPSL, FMT=AFMT) "vol", "surf", "diam", "wt"
!            Do i = 1, runs
!                pcount = Count(psl(3,:,i) > 0.0D0)
!                Do j = 1, pcount
!                    Write(UPSL, FMT=ESFMT) psl(:,j,i), invM0(i)
!                End Do
!            End Do
!            Do i = 1, N
!                Write(UPSL, FMT=ESFMT) Reshape(pcl(:,i,1:runs),(/(mech%ComponentCount+mech%TrackerCount)*runs/))
!            End Do
!            Close(UPSL, IOSTAT=flag)
        End Do
    End Subroutine
    
End Module

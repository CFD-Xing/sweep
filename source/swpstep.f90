! *****************************************************************************
!
! File:				swpstep.f90
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
!   Herein reside the time stepping routines for Sweep.  This module controls
!   the flow of the stochastic alogorithm.
!
! Functions:
!   Run		-	Runs Sweep over the given interval with the current
!				ensemble and gas-phase chemistry.
!   TimeStep	-	Performs a single time step and returns the delta t.
!   CurrentTime	-	Returns the current flow time within the simulation.
!	---------------------------------------------------------------------------
!   ResetCT	-	Resets computation timing variables.
!   MarkCT	-	Sets the current time as the start time for the next
!				computation timing interval.
!   IncCT	-	Increments the given computation timing variable.
! *****************************************************************************

Module SWPSTEP
    ! -------------------------------------------------------
    ! IMPORT PUBLIC MODULE FUNCTIONS.
    ! -------------------------------------------------------
    Use SWPERR ! Error codes.
    ! -------------------------------------------------------

    Implicit None
    Public

    ! Computation timing functions should be private.
    Private :: MarkCT, IncCT, TimeStep, ChooseProcess

    ! -------------------------------------------------------
    ! SOLVER VARIABLES.
    ! -------------------------------------------------------

    ! Simulation time variables.
    Real, Private :: m_t1, m_t2, m_dt, m_tstop, m_dtg, &
                     m_maxts

    ! Computation timing variables.
    Real, Private :: MarkT

    Contains

    ! -------------------------------------------------------

    Subroutine Run(t1, t2, soln, mech, flag)
        ! DESCRIPTION:
        !   Runs Sweep over the given time interval solving the
        !   given solution (ensemble & chemistry) using the
        !   supplied mechanism.

        Use SWPMECH
        Use SWPSOLN
        Use SWPPROCESS
        Use SWPENSEMBLE
        Use SWPCHEM
        Use SWPPARAMS
        Implicit None

        ! ARGUMENTS.
        Real, Intent(INOUT)           :: t1   ! Start and stop times (s).
        Real, Intent(IN)              :: t2
        Type(Solution), Intent(INOUT) :: soln ! Solution to solve.
        Type(Mechanism), Intent(IN)   :: mech ! Mechanism to use.
        Integer, Intent(OUT)          :: flag ! Result flag.

        ! VARIABLES.
        Integer :: N, i, j
        Real    :: chem(soln%Chemistry%Shared%NCHEM), rate_parts(mech%ProcessCount), &
                    jrate, srate
        Real    :: dt_sint, rate, dnp, np, vol !J.Y. Xing
        Character(LEN=100) :: errstr

        ! EXECUTABLE CODE.
        flag = 0

        ! Reset process counters.
        If (RESET_COUNTERS) Then
            Call ResetCounters(soln)
        End If

        ! Copy simulation variables to module level.
        m_t1    = t1
        m_t2    = t2
        m_dt    = 0.0E0
        m_dtg   = (m_t2 - m_t1) / Real(soln%MinStepCount)
        m_maxts = m_dtg / 3.0E0
        m_tstop = 0.0E0
        Do While (m_t1 < m_t2)

            ! ***********************************
            ! CALCULATE SPLITTING STEP LENGTH.
            ! ***********************************

            If (USE_CT) Call MarkCT()

            ! This section works out how long the current time step
            ! should be and puts the end time of the step into
            ! m_tstop.
            If (mech%AnyDeferred) Then
                ! Number of stochastic particles.
                N = ParticleCount(soln%Ensemble)

                ! Get details of the chemical environment.
                chem = GetChem(soln%Chemistry, m_t1, m_t1, flag)
                If (flag < 0) Exit

                ! Now work out the rates of all the different processes.
                rate_parts = RateTerms(m_t1, chem, soln, mech, .True., flag)
                If (flag < 0) Exit

                ! The total rate of the stochastic processes
                jrate = JumpRate(rate_parts, mech)
                ! The total rate of the processes to be split if they
                ! were handled stochastically
                srate = DeferRate(rate_parts, mech)

                ! Work out time to split over.
                IF (n < 20) THEN
                    ! Main process should be inception.
                    m_dt = DataSpacing(soln%Chemistry, m_t1)
                    If (m_dt <= 0.0E0) m_dt = m_dtg
                ELSE
                    ! Splitting time step is expected time for
                    ! split_ratio unsplit events per particle.
                    m_dt = n * soln%SplitRatio / (jrate + 1.0E0)
                END IF

                ! Try to catch extremely long values
                ! caused by srate being very low.
                m_dt    = Min(m_dt, m_dtg)
                m_tstop = Min(m_t1 + m_dt, m_t2)
                m_maxts = m_dt / 3.0E0
            ELSE
                ! No splitting.
                m_tstop = m_t1 + m_dtg
            END IF

            If (USE_CT) Call IncCT(soln%PreSplitCT)

            ! ***********************************
            ! ADVANCE STOCHASTICALLY OVER STEP.
            ! ***********************************

            dt_sint = m_t1
            ! Take the stochastic steps up to step_end_time
            Do While (m_t1 < m_tstop)
                m_dt = TimeStep(soln, mech, flag)
                If (flag < 0) Exit

                ! dt is a module variable that is set by
                ! TimeStep() to the length of the step
                ! just taken.
                m_t1 = m_t1 + m_dt

                ! Keep count of the number of steps.
                soln%StepCount = soln%StepCount + 1
            End Do

            If (flag < 0) Then
                Exit
            Else
                If (m_t1 < m_tstop) Then
                    ! Stochastic steps failed.
                    flag = STEP_END_FAIL_ERR
                    Exit
                End If
            End If
            dt_sint = m_t1 - dt_sint

            ! ***********************************
            ! PERFORM DEFERRED EVENTS.
            ! ***********************************
 
            If (USE_CT) Call MarkCT()

            ! Perform remaining events form deferred processes now.
            If (mech%AnyDeferred) Then
                Call SurfaceIntegrate(0.0E0, m_t1, soln, mech, flag)
                If (flag < 0) Exit
                soln%LPDACount = soln%LPDACount + 1
            End If

            !J.Y. Xing
            If (SINT_ON) Then
                N = ParticleCount(soln%Ensemble)
                Do i = 1, N
                    If (mech%ParticleModel == FRACTAL_MODEL) Then
                        rate = ParticleSinteringRate(soln%Ensemble%Particles(i), mech, chem(soln%Chemistry%Shared%iT))
                        vol = soln%Ensemble%Particles(i)%Properties(iV)
                        np = Max(1.0d0, soln%Ensemble%Particles(i)%Surface**3 / (36*PI*vol**2))
                        !If(np < 2.0d0) Then
                        !    dnp = 3*(np - np**(2.0d0/3.0d0))
                        !Else 
                        !    dnp = 3*(2.0d0 - 2.0d0**(2.0d0/3.0d0))*(np-1.0d0)
                        !End If
                        !np = Max(1.0d0, np - rate*dnp*m_dtg)
                        If(np < 2.0d0) Then
                            np = Max(1.d0, np*exp(-3.d0*(1.d0 - np**(-1.d0/3.d0)) * rate*dt_sint))
                        Else
                            np = Max(1.d0, (np-1.d0)*exp(-3.d0*(2.d0 - 2.d0**(2.d0/3.d0)) * rate*dt_sint) + 1.d0);
                        End If
                        soln%Ensemble%Particles(i)%Surface = (36*PI*np*vol**2)**(1.d0/3.d0)
                        Call CalcProperties(soln%Ensemble%Particles(i), mech)
                    End if
                End Do
            End If

            If (USE_CT) Call IncCT(soln%DeferCT)
        End Do

        If ((flag == 0) .And. (m_t1 < m_t2)) Then
            ! Something must have gone wrong to 
            ! stop us reaching the end.
            flag = RUN_FAIL_ERR
        End If

        ! Set time variables from module level.
        If (flag /= 0) Then
            t1 = m_t1
            If (MSGS) Then
                Print *, "Sweep:\> Fail with error: ", flag
                Call TranslateErr(flag, errstr)
                Print *, "Sweep:\> " // Trim(errstr)
            End If
        Else
            t1 = m_t1
        End If
    End Subroutine

    ! -------------------------------------------------------

    Real Function TimeStep(soln, mech, flag)
        ! DESCRIPTION:
        !   Performs one time step and returns the change
        !   in time (delta t).
        !
        ! RETURNS:
        !   Change in time for the step (s).

        Use SWPMECH
        Use SWPSOLN
        Use SWPCHEM
        Use SWPRNG
        Use SWPPROCESS
        Implicit None

        ! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln ! Solution to solve.
        Type(Mechanism), Intent(IN)   :: mech ! Mechanism to use.
        Integer, Intent(OUT)          :: flag ! Error flag.

        ! VARIABLES.
        Integer :: react_num
        Real    :: chem(soln%Chemistry%Shared%NCHEM), rate_parts(mech%ProcessCount), jrate

        ! EXECUTABLE CODE.
        If (USE_CT) Call MarkCT()

        ! Get details of the chemical environment.
        chem = GetChem(soln%Chemistry, m_t1, m_t1, flag)
        If (flag < 0) Return

        If (USE_CT) Call IncCT(soln%ChemCT)
        If (USE_CT) Call MarkCT()

        !Now work out the rates of all the different processes
        rate_parts = RateTerms(m_t1, chem, soln, mech, .False., flag)
        If (flag < 0) Return

        jrate = JumpRate(rate_parts, mech) 

        IF (jrate > 0) THEN
           TimeStep = - Log(Rnd()) 
           TimeStep = TimeStep / jrate
        ELSE
           ! Do not divide by zero.
           TimeStep = 1.0D+30 ! almost +inf
        END IF

        If (USE_CT) Call IncCT(soln%RatesCT)
        If (USE_CT) Call MarkCT()

        IF (TimeStep > m_maxts) THEN
           ! Exponential rvs are memoryless which means we can
           ! truncate the first waiting time and add on another
           ! with the same rate.  A slight twist here is using
           ! new rates which reflect the change of the conditions
           ! over time.
           TimeStep = m_maxts
           react_num = 0
        ELSE
           ! Ordinary event
           react_num = ChooseProcess(rate_parts, mech%ProcessCount)
        END IF

        If (react_num > 0) Then
           Call DoProcess(react_num, m_t1, chem, soln, mech, flag)
        End If

        If (USE_CT) Call IncCT(soln%ProcessCT)
    End Function

    ! -------------------------------------------------------

    Integer Function ChooseProcess(rates, N)
        ! DESCRIPTION:
        !	Selects a process according to the rates of each and returns
        !	the index of the selected process
        ! RETURNS:
        !	Index of process to perform.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)     :: N
        Real, Intent(IN)        :: rates(N)

        ! VARIABLES.
	Real :: rand

        ! EXECUTABLE CODE.

        ! Generate a uniform random number.
        rand = Rnd() * Sum(rates)

        ! Select a process using DIV.
        ChooseProcess = 1
        rand = rand - rates(1)
        Do While (rand > 0 .AND. ChooseProcess < N)
            ChooseProcess = ChooseProcess + 1
            rand = rand - rates(ChooseProcess)
        End Do
    End Function

    ! -------------------------------------------------------
    ! COMPUTATION TIMING ROUTINES.
    !
    !   These routines are used to profile the coupled
    !   chemistry-soot solver.
    !
    ! -------------------------------------------------------

    Subroutine MarkCT()
        ! DESCRIPTION:
        !	Saves the current time.  Subsequent calls to
        !	"increment time" functions will use MarkT as the
        !	start time.
        Call CPU_Time(MarkT)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine IncCT(tim) 
        ! DESCRIPTION:
        !	Increments a computation time variable using the
        !	MarkT as the start time and the current time
        !	as the end time of the interval.
        Real, Intent(INOUT) :: tim
        Real                :: now
        Call CPU_Time(now)
        tim = tim + (now - MarkT)
    End Subroutine

End Module

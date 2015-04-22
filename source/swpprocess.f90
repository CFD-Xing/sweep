! *****************************************************************************
!
! File:					swpprocess.f90
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
!	Particle dynamics controller for Sweep.  Defines the particle-particle
!	dynamics of coagulation and the surface chemistry of single soot particles.
!	Maintains the soot mechanism currently loaded into Sweep and provides
!	functions for calculating the rates of all processes.
!
!	All inception and surface processes are soft-coded.  Only coagulation is
!	hard-coded into Sweep.  This design decision was taken because the coagulation
!	kernel will change rarely, but the surface reactions may frequently change.
!
! Functions:
!	-- (Rate(s) calculation) -------------------------------------------------
!	RateTerms			-	Calculates the rates of all soot processes and 
!							returns them as an array.
!	GroupedRateTerms	-	Returns rates terms grouped together by chemical
!							reaction.
!	ChemRateTerms		-	Calculates the chemistry-dependent terms for inceptions
!							and surface processes.
!	InceptionRates		-	Calculates the rates of all inception processes.
!	SurfaceRates		-	Calculates the rates of all surface processes.
!	ChemSurfaceRates	-	Calculates the surface process rate terms which only
!							depend on the gas-phase chemistry.
!	SootSurfaceRates	-	Calculates the surface process rate terms which depend
!							on the soot (and the active site fraction).
!	SingleParticleRate	-	Calculates the rate of a given reaction for a given
!							particle.
!	VolChangeRate		-	Calculates the rate of change of a particle's
!							volume due to all surface processes.
!	ParticleDeferRate	-	Calculates the combined rate of all deferred processes
!							for a single particle.
!	JumpRate			-	Calculates the combined rate of all non-deferred
!							processes.
!	DeferRate			-	Calculates the combined rate of all deferred
!							processes.
!	MomentsRates    	-	Calculates the rate of change of the first six moments.
!   ChemChangeRates     -   Calculates the rate of change of chemical species due
!                           to soot processes.
!	-- (Mech Info) ------------------------------------------------------------
!	DeferredProcessCount	-	Returns the number of deferred processes in the
!								mechanism.
!	NonDeferredProcessCount	-	Returns the number of non-deferred processes in
!								the mechanism.
!	-- (DSA) ------------------------------------------------------------------
!	DoProcess			-	Performs a process given by an id number.
!	Incept				-	Performs an inception reaction.
!	Surface				-	Performs a given surface reaction.
!	DSACoagulate		-	Performs a coagulation of two particles.
!	Ficticious			-	Defines a test to see if a reaction is ficticious or
!							not.
!   Inflow              -   Performs an inflow event.
!   Outflow             -   Performs an outflow event.
!	-- (LPDA) -----------------------------------------------------------------
!	UpdateDeferred		-	Performs all deferred events on a stochastic particle.
!	SurfaceIntegrate	-	Updates entire particle ensemble by performing
!							deferred events on all particles.
!	DeferredSurface		-	Performs a deferred surface event on a given particle.
!	-- (Chemistry Interaction) ------------------------------------------------
!	UpdateChemistry(S/I) -	Updates the gas-phase chemistry according to the
!					 	given reaction.
!	-- (Auxiliary routines) ---------------------------------------------------
!	SampleVolume		-	Calculates the sample volume of the particle ensemble
!							for the given temperature.
!	SetScaling			-	Sets the scaling parameters for Sweep.  These are
!							the sample volume and the associated gas-phase
!							temperature for the initial solution point.
!	SetBaseTemperature	-	Sets the base temperature used to calculate the
!							sample volume by Boyle's law.
! *****************************************************************************

Module SWPPROCESS
	! -------------------------------------------------------
	! IMPORT REQUIRED MODULES.
	! -------------------------------------------------------
	Use SWPERR ! Sweep error codes.
	! -------------------------------------------------------

	Implicit None
	Public

    Interface UpdateChemistry
        Module Procedure UpdateChemistryI
        Module Procedure UpdateChemistryS
        Module Procedure UpdateChemistrySV
    End Interface

	Contains

	! -------------------------------------------------------
	! RATE CALCULATION FUNCTIONS.
	!
	!	These routines calculate the rates of different
	!	soot processes.
	!
	! -------------------------------------------------------

	Function RateTerms(t, chem, soln, mech, inclDeferred, flag)
		! DESCRIPTION:
		!	Calculates the rates of all processes and returns
		!	them in an ordered array.
		! RETURNS:
		!	The rates of all processes.

		Use SWPCHEM
		Use SWPENSEMBLE
        Use SWPMECH
        Use SWPSOLN
        Use SWPPARAMS
        Use SWPCOAG
		Implicit None

		! ARGUMENTS.
        Type(Solution), Intent(IN)  :: soln  ! Current solution.
        Type(Mechanism), Intent(IN) :: mech  ! Mechanism to use for rate calculation.
		Real, Intent(IN)     :: t, chem(:)   ! Flow time and chemical conditions.
		Logical, Intent(IN)  :: inclDeferred ! Should deferred process rates also be calculated.
		Integer, Intent(OUT) :: flag         ! Error flag.
		Real ::	RateTerms(mech%ProcessCount) ! Return value.

		! VARIABLES.
		Integer	::	err, N
		Real	::	sums(PROPERTY_COUNT), vol

		! EXECUTABLE CODE.
		flag = 0
		RateTerms = 0.0E0

		! Get required system properties.
		N = ParticleCount(soln%Ensemble)
		sums = GetParticleSums(soln%Ensemble)
        vol  = SampleVolume(soln, chem(soln%Chemistry%Shared%iT))
		
        ! Calculate inception rates (note: inception cannot be switched off).
		RateTerms(1:mech%InceptionCount) = InceptionRates(chem, vol, soln, mech, err)

        ! Calculate surface process rates.
        If (SURF_ON) Then
	    RateTerms(mech%ISR:mech%ISR+mech%ReactionCount-1) = SurfaceRates(chem, N, sums, soln, mech, err)

	    ! If any processes are deferred then the surface rates
	    ! need to be scaled by a majorant factor.
	    If (mech%AnyDeferred) Then
		    RateTerms(mech%ISR:mech%ISR+mech%ReactionCount-1) = &
                RateTerms(mech%ISR:mech%ISR+mech%ReactionCount-1) * SURF_MAJ
	    End If
        End If

		! Calculate coagulation rates.
		If (COAG_ON) Then
			RateTerms(mech%ICG:mech%ICG+NCOAGTERMS-1) = &
            CoagRateTerms(chem(soln%Chemistry%Shared%iT), chem(soln%Chemistry%Shared%iP), &
                          N, sums, vol, flag)

			If (flag < 0) Then
    			! Failed to calculate coagulation terms.
				RateTerms(mech%ICG:mech%ICG+NCOAGTERMS-1) = 0.0E0
			End If
		End If

		! Inflow and outflow rate.
        If (soln%SolveInOut) Then
            RateTerms(mech%IIN)  = 0.0E0
	    RateTerms(mech%IOUT) = N / soln%ResidenceTime
        End If

		If (Any(RateTerms < 0.0E0)) Then
			! Something has gone wrong to give us a
			! negative rate term.
			flag = NEG_RATE_ERR
			RateTerms = 0.0E0
		End If

        ! Clear deferred processes if not required.
		If (.Not. inclDeferred) Where(mech%DeferMask) RateTerms = 0.0E0
	End Function

	! -------------------------------------------------------

	Function GroupedRateTerms(t, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Calculates the rates of all processes, grouping
		!	them together by physical process.  Note that Sweep
        !   splits up processes (e.g. coagulation) into multiple
        !   rate terms.
		! RETURNS:
		!	The grouped rates of all processes.

		Use SWPCHEM
		Use SWPENSEMBLE
        Use SWPSOLN
        Use SWPMECH
        Use SWPPARAMS, only: SURF_MAJ
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN) :: t, chem(:)     ! Current flow time and chemical conditions.
        Type(Solution), Intent(IN)  :: soln    ! Current solution.
        Type(Mechanism), Intent(IN) :: mech    ! Mechanism to use for rate calculation.
		Integer, Intent(OUT) :: flag           ! Error flag.
		Real ::	GroupedRateTerms(mech%ProcessCount) ! Return value.

		! VARIABLES.
		Integer	::	i, j, ir
		Real	::	rates(mech%ProcessCount)

		! EXECUTABLE CODE.
		flag = 0
		ir   = 0
		GroupedRateTerms = 0.0E0
		
		! Get the ungrouped rate terms.
		rates = RateTerms(t, chem, soln, mech, .True., flag)

        If (mech%AnyDeferred) Then
            rates(mech%ISR:mech%ISR+mech%ReactionCount) = &
            rates(mech%ISR:mech%ISR+mech%ReactionCount) / SURF_MAJ
        End If

        If (flag == 0) Then
	    ! Group together the rate terms.
	    Do i = 1, mech%GroupCount
		    Do j = 1,  mech%Groups(i)
			    ir = ir + 1
			    GroupedRateTerms(i) = GroupedRateTerms(i) + rates(ir)
		    End Do
	    End Do
        End If

        GroupedRateTerms = GroupedRateTerms / SampleVolume(soln, chem(soln%Chemistry%Shared%iT))
	End Function

	! -------------------------------------------------------

	Function InceptionRates(chem, vol, soln, mech, flag)
		! DESCRIPTION:
		!	Calculates the rates of all inception processes.
		! RETURNS:
		!	The inception rates.

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
        Use SWPCOAG_MODEL
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN) :: chem(:)     ! Current chemical conditions.
		Real, Intent(IN) :: vol         ! Sample volume of the ensemble.
        Type(Solution), Intent(IN) :: soln  ! Solution to use.
        Type(Mechanism), Intent(IN) :: mech ! Mechanism to use for rate calculation.
		Integer, Intent(OUT) :: flag        ! Error flag.
		Real :: InceptionRates(mech%InceptionCount) ! Return values.

		! VARIABLES.
		Integer	::	i, iSp
		Real	::	fm, sf, sqrtT, T_mu, T_P

		! EXECUTABLE CODE.
		flag = 0

		! Precalculate useful quantities.
		sqrtT = Sqrt(chem(soln%Chemistry%Shared%iT))
		T_mu  = chem(soln%Chemistry%Shared%iT) / Visc(chem(soln%Chemistry%Shared%iT))
		T_P  = chem(soln%Chemistry%Shared%iT) / chem(soln%Chemistry%Shared%iP)

		Do i = 1, mech%InceptionCount
			! Free-molecular and slip-flow kernels.
			fm = sqrtT * mech%Inceptions(i)%kfm
			sf = T_mu * (mech%Inceptions(i)%ksf1 + (T_P * mech%Inceptions(i)%ksf2))

			! Kernel is harmonic mean of free-molecular and slip-flow.  Remember
			! to scale the inception rates to give the number of inception events
            ! in the sample volume, instead of per unit volume.
!			InceptionRates(i) = mech%Inceptions(i)%A * fm * vol
			InceptionRates(i) = mech%Inceptions(i)%A * HalfHarmonicMean(fm, sf) * vol
!   			InceptionRates(i) = 2.0E12 * sqrtT * vol * NA

			! Chemical species terms.
			Do iSp = 1, mech%Inceptions(i)%NREAC
				InceptionRates(i) = InceptionRates(i) * &
									((NA * chem(mech%Inceptions(i)%Species(iSp))) ** mech%Inceptions(i)%Stoich(iSp))
			End Do
		End Do
	End Function

	! -------------------------------------------------------

	Function SurfaceRates(chem, N, sums, soln, mech, flag)
		! DESCRIPTION:
		!	Calculates the rates of all surface processes.
		! RETURNS:
		!	The surface rates.

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
		Use SWPENSEMBLE
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN) ::	chem(:)              ! Current chemical conditions.
        Integer, Intent(IN)  :: N                    ! Particle count.
		Real, Intent(IN) ::	sums(PROPERTY_COUNT) ! Sums of particle properties.
        Type(Solution), Intent(IN) :: soln           ! Solution to use.
        Type(Mechanism), Intent(IN) :: mech          ! Mechanism to use for rate calculation.
		Integer, Intent(OUT) ::	flag                 ! Error flag.
		Real ::	SurfaceRates(mech%ReactionCount)     ! Return values.

		! VARIABLES.
		Integer	:: i, iSp
		Real    :: invRT, lnT

		! EXECUTABLE CODE.
		flag = 0
        SurfaceRates = 0.0E0

		If (N > 0) Then
			! Precalculate useful values.
			invRT = 1.0E0 / (R * chem(soln%Chemistry%Shared%iT))
            lnT   = Log(chem(soln%Chemistry%Shared%iT))

			Do i = 1, mech%ReactionCount
				! Temperature dependent rate constant.
				SurfaceRates(i) = mech%Reactions(i)%A

				! Chemical species terms.
				Do iSp = 1, mech%Reactions(i)%NREAC
					SurfaceRates(i) = SurfaceRates(i) * &
								  (chem(mech%Reactions(i)%Species(iSp)) ** mech%Reactions(i)%Stoich(iSp))
				End Do

				! Temperature power and exponential dependency.
				SurfaceRates(i) = &
                SurfaceRates(i) * Exp((mech%Reactions(i)%n * lnT) - (mech%Reactions(i)%E * invRT))

                If (mech%UseHACA) Then
			    If (mech%Reactions(i)%ID == iAS) Then
				    ! Active site fraction (if part of the chemistry profile).  If 
				    ! the active surface model gets alpha per particle (ACT_SURF_PARTICLE)
				    ! then there is no need to calculate alpha as it has
				    ! already be applied to each particle.
				    If (mech%ActSurfModel == ACT_SURF_CONST) Then
					    SurfaceRates(i) = SurfaceRates(i) * mech%ConstActSurf * chem(soln%Chemistry%Shared%iRad)
				    ElseIf (mech%ActSurfModel == ACT_SURF_ABF) Then
					    SurfaceRates(i) = SurfaceRates(i) * &
                                            ABFAlpha(chem(soln%Chemistry%Shared%iT), sums(iSz)/N) * chem(soln%Chemistry%Shared%iRad)
				    ElseIf (mech%ActSurfModel == ACT_SURF_PROFILE) Then
					    SurfaceRates(i) = SurfaceRates(i) * chem(soln%Chemistry%Shared%iAlpha) * chem(soln%Chemistry%Shared%iRad)
                        ElseIf (mech%ActSurfModel == ACT_SURF_PARTICLE) Then
					    SurfaceRates(i) = SurfaceRates(i) * chem(soln%Chemistry%Shared%iRad)
				    End If
			    End If
                End If

				! Soot property terms.
				If (mech%Reactions(i)%ID > 0) Then
					SurfaceRates(i) = SurfaceRates(i) * sums(mech%Reactions(i)%ID)
				ElseIf (mech%Reactions(i)%ID == UNIFORM_ID) Then
					SurfaceRates(i) = SurfaceRates(i) * N
!                    SurfaceRates(i) = SurfaceRates(i) * TypedParticleCount(soln%Ensemble, mech%Reactions(i)%PTypeIn)
				ElseIf (mech%Reactions(i)%ID == CUSTOM_ID) Then
					SurfaceRates(i) = &
                    SurfaceRates(i) * CustomParticleSum(soln%Ensemble, &
                                                        mech%Reactions(i)%P(1),mech%Reactions(i)%P(2),&
                                                        mech%Reactions(i)%P(3),mech%Reactions(i)%P(4))
				End If
			End Do
		Else
			SurfaceRates = 0.0E0
		End If
	End Function

	! -------------------------------------------------------

	Real Function SingleParticleRate(rxn, sp, chem, soln, mech, maj, flag)
		! DESCRIPTION:
		!	Calculates the rate of the given reaction for
		!	a single particle using the given chemical
		!	conditions.
		! RETURNS:
		!	The reaction rate for a single particle.

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
		Use SWPPART
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Type(Reaction), Intent(IN) :: rxn      ! Index of the reaction for which 
                                               ! to calculate rate.
		Type(StochParticle), Intent(IN)	:: sp  ! Stochastic particle for which
                                               ! to calculate the rate.
		Real, Intent(IN) :: chem(:)        ! Current chemical conditions.
        Type(Solution), Intent(IN) :: soln     ! Solution to use.
        Type(Mechanism), Intent(IN) :: mech    ! Mechanism to use.
		Logical, Intent(IN) :: maj           ! Flag indicating whether or not to calculate
	                                       ! the majorant rate.
		Integer, Intent(OUT) :: flag           ! Error flag.

		! VARIABLES.
		Integer	:: i, err

		! EXECUTABLE CODE.
		flag = 0

		If ((rxn%A > 0.0E0) .And. (rxn%PTypeIn == sp%TypeID)) Then
			! Rate consant.
			SingleParticleRate = rxn%A

			! Chemical species terms.
			Do i = 1, rxn%NREAC
				SingleParticleRate = SingleParticleRate * (chem(rxn%Species(i)) ** rxn%Stoich(i))
			End Do

			! Temperature dependent terms.
			SingleParticleRate = SingleParticleRate * (chem(soln%Chemistry%Shared%iT)**rxn%n)
			SingleParticleRate = SingleParticleRate * Exp(-rxn%E / (R * chem(soln%Chemistry%Shared%iT)))

            If (mech%UseHACA) Then
		    If (rxn%ID == iAS) Then
			    ! Active site fraction (if part of the chemistry profile).
			    If (mech%ActSurfModel == ACT_SURF_CONST) Then
				    SingleParticleRate = SingleParticleRate * mech%ConstActSurf * chem(soln%Chemistry%Shared%iRad)
			    ElseIf (mech%ActSurfModel == ACT_SURF_PROFILE) Then
				    SingleParticleRate = SingleParticleRate * chem(soln%Chemistry%Shared%iAlpha) * chem(soln%Chemistry%Shared%iRad)
			    ElseIf (mech%ActSurfModel == ACT_SURF_ABF) Then
				    SingleParticleRate = SingleParticleRate * &
                                    ABFAlpha(chem(soln%Chemistry%Shared%iT), Real(SizeP(sp))) * chem(soln%Chemistry%Shared%iRad)
			    ElseIf (mech%ActSurfModel == ACT_SURF_PARTICLE) Then
				    SingleParticleRate = SingleParticleRate * chem(soln%Chemistry%Shared%iRad)
			    End If
		    End If
            End If

			! Soot terms.
			If (rxn%ID > 0) Then
				SingleParticleRate = SingleParticleRate * sp%Properties(rxn%ID)
			ElseIf (rxn%ID == CUSTOM_ID) Then
				SingleParticleRate = &
                SingleParticleRate * CalcParticleWeight(sp, rxn%P(1), rxn%P(2), rxn%P(3), rxn%P(4))
			End If

			! If required, scale the rate using a majorant factor.
			If (maj) SingleParticleRate = SingleParticleRate * SURF_MAJ
		Else
			SingleParticleRate = 0.0E0
		End If
   	End Function

	! -------------------------------------------------------

	Real Function VolChangeRate(sp, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Calculates the volume change rate for
		!	a single particle.
		!
		! RETURNS:
		!	The volume change rate for the given particle.

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
		Use SWPPART
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(IN)	:: sp   ! Particle for which to calculate the rate.
		Real, Intent(IN)				:: chem(:) ! Current chemical conditions.
        Type(Solution), Intent(IN)      :: soln    ! Solution to use.
        Type(Mechanism), Intent(IN)     :: mech    ! Mechanism used for rate calculation.
		Integer, Intent(OUT)			:: flag   ! Error flag.

		! VARIABLES.
		Integer :: i

		! EXECUTABLE CODE.
		flag = 0
		VolChangeRate = 0.0E0

		If (SURF_ON) Then
			Do i = 1, mech%ReactionCount
				If (mech%DeferMask(mech%ISR+i-1)) Then
					VolChangeRate = &
                    VolChangeRate + SingleParticleRate(mech%Reactions(i), sp, chem, soln, mech, .False., flag) * &
							    Volume(mech%Reactions(i)%dComp, mech%Components, mech%ComponentCount)

					! Failed to calculate the single-particle
					! rate for a surface reaction.
					If (flag < 0) Exit
				End If
			End Do
		End If
	End Function

	! -------------------------------------------------------

	Real Function ParticleDeferRate(sp, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Calculates the total rate of all deferred
		!	processes for a single particle.
		!
		! RETURNS:
		!	Total deferred rate for given particle..

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
		Use SWPPART
        USe SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(IN)	:: sp      ! Particle for which to calculate the rate.
		Real, Intent(IN)            :: chem(:) ! Current chemical conditions.
        Type(Solution), Intent(IN)      :: soln    ! Solution to use.
        Type(Mechanism), Intent(IN)     :: mech    ! Mechanism used for rate calculation.
		Integer, Intent(OUT)            :: flag   ! Error flag.

		! VARIABLES.
		Integer :: i

		! EXECUTABLE CODE.
		flag = 0
		ParticleDeferRate = 0.0E0

		If (SURF_ON) Then
			Do i = 1, mech%ReactionCount
				If (mech%DeferMask(mech%ISR+i-1)) Then
					ParticleDeferRate = &
                    ParticleDeferRate * SingleParticleRate(mech%Reactions(i), sp, chem, soln, mech, .False., flag)

					! Failed to calculate the single-particle
					! rate for a surface reaction.
					If (flag < 0) Exit
				End If
			End Do
		End If
	End Function

	! -------------------------------------------------------

	Real Function ParticleSinteringRate(sp, mech, T)
		! DESCRIPTION:
		!	Calculates the sintering rate for a single particle.
		! RETURNS:
		!	Sintering rate for given particle.

        Use SWPMECH
		Use SWPPART
        USe SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(IN)	:: sp   ! Particle for which to calculate the rate.
        Type(Mechanism), Intent(IN)     :: mech ! Mechanism used for rate calculation.
		Real, Intent(IN)            :: T    ! Current temperature (K).

		! VARIABLES.
        Real :: dS, So, tau

		! EXECUTABLE CODE.
		If (SINT_ON) Then
            So = EquivSphereSurface(sp, mech%Components, mech%ComponentCount)
            dS = SINT_SCALE * So

            If (sp%Surface >= So+dS) Then
                tau = mech%SintParams(1) * T * (6.0E0 * sp%Properties(iV) / sp%Surface)**4 * exp(mech%SintParams(2) / T)
                ParticleSinteringRate = (sp%Surface - So) / (dS * tau)
            Else
    	    ParticleSinteringRate = 0.0E0
            End If
        Else
    		ParticleSinteringRate = 0.0E0
		End If
	End Function

	! -------------------------------------------------------

	Subroutine MomentsRates(rates, N, sums, vol, soln, mech, mrates)
		! DESCRIPTION:
		!	Calculates the rates of change of the first 6
        !   moments given a rate terms array.

        Use SWPMECH
        Use SWPSOLN
		Use SWPPART
        Use SWPCOAG
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)              :: rates(:)    ! Process rates.
        Real, Intent(IN)              :: sums(:)     ! Particle property sums.
        Real, Intent(IN)              :: vol         ! Sample volume associated with rates.
        Integer, Intent(IN)           :: N           ! Particle count.
        Type(Solution), Intent(IN)    :: soln        ! Solution to use.
        Type(Mechanism), Intent(IN)   :: mech        ! Mechanism used for rate calculation.
        Double Precision, Intent(OUT) :: mrates(0:5) ! Return moments rates.

		! EXECUTABLE CODE.
        Integer :: i
        Double Precision :: mean, truerates(mech%ProcessCount), dx

        truerates = Dble(rates) ! Convert rates to double.

		! dM0/dt = Rincept - Rcoag + Rinflow - Routflow
		mrates(0) = Sum(truerates(1:mech%InceptionCount)) - Sum(truerates(mech%ICG:mech%ICG+NCOAGTERMS-1))
        If (soln%SolveInOut) Then
        mrates(0) = mrates(0) + truerates(mech%IIN)
            mrates(0) = mrates(0) - truerates(mech%IOUT)
        End If
        
		! dMn/dt =  Sum(Rincept * Xincept^n) + Sum(Rsurf  * Xchange^n) + 
        !          (Rinflow * Xmean_inf^n) - (Routflow * Xmean^n)
        mrates(1) = 0.0E0
		Do i = 1, mech%InceptionCount
            dx = Dble(SizeP(mech%Inceptions(i)%dComp))
			mrates(1) = mrates(1) + (truerates(i) * dx)
			mrates(2) = mrates(2) + (truerates(i) * Sign(dx**2, dx))
			mrates(3) = mrates(3) + (truerates(i) * Sign(dx**3, dx))
			mrates(4) = mrates(4) + (truerates(i) * Sign(dx**4, dx))
			mrates(5) = mrates(5) + (truerates(i) * Sign(dx**5, dx))
		End Do
        If (N > 0) Then
            ! Surface reactions.
	    Do i = 1, mech%ReactionCount
                dx = Dble(SizeP(mech%Reactions(i)%dComp))
		    mrates(1) = mrates(1) + (truerates(mech%ISR+i-1) * dx)
		    mrates(2) = mrates(2) + (truerates(mech%ISR+i-1) * Sign(dx**2, dx))
		    mrates(3) = mrates(3) + (truerates(mech%ISR+i-1) * Sign(dx**3, dx))
		    mrates(4) = mrates(4) + (truerates(mech%ISR+i-1) * Sign(dx**4, dx))
		    mrates(5) = mrates(5) + (truerates(mech%ISR+i-1) * Sign(dx**5, dx))
	    End Do
            ! Outflow.
            If (soln%SolveInOut) Then
                mean = Dble(sums(iV)) / Dble(N)
	        mrates(1) = mrates(1) - (truerates(mech%IOUT) * mean)
	        mrates(2) = mrates(2) - (truerates(mech%IOUT) * mean**2)
	        mrates(3) = mrates(3) - (truerates(mech%IOUT) * mean**3)
	        mrates(4) = mrates(4) - (truerates(mech%IOUT) * mean**4)
                mrates(5) = mrates(5) - (truerates(mech%IOUT) * mean**5)
            End If
        End If

        ! Scale to a unit volume.
        mrates = mrates / Dble(vol)
	End Subroutine


	! -------------------------------------------------------

	Function ChemChangeRates(chem, soln, mech) Result(rates)
		! DESCRIPTION:
		!	Calculates the rate of change of chemical species
        !   due to soot processes based on the current soot
        !   population.
		! RETURNS:
		!	Rates of change of chemical species due to soot
        !   processes.  The rates are in units of moles/smp vol.

        Use SWPMECH
        Use SWPSOLN
        Use SWPENSEMBLE
        Use SWPCHEM
		Implicit None
		
        ! ARGUMENTS.
        Real, Intent(IN)            :: chem(:)  ! Chemistry to use
        Type(Solution), Intent(IN)  :: soln     ! Particle solution.
        Type(Mechanism), Intent(IN) :: mech     ! Mechanism used for rate calculation.
        Real :: rates(soln%Chemistry%Shared%NCHEMEXT) ! Chemistry change rates.

        ! VARIABLES.
        Integer :: err, i, isp
        Real    :: sums(PROPERTY_COUNT), pr(mech%InceptionCount+mech%ReactionCount)
        Real    :: vol

        ! EXECUTABLE CODE.
        rates = 0.0E0

        ! Get prerequisite values.
        sums = GetParticleSums(soln%Ensemble)
        vol  = SampleVolume(soln, chem(soln%Chemistry%Shared%iT))

        ! Get the inception and surface rates.
        pr(1:mech%InceptionCount)  = InceptionRates(chem, vol, soln, mech, err)
        If (err /= 0) Return
        pr(mech%InceptionCount+1:) = SurfaceRates(chem, ParticleCount(soln%Ensemble), sums, soln, mech, err)
        If (err /= 0) Return

        ! Calculate rates of change.
        Do i = 1, mech%InceptionCount
            Do isp = 1, mech%Inceptions(i)%NREAC
                rates(mech%Inceptions(i)%Species(isp)) = &
                rates(mech%Inceptions(i)%Species(isp)) - (mech%Inceptions(i)%Stoich(isp) * pr(i))
            End Do
            Do isp = mech%Inceptions(i)%NREAC + 1, mech%Inceptions(i)%NREAC + mech%Inceptions(i)%NPROD
                rates(mech%Inceptions(i)%Species(isp)) = &
                rates(mech%Inceptions(i)%Species(isp)) + (mech%Inceptions(i)%Stoich(isp) * pr(i))
            End Do
        End Do
        Do i = 1, mech%ReactionCount
            Do isp = 1, mech%Reactions(i)%NREAC
                rates(mech%Reactions(i)%Species(isp)) = &
                rates(mech%Reactions(i)%Species(isp)) - (mech%Reactions(i)%Stoich(isp) * &
                                                         pr(i+mech%InceptionCount))
            End Do
            Do isp = mech%Reactions(i)%NREAC + 1, mech%Reactions(i)%NREAC + mech%Reactions(i)%NPROD
                rates(mech%Reactions(i)%Species(isp)) = &
                rates(mech%Reactions(i)%Species(isp)) + (mech%Reactions(i)%Stoich(isp) * &
                                                         pr(i+mech%InceptionCount))
            End Do
        End Do

        ! Converts to mol/cm3s.
        rates = rates / (NA * vol)
	End Function

	! -------------------------------------------------------

	Real Function JumpRate(rates, mech)
		! DESCRIPTION:
		!	Calculates combined rate of all non-deferred
		!	processes.
		! RETURNS:
		!	Total non-deferred rate.
        Use SWPMECH
		Implicit None
		Real, Intent(IN)            :: rates(:) ! All process rates.
        Type(Mechanism), Intent(IN) :: mech     ! Mechanism used for rate calculation.
		JumpRate = Sum(rates, .Not. mech%DeferMask)
	End Function

	! -------------------------------------------------------

	Real Function DeferRate(rates, mech)
		! DESCRIPTION:
		!	Calculates combined rate of all deferred
		!	processes.
		! RETURNS:
		!	Total deferred rate.
        Use SWPMECH
		Implicit None
		Real, Intent(IN)            :: rates(:) ! All process rates.
        Type(Mechanism), Intent(IN) :: mech     ! Mechanism used for rate calculation.
		DeferRate = Sum(rates, mech%DeferMask)
	End Function

	! -------------------------------------------------------
	! DIRECT SIMULATION EVENTS.
	!
	!	These routines enact directly different soot
	!	processes, including inception, surface reaction and
	!	coagulation.
	!
	! -------------------------------------------------------

	Subroutine DoProcess(index, t, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Performs a soot process given by the index
		!	argument.

        Use SWPMECH
        Use SWPSOLN
        Use SWPCOAG
		Use SWPCHEM
		Use SWPPART
		Implicit None

		! ARGUMENTS.
		Integer, Intent(IN)          :: index    ! Index of process to perform.
		Real, Intent(IN)          :: t, chem(:) ! Current time and chemical conditions.
        Type(Solution), Intent(INOUT) :: soln       ! Particle solution.
        Type(Mechanism), Intent(IN)   :: mech       ! Mechanism to use.
		Integer, Intent(OUT)          :: flag       ! Error flag.

		! VARIABLES.
		Integer :: i

		! EXECUTABLE CODE.
		flag = 0

		If (index <= mech%InceptionCount) Then
			! Process is an inception.		
			Call Incept(index, t, chem, soln, mech, flag)
		ElseIf (index < mech%ICG) Then
			! Process is a surface reaction.
			Call Surface(index, t, chem, soln, mech, flag)
		ElseIf (index <= mech%ICG+NCOAGTERMS-1) Then
			! Process is a coagulation.
			i = index - mech%ICG
			Select Case (i)
				Case (0)
					! Slip-flow; both particles chosen uniformly.
					Call DSACoagulate(0, 0, t, chem, soln, mech, SLIPMAJ, flag)
				Case (1)
					! Slip-flow; one by d, the other by d^-1.
					Call DSACoagulate(iD, iD_1, t, chem, soln, mech, SLIPMAJ, flag)
				Case (2)
					! Slip-flow; one uniformly, the other by d^-1.
					Call DSACoagulate(0, iD_1, t, chem, soln, mech, SLIPMAJ, flag)
				Case (3)
					! Slip-flow; one by d, the other by d^-2.
					Call DSACoagulate(iD, iD_2, t, chem, soln, mech, SLIPMAJ, flag)
				Case (4)
					! Free-molecular; one uniformly, the other by d^2*m^1/2.
					Call DSACoagulate(0, iD2M_1_2, t, chem, soln, mech, FREEMAJ, flag)
				Case (5)
					! Free-molecular; one by d^2, the other by m^-1/2.
					Call DSACoagulate(iD2, iM_1_2, t, chem, soln, mech, FREEMAJ, flag)
				Case (6)
					! Free-molecular; not currently used.
                    ! In use for HD coagulation term.
!					Call DSACoagulate(0, 0, t, chem, soln, mech, NOMAJ, flag)
				Case Default
					! Unknown process index, alert calling code.
					flag = PROCESS_RANGE_ERR
			End Select
		ElseIf (index == mech%IIN) Then
            ! Process is inflow.
            Print *, "Inflow?"
			Call Inflow(soln, flag)
        ElseIf (index == mech%IOUT) Then
            ! Process is outflow.
			Call Outflow(t, soln, mech, flag)
        Else
            flag = PROCESS_RANGE_ERR
		End If

		If (flag >= 0) Then
			! Count the number of times each process
			! is performed.
			If (flag /= FICTICIOUS_EVENT) Then
				soln%ProcessCounter(index) = soln%ProcessCounter(index) + 1
			Else
				soln%FicticiousCounter(index) = soln%FicticiousCounter(index) + 1
			End If
		End If
	End Subroutine

	! -------------------------------------------------------
	
	Subroutine Incept(i, t, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Performs a particle inception of the given
		!	size.

        Use SWPMECH
        Use SWPSOLN
		Use SWPENSEMBLE
		Use SWPCHEM
        Use SWPPARAMS, only: INITIAL_ACT_SURF
		Implicit None

		! ARGUMENTS.	
		Integer, Intent(IN)        :: i          ! Reaction being done.
		Real, Intent(IN)        :: t, chem(:) ! Current time and chemical conditions.
        Type(Solution), Intent(INOUT) :: soln     ! Particle solution.
        Type(Mechanism), Intent(IN) :: mech       ! Mechanism to use.
		Integer, Intent(OUT)        :: flag	  ! Error flag.

		! VARIABLES.
		Integer				:: err
		Type(StochParticle) :: sp

		! EXECUTABLE CODE.
		flag = 0

		! Create a new particle.
        sp = CreateParticle(mech%Inceptions(i)%PType, mech%Inceptions(i)%dComp, &
                            mech%Inceptions(i)%dS, t, t, INITIAL_ACT_SURF, &
                            mech%Inceptions(i)%dTrack, mech)

		If (AddParticle(soln%Ensemble, sp) >= 0) Then
			! Update chemistry.
            If (CHEM_TYPE == VARIABLE_CHEM) Then
		    Call UpdateChemistry(t, mech%Inceptions(i), 1, soln%Chemistry, SampleVolume(soln, chem(soln%Chemistry%Shared%iT)), flag)
            End If
		Else
			! Failed to add particle to ensemble.
			flag = ADD_PARTICLE_ERR
		End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine Surface(irxn, t, chem, soln, mech, flag)
		! DESCRIPTION:
		!	Performs a surface event on a particle in the
		!	ensemble.

        Use SWPMECH
        Use SWPSOLN
		Use SWPENSEMBLE
		Use SWPCHEM
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Integer, Intent(IN)        :: irxn       ! Reaction being done.
		Real, Intent(IN)        :: t, chem(:) ! Current time and chemical conditions.
        Type(Solution), Intent(INOUT) :: soln     ! Particle solution.
        Type(Mechanism), Intent(IN) :: mech       ! Mechanism to use.
		Integer, Intent(OUT)        :: flag	  ! Error flag.

		! VARIABLES.
		Integer				:: i, err
		Type(StochParticle) :: sp
		Real				:: maj, true
        Type(Reaction)      :: rxn

		! EXECUTABLE CODE.
		flag = 0
        rxn = mech%Reactions(irxn-mech%ISR+1)
        
		! Select a particle from the ensemble using the
		! given property index for weighting.
		If (rxn%ID < 0) Then
			i = CustomChooseParticle(soln%Ensemble, rxn%P(1),rxn%P(2),rxn%P(3),rxn%P(4))
		Else
			i = ChooseParticle(soln%Ensemble, rxn%ID)
!			i = ChooseTypedParticle(soln%Ensemble, rxn%ID, rxn%PTypeIn)
		End If

		If (i < 0) Then
			! Failed to choose a particle.
			flag = CHOOSE_PARTICLE_ERR
		Else If (GetParticle(soln%Ensemble, i, sp) < 0) Then
			! Could not get details of chosen particle.
			flag = GET_PARTICLE_ERR
		Else
			! Calculate the majorant rate for this
			! particle (used to check for ficticious
			! events).
			maj = SingleParticleRate(rxn, sp, chem, soln, mech, mech%AnyDeferred, flag)
			If (flag < 0) Return

			! Perform the deferred events.
			If (mech%AnyDeferred) Then
				Call UpdateDeferred(sp, t, soln, mech, flag)
				If (flag < 0) Return
			End If

			If (SizeP(sp) < MIN_PARTICLE_SIZE) Then 
				! Particle has been oxidised out.
				If (RemoveParticle(soln%Ensemble, i) < 0) Then
					! Failed to remove particle from ensemble.
					flag = REMOVE_PARTICLE_ERR
				End If
			ElseIf (maj > 0.0E0) Then
				! Calculate the new rate after deferred events have been done.
				true = SingleParticleRate(rxn, sp, chem, soln, mech, .False., flag)
				If (flag < 0) Return
		
				! Check for a ficticious event.
				If (.Not. Ficticious(true, maj)) Then
					! Reaction goes ahead.
                    Call AdjustParticleComponents(sp, rxn%dComp, rxn%dSID, mech)
					sp%TypeID = rxn%PTypeOut
                    sp%Track = sp%Track + rxn%dTrack

					If (SizeP(sp) < MIN_PARTICLE_SIZE) Then
						! Particle has been oxidised out.
						If (RemoveParticle(soln%Ensemble, i) < 0) Then
							! Failed to remove particle from ensemble.
							flag = REMOVE_PARTICLE_ERR
						End If
					Else
						If (ReplaceParticle(soln%Ensemble, i, sp) < 0) Then
							! Failed to replace the particle.
							flag = REPLACE_PARTICLE_ERR
							Return
						End If
					End If

					! Update chemistry.
                    If (CHEM_TYPE == VARIABLE_CHEM) Then
				    Call UpdateChemistry(t, rxn, 1, soln%Chemistry, SampleVolume(soln, chem(soln%Chemistry%Shared%iT)), flag)
                    End If
				Else
					flag = FICTICIOUS_EVENT
				End If
            Else
				If (ReplaceParticle(soln%Ensemble, i, sp) < 0) Then
					! Failed to replace the particle.
					flag = REPLACE_PARTICLE_ERR
					Return
                Else
                    flag = FICTICIOUS_EVENT
				End If                
			End If
		End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine DSACoagulate(rateID1, rateID2, t, chem, soln, mech, maj, flag)
		! DESCRIPTION:
		!	Coagulates two particles under direct simulation, 
		!	returns true on success, false otherwise.  The rate
		!	indices control which set of weights are used to pick
		!	a particle in the binary tree, the values are passed
		!	directly to choose_particle.  maj is passed to
		!	fictitious_jump to give some information about the 
		!	majorant kernel to use.  chem is a structure containing
		!	data on the conditions.

        Use SWPMECH
        Use SWPSOLN
		Use SWPENSEMBLE
        Use SWPCOAG
		Implicit None

		! ARGUMENTS.
		! Rate indices controlling the weights used to 
		! select particles from the ensemble.
		Integer, Intent(IN)          :: rateID1, rateID2
		Real, Intent(IN)          :: t, chem(:) ! Current time and chemical conditions.
        Type(Solution), Intent(INOUT) :: soln       ! Particle solution.
        Type(Mechanism), Intent(IN)   :: mech       ! Mechanism to use.
        Integer, Intent(IN)           :: maj        ! Majorant to use ID.
		Integer, Intent(OUT)          :: flag		! Error flag.

		! VARIABLES.
		Integer				::	guard, i1, i2, err, i
		Type(StochParticle)	::	sp1, sp2, sp1old, sp2old
		Logical				::	p1_oxidised, p2_oxidised
		Real				::	majk, truek, Tmp, P

		! EXECUTABLE CODE.
		p1_oxidised	= .False.
		p2_oxidised	= .False.
		flag		= 0

		! Save the current temperature and pressure, in
		! order to calculate the majorant kernel after
		! deferred processes have been performed.
		Tmp = chem(soln%Chemistry%Shared%iT)
		P   = chem(soln%Chemistry%Shared%iP)

		! ***********************************
		! GET FIRST PARTICLE.
		! ***********************************

		! Select first particle.
		i1 = ChooseParticle(soln%Ensemble, rateID1)

		If (i1 > 0) Then
			If (GetParticle(soln%Ensemble, i1, sp1) == 0) Then
				sp1old = sp1

				Call UpdateDeferred(sp1, t, soln, mech, flag) 
				If (flag < 0) Return

				If (SizeP(sp1) < MIN_PARTICLE_SIZE) Then
					! Particle has been oxidised out.
					If (RemoveParticle(soln%Ensemble, i1) >= 0) Then
						i1 = -1
						p1_oxidised = .True.
					Else
						! Failed to remove particle from
						! ensemble.
						flag = REMOVE_PARTICLE_ERR
					End If
				End If
			Else
				! Could not get first particle.
				flag = GET_PARTICLE_ERR
			End If
		Else
			! Could not choose first particle.
			flag = CHOOSE_PARTICLE_ERR
		End If

		! ***********************************
		! GET SECOND PARTICLE.
		! ***********************************

		i2  = i1
		guard = 0

		Do While (i2 == i1 .And. guard < 1000)
			! Choose a second particle distinct from the first.
			i2 = ChooseParticle(soln%Ensemble, rateID2)
            guard = guard + 1
        End Do

		If ((i2 > 0) .And. (i2 /= i1)) Then
			If (GetParticle(soln%Ensemble, i2, sp2) == 0) THEN
				sp2old = sp2

				Call UpdateDeferred(sp2, t, soln, mech, flag)
				If (flag < 0) Return

				If (SizeP(sp2) < MIN_PARTICLE_SIZE) Then
					! Particle has been oxidised out.

					If (.Not. p1_oxidised) Then
						! We need to put the updated particle 1 
						! back in the tree before we invalidate
						! its index, p1, by removing particle 2
						! from the tree.
						If (ReplaceParticle(soln%Ensemble, i1, sp1) < 0) Then
							! Failed to replace particle in the
							! ensemble.
							flag = REPLACE_PARTICLE_ERR
							Return
						End If
					End If
					
					If (RemoveParticle(soln%Ensemble, i2) >= 0) Then
						i2 = 0
						p2_oxidised = .True.
					Else
						! Failed to remove particle from
						! ensemble.
						flag = REMOVE_PARTICLE_ERR
						Return
					End If
				End If
			Else
				! Could not get second particle.
				flag = GET_PARTICLE_ERR
				Return
			End If
		Else
			! Could not choose second particle.
			flag = CHOOSE_PARTICLE_ERR
			Return
		End If

		! ***********************************
		! COAGULATE PARTICLES.
		! ***********************************

		If (flag >= 0) Then
			If (.Not. (p1_oxidised .Or. p2_oxidised)) Then
				! Calculate the majorant and true kernels.
				majk  = CoagKernel(sp1old, sp2old, Tmp, P, maj)
				truek = CoagKernel(sp1, sp2, chem(soln%Chemistry%Shared%iT), chem(soln%Chemistry%Shared%iP), NOMAJ)

				If (Ficticious(truek, majk)) Then
					! This is a fictitious jump.
					If (ReplaceParticle(soln%Ensemble, i1, sp1) < 0 .Or. &
						ReplaceParticle(soln%Ensemble, i2, sp2) < 0) Then
						! Failed to replace one of the
						! particles in the ensemble.
						flag = REPLACE_PARTICLE_ERR
						Return
					End If

					flag = FICTICIOUS_EVENT
				Else
                    If (mech%CoagRuleCount > 0) Then
				    ! Coagulation (rule dependent).
                        Do i = 1, mech%CoagRuleCount
                            If (((sp1%TypeID == mech%CoagRules(i)%In1) .And. (sp2%TypeID == mech%CoagRules(i)%In2)) .Or. &
                                ((sp1%TypeID == mech%CoagRules(i)%In2) .And. (sp2%TypeID == mech%CoagRules(i)%In1))) Then
                                sp1old = sp1
				            sp1 = Combine(sp1, sp2, mech)
                                sp1%TypeID = mech%CoagRules(i)%Out

				            If (ReplaceParticle(soln%Ensemble, i1, sp1) == 0) Then
					            If (RemoveParticle(soln%Ensemble, i2) < 0) Then
						            ! Failed to remove particle.
						            flag = REMOVE_PARTICLE_ERR
					            End If
				            Else
					            ! Failed to store coagulated particle.
					            flag = REPLACE_PARTICLE_ERR
				            End If
                               Exit
                            End If
                        End Do
                        If (i > mech%CoagRuleCount) flag = FICTICIOUS_EVENT
                    Else
                        ! No coagulation rules.
                        sp1old = sp1
				    sp1 = Combine(sp1, sp2, mech)

				    If (ReplaceParticle(soln%Ensemble, i1, sp1) == 0) Then
					    If (RemoveParticle(soln%Ensemble, i2) < 0) Then
						    ! Failed to remove particle.
						    flag = REMOVE_PARTICLE_ERR
					    End If
				    Else
					    ! Failed to store coagulated particle.
					    flag = REPLACE_PARTICLE_ERR
				    End If
                    End If
				End If  
			Else
				If (p1_oxidised .And. .Not. p2_oxidised) Then
					If (ReplaceParticle(soln%Ensemble, i2, sp2) < 0) Then
						! Failed to replace particle.
						flag = REPLACE_PARTICLE_ERR
					End If
				End If

				! We do not need to handle the following case because if this arises the updating
				! of particle 1 is performed prior to the removal of particle 2 because the removal
				! process invalidates the index p1 which gives the location of p1 in the tree.
				! IF(p2_oxidised .AND. .NOT. p1_oxidised) THEN
			End If
		End If
	End Subroutine

	! -------------------------------------------------------

	Logical Function Ficticious(truek, majk)
		! DESCRIPTION:
		!	Takes a true rate and a majorant
		!	rate calculated for a single particle and
		!	determines whether or not they represent
		!	a ficticious event.
		! RETURNS:
		!	True if event is ficticious, otherwise
		!	False.
		Use SWPRNG
		Implicit None
		Real, Intent(IN)::	truek, majk ! True and majorant rates.
		Ficticious = .Not. (majk * Rnd() < truek)
	End Function

	! -------------------------------------------------------
	
	Subroutine Inflow(soln, flag)
		! DESCRIPTION:
		!	Chooses a particle from the inflowing particle
		!	list and adds it to the ensemble.

        USe SWPSOLN
		Implicit None

		! ARGUMENTS.
        Type(Solution), Intent(INOUT) :: soln ! Solution to modify.
		Integer, Intent(OUT)          :: flag ! Error flag.

		! VARIABLES.

		! EXECUTABLE CODE.
		flag = 0

		! THIS FUNCTION IS NOT YET COMPLETE!!
	End Subroutine

	! -------------------------------------------------------
	
	Subroutine Outflow(t, soln, mech, flag)
		! DESCRIPTION:
		!	Chooses a particle from the particle ensemble
		!	uniformly and removes it.

        USe SWPMECH
        Use SWPSOLN
		Use SWPENSEMBLE
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)	      :: t    ! Current time.
        Type(Solution), Intent(INOUT) :: soln ! Solution to modify.
        Type(Mechanism), Intent(IN)   :: mech ! Mechanism to use.
		Integer, Intent(OUT)      :: flag ! Error flag.

		! VARIABLES.
        Real                :: N
		Integer				:: i, err
		Type(StochParticle)	:: sp

		! EXECUTABLE CODE.
		flag = 0
        N = Real(ParticleCount(soln%Ensemble))

        If (OUTFLOW_METHOD == OUTFLOW_SCALE_METHOD) Then
            ! Do not actually remove particle, just scale
            ! sample volume so that it appears a particle
            ! was removed.
            If (N > 1) Then
                soln%BaseVolume = soln%BaseVolume * N / (N-1.0E0)
            Else
                ! Remove the only particle in the system.
                ! Note, it should be impossible to get here
                ! if N==0 as the outflow rate contains N.
	        If (RemoveParticle(soln%Ensemble, 1) < 0) Then
		        flag = REMOVE_PARTICLE_ERR
	        End If
            End If
        ElseIf (OUTFLOW_METHOD == OUTFLOW_PARTICLE_METHOD) Then
            ! Actually remove a particle from the system, chosen
            ! uniformly.
	    i = ChooseParticle(soln%Ensemble, 0)

	    If (i < 0) Then
		    ! Failed to choose a particle.
		    flag = CHOOSE_PARTICLE_ERR
	    Else
                If (GetParticle(soln%Ensemble, i, sp) < 0) Then
		        ! Could not get details of chosen particle.
		        flag = GET_PARTICLE_ERR
	        Else
		        ! Perform the deferred events in order to update
		        ! the chemistry.
		        If (SURF_ON .And. mech%AnyDeferred) Then
			        Call UpdateDeferred(sp, t, soln, mech, flag)
			        If (flag < 0) Return
		        End If

	            If (RemoveParticle(soln%Ensemble, i) < 0) Then
		            flag = REMOVE_PARTICLE_ERR
	            End If
	        End If
            End If
        End If
	End Subroutine

	! -------------------------------------------------------
	! DEFERRED PROCESSES.
	!
	!	These routines describe the process deferment (LPDA)
	!	mechanism.
	!
	! -------------------------------------------------------

	Subroutine SurfaceIntegrate(t1, t2, soln, mech, flag)
		! DESCRIPTION:
		!	Retrieves the particle list from the ensemble,
		!	performs surface growth by calling an integration
		!	routine and then puts any remaining particles back into
		!	the tree.
		!
		!	The return value is negative on error and if >= 0 is
		!	the number of particles oxidised out of the system
		!	during the integration.

        Use SWPMECH
        Use SWPSOLN
		Use SWPENSEMBLE
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)	      :: t1, t2 ! Start and stop times of integration (s).
        Type(Solution), Intent(INOUT) :: soln   ! Solution to modify.
        Type(Mechanism), Intent(IN)   :: mech   ! Mechanism to use.
		Integer, Intent(OUT)      :: flag   ! Error flag.

		! VARIABLES.
		Type(StochParticle) :: particles(soln%Ensemble%FirstSpace-1), sp
		Integer				:: N, N2, i, j

		! EXECUTABLE CODE.
		flag = 0
		If (.Not. SURF_ON) Return

		! Get the list of particles from the
		! ensemble.
		Call GetTypedList(soln%Ensemble, N, particles)

		If (N == 0) Then
			! There are no particles to
			! integrate.
			Return
		End If

		! Loop over all the particles and update
		! deferred events.
		Do i = 1, N
			Call UpdateDeferred(particles(i), t2, soln, mech, flag)
			If (flag < 0) Return
		End Do

		! Save scaling factor as it will be invalidated
		! when updated particles are replaced into ensemble.
        soln%LastKnownScale = ScalingFactor(soln%Ensemble)
		soln%BaseVolume = soln%BaseVolume * soln%LastKnownScale

		! Count how many particles still have at
		! least 32 C atoms - i.e. have not been oxidised
		! out of the model.
		N2 = Count(SizeP(particles) >= MIN_PARTICLE_SIZE)

		If (N2 > 0) Then
			If (N2 < N) Then
				! Some particles have been oxidised out.

				! This section is intended to move all the
				! particles that have been oxidised out of the
				! model to the far end of the array so we are left
				! with a contiguous section at the low end of particles
				! to put back into the binary tree.  The idea is to
				! search from the beginning for oxidised particles to move
				! to the end and to search from the end for non-oxidised
				! particles that can be exchanged with oxidised particles
				! found nearer the start.
				!
				! - Basically it is a quick sort.

				i = 1
				j = N

				Do While (i < j)
					Do While (SizeP(particles(i)) >= MIN_PARTICLE_SIZE .AND. i <= N2) 
						i = i + 1
					End Do
					Do While (SizeP(particles(j)) < MIN_PARTICLE_SIZE .AND. j > 0)
						j = j - 1
					End Do

					If (i < j) Then
						! Swap particles i and j.
						sp = particles(j)
						particles(j) = particles(i)
						particles(i) = sp
					Else
						! All the oxidised particles have been
						! moved to the end.
						Exit
					End If
				End Do
			End If

			! Any oxidised particles will be at the end of the
			! array now, so we can load the remaining particles
			! back into the binary tree.
			Call SetTypedList(N2, particles, soln%Ensemble, flag)

			If (flag == N2) Then
                flag = 0 ! N - N2
			Else
				! Could not load particles back into
				! tree after integration.
				flag = SET_LIST_ERR
			End If
		Else
			! No particles remain after integration.
			Call ClearEnsemble(soln%Ensemble)
		End If
	End Subroutine
	
	! -------------------------------------------------------
    
    Subroutine UpdateDeferred(sp, t, soln, mech, flag)
		! DESCRIPTION:
		!	Performs all deferred events on a stochastic
		!	particle.

        Use SWPMECH
        Use SWPSOLN
		Use SWPCHEM
		Use SWPENSEMBLE
		Use SWPRNG
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp   ! Particle to update.
		Real, Intent(IN)               :: t    ! Current time.
        Type(Solution), Intent(INOUT)      :: soln ! Solution to modify.
        Type(Mechanism), Intent(IN)        :: mech ! Mechanism to use.
		Integer, Intent(OUT)               :: flag ! Error flag.

		! VARIABLES.
		Integer	::	num, i, loopCount, err
		Real	::	chem(soln%Chemistry%Shared%NCHEM), dt, rate, tup1

		! EXECUTABLE CODE.
		flag = 0
		loopCount = 0
        tup1 = sp%LastUpdate

		If (.Not. (SURF_ON .Or. SINT_ON)) Return

		! Loop through time steps.
		Do While (sp%LastUpdate < t .And. SizeP(sp) >= MIN_PARTICLE_SIZE)

			! Get chemistry (it might have been changed by other
			! surface reactions.
			chem = GetChem(soln%Chemistry, sp%LastUpdate, t, flag)
			If (flag < 0) Return

            If (LPDA_TWO_LEVELS) Then
		    rate = Abs(VolChangeRate(sp, chem, soln, mech, err))

		    If (rate > 1.0E-6) Then
			    ! Take a step that does not change the
			    ! particle volume by more than about 10%.
			    dt = sp%Properties(iV) * 0.1E0 / rate

			    If (dt > 5.0E0 * DataSpacing(soln%Chemistry, sp%LastUpdate)) Then
				    dt = 5.0E0 * DataSpacing(soln%Chemistry, sp%LastUpdate)
			    End If
		    Else
			    ! If the volume change rate is approximately
			    ! 0, either due to nothing much happening or
			    ! lots of things cancelling each other out
			    ! take a small step.
			    dt = DataSpacing(soln%Chemistry, sp%LastUpdate)
		    End If

		    ! Also check that the time step does not go past stop_time and truncate
		    ! if this is a problem.  In the case of truncation it is important that
		    ! stop_time is assigned to spart%last_update in such a way as to ensure
		    ! the test at the top of the loop will return false when it is next
		    ! evaluated.  This should not be a problem if they are both of the same
		    ! precision.
		    If (t - sp%LastUpdate < dt) Then
			    dt = t - sp%LastUpdate
			    sp%LastUpdate = t
		    Else
			    sp%LastUpdate = sp%LastUpdate + dt
		    End If
            Else
		    dt = t - sp%LastUpdate
		    sp%LastUpdate = t
            End If

			! Now calculate the number of each event and perform them.  Check
			! particle volumes do not drop below 0 to avoid problems with 
			! exponentiation in the rate calculations and any other numerical
			! traps.  Such particles are unlikely to grow enough to stay in the
			! ensemble anyway.
			Do i = 1, mech%ReactionCount
				If (SizeP(sp) < 0) Exit

				If (mech%DeferMask(mech%ISR+i-1)) Then
					! Calculate process rate for this particle.
					rate = SingleParticleRate(mech%Reactions(i), sp, chem, soln, mech, .False., flag)
					If (flag < 0) Return

					If (rate > 0.0E0) Then
                        ! Poisson process.
						num = ignpoi(rate * dt)

						If (num > 0) Then
							Call DeferredSurface(t, i, sp, num, soln, mech, &
                                                 SampleVolume(soln, chem(soln%Chemistry%Shared%iT)), flag)

					    If (flag < 0) Then
						    Return
					    Else
						    soln%ProcessCounter(mech%ISR+i-1) = &
                                soln%ProcessCounter(mech%ISR+i-1) + num
					    End If
						End If
					End If
				End If
			End Do

            ! Deferred sintering.
            If (SINT_ON) Then
                rate = ParticleSinteringRate(sp, mech, chem(soln%Chemistry%Shared%iT))
                num  = ignpoi(rate * dt)
                rate = EquivSphereSurface(sp, mech%Components, mech%ComponentCount)
                sp%Surface = Max(sp%Surface - (Real(num) * SINT_SCALE * rate), rate)
            End If

			! Increment the loop count.
			loopCount = loopCount + 1
		End Do

		If (SizeP(sp) >= MIN_PARTICLE_SIZE) Then
			! No point in doing anything if particle
			! is to be removed.
			Call CalcProperties(sp, mech)
		End If 
	End Subroutine
	
	! -------------------------------------------------------

	Subroutine DeferredSurface(t, i, sp, repeat, soln, mech, vol, flag)
		! DESCRIPTION:
		!	Performs a deferred surface event a number of
		!	times on a given particle.

        Use SWPMECH
        Use SWPSOLN
		Use SWPPART
        Use SWPCHEM
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)               :: t     ! Current flow time (s).
		Integer, Intent(IN)               :: i      ! The surface reaction being performed.
		Type(StochParticle), Intent(INOUT) :: sp     ! The particle to update.
		Integer, Intent(IN)               :: repeat ! The number of times to apply this reaction.
        Type(Solution), Intent(INOUT)      :: soln   ! Solution to modify.
        Type(Mechanism), Intent(IN)        :: mech   ! Mechanism to use.
        Real, Intent(IN)                   :: vol    ! Current sample volume.
		Integer, Intent(OUT)               :: flag   ! Error flag.

		! EXECUTABLE CODE.
		flag = 0
	
		! Apply reaction to soot particle.
        Call AdjustParticleComponents(sp, mech%Reactions(i)%dComp * repeat, &
                                      mech%Reactions(i)%dSID, mech)

		If (SizeP(sp) > 0.0E0) Then
	    sp%TypeID = mech%Reactions(i)%PTypeOut
            sp%Track  = sp%Track + (mech%Reactions(i)%dTrack * Real(repeat))
        ElseIf (SizeP(sp) < 0.0E0) Then
            Print *, "waaah!"
        End If

		! Update chemistry.
        If (CHEM_TYPE == VARIABLE_CHEM) Then
	    Call UpdateChemistry(t, mech%Reactions(i), repeat, soln%Chemistry, vol, flag)
        End If
	End Subroutine
	
	! -------------------------------------------------------
	! CHEMISTRY FUNCTIONS.
	!
	!	These routines control how the mechanism interacts
	!	with the gas-phase chemistry.
	!
	! -------------------------------------------------------

	Subroutine UpdateChemistryI(t, rxn, repeat, chem, volume, flag)
		! DESCRIPTION:
		!	Applies the effects of a soot inception to the
		!	gas-phase chemistry stored in the chemistry
		!	module.

        Use SWPMECH
		Use SWPCHEM
		Use SWPPARAMS, only: NA
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)		  :: t      ! Time of the change (s).
		Type(Inception), Intent(IN)  :: rxn  ! The inception being performed.
		Integer, Intent(IN)		  :: repeat ! The number of times to perform the update.
        Type(ChemData), Intent(INOUT) :: chem   ! Chemistry to modify.
		Real, Intent(IN)		  :: volume ! Current sample volume.
		Integer, Intent(OUT)	  :: flag   ! Error flag.

		! VARIABLES.
		Integer :: i
		Real	:: scale

		! EXECUTABLE CODE.
		flag = 0

        If (CHEM_TYPE == VARIABLE_CHEM) Then
	    ! Calculate the scaling from ensemble level to molar
	    ! concentration level.
	    scale = Real(repeat) / (volume * NA)

	    ! Update chemistry accordingly.
	    Do i = 1, rxn%NREAC
		    Call ChangeChem(chem, t, rxn%Species(i), - rxn%Stoich(i) * scale)
	    End Do
	    Do i = rxn%NREAC + 1, rxn%NREAC + rxn%NPROD
		    Call ChangeChem(chem, t, rxn%Species(i), + rxn%Stoich(i) * scale)
	    End Do
        End If
	End Subroutine
	
	! -------------------------------------------------------

	Subroutine UpdateChemistryS(t, rxn, repeat, chem, volume, flag)
		! DESCRIPTION:
		!	Applies the effects of a soot reaction to the
		!	gas-phase chemistry stored in the chemistry
		!	module.

        Use SWPMECH
		Use SWPCHEM
		Use SWPPARAMS, only: NA
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)		  :: t      ! Time of the change (s).
		Type(Reaction), Intent(IN)    :: rxn	! The reaction being performed.
		Integer, Intent(IN)		  :: repeat ! The number of times to perform the update.
        Type(ChemData), Intent(INOUT) :: chem   ! Chemistry to modify.
		Real, Intent(IN)		  :: volume ! Current sample volume.
		Integer, Intent(OUT)	  :: flag   ! Error flag.

		! VARIABLES.
		Integer :: i
		Real	:: scale

		! EXECUTABLE CODE.
		flag = 0

        If (CHEM_TYPE == VARIABLE_CHEM) Then
	    ! Calculate the scaling from ensemble level to molar
	    ! concentration level.
	    scale = Real(repeat) / (volume * NA)

	    ! Update chemistry accordingly.
	    Do i = 1, rxn%NREAC
		    Call ChangeChem(chem, t, rxn%Species(i), - rxn%Stoich(i) * scale)
	    End Do
	    Do i = rxn%NREAC + 1, rxn%NREAC + rxn%NPROD
		    Call ChangeChem(chem, t, rxn%Species(i), + rxn%Stoich(i) * scale)
	    End Do
        End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine UpdateChemistrySV(t, rxn, repeat, chem, volume, flag)
		! DESCRIPTION:
		!	Applies the effects of a soot reaction to the
		!	gas-phase chemistry stored in the chemistry
		!	module.

        Use SWPMECH
		Use SWPCHEM
		Use SWPPARAMS, only: NA
		Implicit None

		! ARGUMENTS.
		Real, Intent(IN)		  :: t       ! Time of the change (s).
		Type(Reaction), Intent(IN)    :: rxn ! The reaction being performed.
		Integer, Intent(IN)		  :: repeat  ! The number of times to perform the update.
        Real, Intent(INOUT)           :: chem(:) ! Chemistry to modify.
		Real, Intent(IN)		  :: volume  ! Current sample volume.
		Integer, Intent(OUT)	  :: flag    ! Error flag.

		! VARIABLES.
		Integer :: i
		Real	:: scale

		! EXECUTABLE CODE.
		flag = 0

        If (CHEM_TYPE == VARIABLE_CHEM) Then
	    ! Calculate the scaling from ensemble level to molar
	    ! concentration level.
	    scale = Real(repeat) / (volume * NA)

	    ! Update chemistry accordingly.
	    Do i = 1, rxn%NREAC
                chem(rxn%Species(i)) = Max(0.0E0, chem(rxn%Species(i)) - (rxn%Stoich(i) * scale))
	    End Do
	    Do i = rxn%NREAC + 1, rxn%NREAC + rxn%NPROD
                chem(rxn%Species(i)) = chem(rxn%Species(i)) + (rxn%Stoich(i) * scale)
	    End Do
        End If
	End Subroutine
End Module

! *****************************************************************************
!
! File:					swppart.f90
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
!	Definition of a stochastic particle as they are known to Sweep.
!
! Functions:
!	-- (Operator routines) ----------------------------------------------------
!   SetComponents       -   Sets the components of a stochastic particle.
!   SetAttributes       -   Sets the particle properties which are not functions
!                           of other properties.
!   SetCoordsMulti      -   As SetCoords() except builds multiple particles with the
!                           same coordinates.
!   ParticleFromArray   -   Builds a stochastic particle from a Real array.
!   ParticleToArray     -   Creates a Real array which fully describes a particle.
!   Combine             -   A mechanism independent method for combining two
!                           particles (coagulation).
!   AdjustParticleComponents  -   Adjusts a particle with the given changes to
!                                 to particle components.
!	-- (Particle property routines) -------------------------------------------
!   GetPreCalcs         -   Returns the pre-calculated properties of a particle.
!   CreateParticle      -   Creates a new particle with the given unique
!                           properties.  Also sets default values for other coordinates.
!   CalcProperties      -   Updates the pre-calculated values stored per particle
!                           from the unique coordinates.
!	-- (Particle properties) --------------------------------------------------
!   Mass                -   Returns the mass of a particle.
!   Volume              -   Returns the volume of a particle.
!   SurfaceArea         -   Returns the surface area of a particle.
!   EquivSphereSurface  -   Returns the surface area of an equivalent volume
!                           spherical particle.
!   EquivSphereDiameter -   Returns the diameter of an equivalent volume spherical
!                           particle.
!   CollisionDiameter   -   Returns the collision diameter of a particle.
!   GrowthRadius        -   Returns the radius used to calculate dS/dV for growth
!                           processes.
!   OxidationRadius     -   Returns the radius used to calculate dS/dV for oxidation/loss
!                           processes.
!	-- (Process rate helper routines) -----------------------------------------
!   CalcParticleWeight  -   Calculates a custom particle rate weighting given the
!                           property exponents.
!   GetPropertyIndex    -   Returns the index of a particle property based on the
!                           powers of standard values (mass, surface etc).
!	-- (Statistics routines) --------------------------------------------------
!   GetStatItems        -   Returns an array of useful statistical properties of
!                           a particle.
! *****************************************************************************

! Indices of particle statistical items, returned
! by the GetStatItems() routine.
Module SWPPART_STATS
	Implicit Integer (i,C), Character*50 (N)
	Public
	Parameter (Count=12,iM0=1,iM1=2,iM2=3,iM3=4,iM4=5,iM5=6,&
			   iM6=7,iV=8,iMass=9,iSurf=10,iActSurf=11,iDiam=12)
	Character*50, Parameter :: Names(Count)=(/"M0","M1","M2","M3","M4","M5", &
			                                  "M6","Volume (cm3)","Mass (g)", &
                                              "Surface (cm2)","Alpha", "Diameter (cm)"/)
End Module

! =======================================================

Module SWPPART
    Use SWPMECH_TYPES, only: MAX_COMP, MAX_TRACK
	Implicit None
	Public

	! -------------------------------------------------------
	! PARTICLE DEF'N PARAMETERS.
	! -------------------------------------------------------

	! Property indices - These indices are used to get specific
	! properties from a particles Properties data member.
    ! NOTE:  The last properties are the unique coordinates, and
    !        therefore these indices only become valid if you get
    !        the values using the GetPreCalcs() function.
	Integer, Parameter, Private :: CACHE_COUNT    = 10
	Integer, Parameter, Public  :: PROPERTY_COUNT = CACHE_COUNT + 1
	Integer, Parameter :: iAS      = 1 ! Active surface area.
    Integer, Parameter :: iD       = 2 ! Collision diameter.
	Integer, Parameter :: iD2      = 3 ! Collision diameter squared.
	Integer, Parameter :: iD_1     = 4 ! 1 / Collision diameter.
	Integer, Parameter :: iD_2     = 5 ! 1 / Collision diameter squared.
	Integer, Parameter :: iM_1_2   = 6 ! Mass to the -1/2.
	Integer, Parameter :: iD2M_1_2 = 7 ! Collision diameter squared * mass to the -1/2.
    Integer, Parameter :: iV       = 8 ! Volume (required to get total volume).
    Integer, Parameter :: iM       = 9 ! Mass (required to get total mass).
    Integer, Parameter :: iSz      = 10 ! Particle "size".
	Integer, Parameter :: iS	   = 11 ! Surface area.


    ! Length of a particle array were it converted into a list
    ! of reals.
    Integer, Parameter, Public :: PARTICLE_ARRAY_SIZE = MAX_COMP + CACHE_COUNT + MAX_TRACK + 5

	! -------------------------------------------------------
	! PARTICLE DEF'N.
	! -------------------------------------------------------

	Type StochParticle
        ! Particle type ID.
        Integer :: TypeID = 0
        ! Counts of smallest particle components.
        Integer :: Components(MAX_COMP) = 0
        ! Particle attributes: Properties which are not
        ! functions other particle properties.
        Real :: Surface    = 0.0E0 ! Surface area (for surface-volume model).
        Real :: CreateTime = 0.0E0 ! Time at which particle was incepted.
		Real :: LastUpdate = 0.0E0 ! Time at which particle was last changed.
		Real :: ActSurf    = 1.0E0 ! Fraction active sites; Alpha (ABF).
        ! Properties which depend entirely on other
        ! particle properties.  Pre-calculated to
        ! improve program run time.
        Real :: Properties(CACHE_COUNT) = 0.0E0
        ! Custom tracking variables defined in
        ! the mechanism file.
		Real :: Track(MAX_TRACK) = 0.0
	End Type StochParticle

	! -------------------------------------------------------
	! PARTICLE OPERATORS.
	! -------------------------------------------------------

	! Assignment operator.
	Interface Assignment(=)
		Module Procedure ParticleFromArray
		Module Procedure ParticlesFromArray
		Module Procedure ParticleToArray
	End Interface

	! -------------------------------------------------------
	! ROUTINE OVERLOADS.
	! -------------------------------------------------------

    Interface SetAttributes
        Module Procedure SetAttributes_Array
        Module Procedure SetAttributes_Piecewise
    End Interface

    Interface SizeP
        Module Procedure SizeP_Part
        Module Procedure SizeP_PartArray
        Module Procedure SizeP_Comp
    End Interface

    Interface Mass
        Module Procedure Mass_Part
        Module Procedure Mass_Comp
    End Interface

    Interface Volume
        Module Procedure Volume_Part
        Module Procedure Volume_Comp
    End Interface

    Interface EquivSphereSurface
        Module Procedure EquivSphereSurface_Part
        Module Procedure EquivSphereSurface_Comp
        Module Procedure EquivSphereSurface_Vol
    End Interface

    Interface EquivSphereDiameter
        Module Procedure EquivSphereDiameter_Part
        Module Procedure EquivSphereDiameter_Vol
    End Interface

    Interface CollisionDiameter
        Module Procedure CollisionDiameter_Part
        Module Procedure CollisionDiameter_SurfVol
    End Interface

    Interface GrowthRadius
        Module Procedure GrowthRadius_Part
        Module Procedure GrowthRadius_Surf
    End Interface

    Interface OxidationRadius
        Module Procedure OxidationRadius_Part
        Module Procedure OxidationRadius_SurfVol
    End Interface

	Contains

	! -------------------------------------------------------
	! PARTICLE OPERATOR ROUTINES.
	!
	!	The following routines define the programmatical
	!	operators for stochastic particles.
	!
	! -------------------------------------------------------

	Subroutine SetComponents(sp, values, mech)
		! DESCRIPTION:
		!	Sets the components of a stochastic particle.
        !   Also sets values for other particle attributes
        !   and properties under the assumption that the
        !   particle is spherical and homogeneous.
		!
		!	This routine is used as an interface for
		!	the assignment operator for stochastic
		!	particles.

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp ! Particle to set.
		Integer, Intent(IN)         :: values(:) ! Values to set.
        Type(Mechanism), Intent(IN) :: mech      ! Mechanism in use.

        ! VARIABLES
        Integer :: N

		! EXECUTABLE CODE.
		N = Min(MAX_COMP, Size(values))
        sp%Components(1:N) = values(1:N)

        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
            sp%Surface = EquivSphereSurface(sp, mech%Components, mech%ComponentCount)
        End If

		Call CalcProperties(sp, mech)
	End Subroutine

	! -------------------------------------------------------

	Subroutine SetAttributes_Array(sp, attr, mech)
		! DESCRIPTION:
		!	Sets particle attributes other than composition..

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp
		Real, Intent(IN)                   :: attr(:)
        Type(Mechanism), Intent(IN)        :: mech ! Mechanism in use.

		! EXECUTABLE CODE.
        If (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
		    sp%Surface = attr(1)
        End If
        sp%CreateTime  = attr(2)
		sp%LastUpdate  = attr(3)
		sp%ActSurf     = attr(4)
		Call CalcProperties(sp, mech)
	End Subroutine

	! -------------------------------------------------------

	Subroutine SetAttributes_Piecewise(sp, surf, createt, lastt, actsurf, mech)
		! DESCRIPTION:
		!	Sets particle attributes other than composition..

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp
		Real, Intent(IN)                   :: surf, createt, lastt, actsurf
        Type(Mechanism), Intent(IN)        :: mech ! Mechanism in use.

		! EXECUTABLE CODE.
        If (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
		    sp%Surface = surf
        End If
        sp%CreateTime  = createt
		sp%LastUpdate  = lastt
		sp%ActSurf     = actsurf
		Call CalcProperties(sp, mech)
	End Subroutine

	! -------------------------------------------------------

	Subroutine ParticleFromArray(sp, arr)
		! DESCRIPTION:
		!	Builds a particle from a Real array.
		!
		!	This routine is used as an interface for
		!	the assignment operator for stochastic
		!	particles.

		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(OUT) :: sp
		Real, Intent(IN) :: arr(:)

		! EXECUTABLE CODE.
        If (Size(arr) >= PARTICLE_ARRAY_SIZE) Then
            sp%TypeID     = Int(arr(1))
            sp%Components = Int(arr(2:MAX_COMP+1))
            sp%Surface    = arr(MAX_COMP+2)
            sp%CreateTime = arr(MAX_COMP+3)
            sp%LastUpdate = arr(MAX_COMP+4)
            sp%ActSurf    = arr(MAX_COMP+5)
            sp%Properties = arr(MAX_COMP+6:MAX_COMP+CACHE_COUNT+5)
            sp%Track      = arr(MAX_COMP+CACHE_COUNT+6:MAX_COMP+CACHE_COUNT+MAX_TRACK+5)
        End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine ParticlesFromArray(sp, arr)
		! DESCRIPTION:
		!	Builds a particle from a Real array.
		!
		!	This routine is used as an interface for
		!	the assignment operator for stochastic
		!	particles.

		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(OUT) :: sp(:)
		Real, Intent(IN) :: arr(:)

        ! VARIABLES.
        Integer :: i, N

		! EXECUTABLE CODE.
        If (Size(arr) >= PARTICLE_ARRAY_SIZE) Then
            N = Size(sp)
            Do i = 1, N
                Call ParticleFromArray(sp(i), arr)
            End Do
        End If
	End Subroutine

	! -------------------------------------------------------

	Subroutine ParticleToArray(arr, sp)
		! DESCRIPTION:
		!	Gets unique particle coordinates.
		!
		!	This routine is used as an interface for
		!	the assignment operator for stochastic
		!	particles.

		Implicit None

		! ARGUMENTS.
		Real, Intent(OUT) :: arr(:)
		Type(StochParticle), Intent(IN) :: sp

		! EXECUTABLE CODE.
        If (Size(arr) >= PARTICLE_ARRAY_SIZE) Then
            arr(1) = Real(sp%TypeID)
            arr(2:MAX_COMP+1) = Real(sp%Components)
            arr(MAX_COMP+2) = sp%Surface
            arr(MAX_COMP+3) = sp%CreateTime
            arr(MAX_COMP+4) = sp%LastUpdate
            arr(MAX_COMP+5) = sp%ActSurf
            arr(MAX_COMP+6:MAX_COMP+CACHE_COUNT+5) = sp%Properties
            arr(MAX_COMP+CACHE_COUNT+6:MAX_COMP+CACHE_COUNT+MAX_TRACK+5) = sp%Track
        Else
            arr = 0.0E0
        End If
	End Subroutine

	! -------------------------------------------------------

	Function Combine(sp1, sp2, mech) Result (sp)
		! DESCRIPTION:
		!	Combines two stochastic particles.  Gives a
		!	particle-independent method for coagulation.
		!
		!	This function is an interface for the
		!	stochastic particle addition (+) operator.
        !
        !   Note this function invalidates the new
        !   particle's property list.
		! RETURNS:
		!	The combined particle.
        
        Use SWPPARAMS
        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle)	:: sp
		Type(StochParticle), Intent(IN)	:: sp1, sp2
        Type(Mechanism), Intent(IN) :: mech

        ! VARIABLES.
        Integer :: i
        Real :: sc2, vc2, theta, v1, vol, surf, ssph

		! EXECUTABLE CODE.
		sp%Components = sp1%Components + sp2%Components
        sp%CreateTime = (sp1%CreateTime + sp2%CreateTime) * ONE_HALF
		sp%LastUpdate = (sp1%LastUpdate * sp1%Surface + &
						 sp2%LastUpdate * sp2%Surface) / &
						(sp1%Surface + sp2%Surface)
		sp%ActSurf    = (sp1%properties(iAS) + sp2%properties(iAS)) / &
						(sp1%Surface + sp2%Surface)
        sp%Track      = sp1%Track + sp2%Track

        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
            ! Total coalescence.
            sp%Surface = EquivSphereSurface_Vol(sp1%Properties(iV) + sp2%Properties(iV))
        ElseIf (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
            ! Point contact.
            sp%Surface = sp1%Surface + sp2%Surface
        End If

        v1 = 10**(2.6) * 1.0e-21
        vc2 = 2.1443e-23
        sc2 = PI * ( 6.0*vc2/PI )**(2.0/3.0)
        vol=sp1%Properties(iV) + sp2%Properties(iV)
        if(vol<v1) then
          theta = 2.0
        else
          theta = 3.0 * ( log(vol/v1) + 2.0/3.0 * log(v1/vc2) ) / log(vol/vc2)
        endif

        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
           theta=2.0
        End if
        ssph = (vol/vc2)**(2.0/3.0) * sc2
        surf = (vol/vc2)**(theta/3.0) * sc2
! MODIFICATION SURFACE POUR IMPOSER LA LOI DE RODRIGUES ET AL 2018
!        sp%Surface = surf
!        sp%Surface = Max(sp%Surface, ssph)


        ! Now use coagulation rules to set particle type.
		sp%TypeID = 0 ! Type reset on combination!
        Do i = 1, mech%CoagRuleCount
            If (((sp1%TypeID == mech%CoagRules(i)%In1) .And. (sp2%TypeID == mech%CoagRules(i)%In2)) .Or. &
                ((sp1%TypeID == mech%CoagRules(i)%In2) .And. (sp2%TypeID == mech%CoagRules(i)%In1))) Then
               sp%TypeID = mech%CoagRules(i)%Out
               Exit
            End If
        End Do

        Call CalcProperties(sp, mech)
	End Function

	! -------------------------------------------------------

	Subroutine AdjustParticleComponents(sp, dComp, rid, mech)
		! DESCRIPTION:
		!	Adjusts the particle composition and surface area
        !   using the given change in composition and radius ID.

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp ! Particle to adjust.
		Integer, Intent(IN) :: dComp(:)          ! Change in particle components.
        Integer, Intent(IN) :: rid               ! Type of radius to use for surface-volume.
        Type(Mechanism), Intent(IN) :: mech      ! Mechanism in use.

        ! VARIABLES.
        Real :: dvol, vold, ssph
        Real :: sc2, vc2, theta, v1, vol, surf

		! EXECUTABLE CODE.
        vold = Volume(sp%Components, mech%Components, mech%ComponentCount)
        dvol = Volume(dComp, mech%Components, mech%ComponentCount)
        ssph = EquivSphereSurface(vold+dvol)
        sp%Components = sp%Components + dComp
        Where(sp%Components < 0) sp%Components = 0

        ! Adjust surface area based on particle model.
        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
            sp%Surface = ssph
        ElseIf (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
            Select Case (rid) 
                Case (GROWTH_RADIUS)
                    sp%Surface = sp%Surface + (2.0 * dvol / GrowthRadius(sp%Surface))
                !    sp%Surface = sp%Surface + 2.0/3.0* (1.0/36.0/3.14)**(-0.2043) * dvol * vold**(2.0*0.2043-1) * sp%Surface**(-0.2043*3.0+1)
                Case (OXIDATION_RADIUS)
                     sp%Surface = sp%Surface + (2.0 * dvol / OxidationRadius(vold, sp%Surface))
               !     sp%Surface = sp%Surface + 2.0/3.0 * dvol / vold * sp%Surface
                Case (EQUIV_SPHERE_DIAMETER)
                    sp%Surface = sp%Surface + (4.0 * dvol / EquivSphereDiameter(vold))
                Case (COLLISION_DIAMETER)
                    sp%Surface = sp%Surface + (4.0 * dvol / CollisionDiameter(vold, sp%Surface, mech))
            End Select
            sp%Surface = Max(sp%Surface, ssph)
        End If

        v1 = 10**(2.6) * 1.0e-21
        vc2 = 2.1443e-23
        sc2 = PI * ( 6.0*vc2/PI )**(2.0/3.0)
        vol=vold+dvol
        if(vol<v1) then
          theta = 2.0
        else
          theta = 3.0 * ( log(vol/v1) + 2.0/3.0 * log(v1/vc2) ) / log(vol/vc2)
        endif

        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
           theta=2.0
        End if
        surf = (vol/vc2)**(theta/3.0) * sc2
!*****************************************************************
! MODIFICATION SURFACE POUR IMPOSER LA LOI DE RODRIGUES ET AL 2018
!*****************************************************************
!        sp%Surface = surf
!        sp%Surface = Max(sp%Surface, ssph)


        ! Recalc properties.
		Call CalcProperties(sp, mech)
	End Subroutine

	! -------------------------------------------------------
	! PARTICLE PROPERTY ROUTINES.
	!
	!	The following routines calculate properties of
	!	stochastic particles that are functions of the unique
	!	particle coordinates.
	!
	! -------------------------------------------------------

	Function GetPreCalcs(sp)
		! DESCRIPTION:
		!	Gets pre-calculated properties defined for a
        !   stochastic particle.

		Implicit None

		! ARGUMENTS.
		Real :: GetPreCalcs(PROPERTY_COUNT)    ! Return value.
		Type(StochParticle), Intent(IN)	::	sp ! Particle for which to acquire properties.

		! EXECUTABLE CODE.
        GetPreCalcs(1:CACHE_COUNT) = sp%Properties
        GetPreCalcs(CACHE_COUNT+1) = sp%Surface
	End Function

	! -------------------------------------------------------

	Function CreateParticle(typeid, initComp, surf, createt, updatet, &
                            actsurf, track, mech) Result(sp)
		! DESCRIPTION:
		!	Creates a new particle with the given properties for
        !   the given mechanism.

        Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle) :: sp           ! Particle returned.
        Integer, Intent(IN) :: typeid       ! Particle type.
        Integer, Intent(IN) :: initComp(:)  ! Initial components
        Real, Intent(IN)    :: surf         ! Surface area.
        Real, Intent(IN)    :: createt      ! Create time.
        Real, Intent(IN)    :: updatet      ! Last update time.
        Real, Intent(IN)    :: actsurf      ! Active surface fraction.
        Real, Intent(IN)    :: track(:)     ! Tracking variables.
        Type(Mechanism), Intent(IN) :: mech ! Mechanism in use.

		! EXECUTABLE CODE.
		sp%Components = initComp
		sp%TypeID	  = typeid
        sp%CreateTime = createt
		sp%LastUpdate = updatet
		sp%ActSurf    = actsurf
		sp%Track      = track

        ! Surface area based on particle model.
        If (mech%ParticleModel == SPHERICAL_PARTICLE_MODEL) Then
            sp%Surface = EquivSphereSurface(sp, mech%Components, mech%ComponentCount)
        ElseIf (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
            sp%Surface = surf
        End If
        
		Call CalcProperties(sp, mech)
	End Function

	! -------------------------------------------------------

	Subroutine CalcProperties(sp, mech)
		! DESCRIPTION:
		!	Sets a particle's properties assuming that
		!	the unique particle coordinates have already
		!	been set.

		Use SWPMECH_TYPES
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(INOUT) :: sp ! Particle to update.
        Type(Mechanism), Intent(IN)        :: mech ! Mechanism in use.

		! EXECUTABLE CODE.
        If (Sum(sp%Components) > 0) Then
		    sp%Properties(iAS)		= sp%ActSurf * sp%Surface
		    sp%Properties(iS)		= sp%Surface
            sp%Properties(iD)       = CollisionDiameter(sp, mech)
		    sp%Properties(iD2)	    = sp%Properties(iD) * sp%Properties(iD)
		    sp%Properties(iD_1)	    = 1.0E0 / sp%Properties(iD)
		    sp%Properties(iD_2)		= 1.0E0 / sp%Properties(iD2)
            sp%Properties(iM)       = Mass(sp, mech%Components, mech%ComponentCount)
		    sp%Properties(iM_1_2)	= 1.0E0 / Sqrt(sp%Properties(iM))
		    sp%Properties(iD2M_1_2) = sp%Properties(iD2) * sp%Properties(iM_1_2)
            sp%Properties(iV)       = Volume(sp, mech%Components, mech%ComponentCount)
            sp%Properties(iSz)      = SizeP(sp)
        Else
		    sp%Properties = 0.0E0
        End If
	End Subroutine
	
	! -------------------------------------------------------
	! PARTICLE PROPERTIES.
	!
	!	These routines calculate additional properties of
    !   stochastic particles using the unique particle
    !   coordinates.
	!
	! -------------------------------------------------------

	Pure Integer Function SizeP_Part(sp)
		! DESCRIPTION:
		!	Returns the "size" of a stochastic particle.
		! RETURNS:
		!	"Size" of particle.
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        SizeP_Part = Sum(sp%Components)
	End Function

	! -------------------------------------------------------

	Pure Function SizeP_PartArray(sp) result (sizes)
		! DESCRIPTION:
		!	Returns the "sizes" of some stochastic particles.
		! RETURNS:
		!	"Sizes" of particles.
		Implicit None
		Type(StochParticle), Intent(IN) :: sp(:)
        Integer :: sizes(Size(sp)), i
        Do i = 1, Size(sp)
            sizes(i) = Sum(sp(i)%Components)
        End Do
	End Function

	! -------------------------------------------------------

	Pure Integer Function SizeP_Comp(comp)
		! DESCRIPTION:
		!	Returns the "size" of a stochastic particle.
		! RETURNS:
		!	"Size" of particle.
		Implicit None
        Integer, Intent(IN) :: comp(:)
        SizeP_Comp = Sum(comp)
	End Function

	! -------------------------------------------------------

	Pure Real Function Mass_Part(sp, comps, ncomps)
		! DESCRIPTION:
		!	Returns the mass of a stochastic particle.
		! RETURNS:
		!	Mass of particle (g).

        Use SWPMECH_TYPES
		Implicit None

        ! ARGUMENTS.
		Type(StochParticle), Intent(IN) :: sp
        Type(Component), Intent(IN)    :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
		
        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        Mass_Part = 0.0E0
        Do i = 1, ncomps
            Mass_Part = Mass_Part + (Real(sp%Components(i)) * comps(i)%MolWt)
        End Do
        Mass_Part = Mass_Part / NA
	End Function

	! -------------------------------------------------------

	Pure Real Function Mass_Comp(compValues, comps, ncomps)
		! DESCRIPTION:
		!	Returns the mass of a stochastic particle.
		! RETURNS:
		!	Mass of particle (g).

        Use SWPMECH_TYPES
		Implicit None

        ! ARGUMENTS.
		Integer, Intent(IN)         :: compValues(:)
        Type(Component), Intent(IN) :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
		
        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        Mass_Comp = 0.0E0
        Do i = 1, ncomps
            Mass_Comp = Mass_Comp + (Real(compValues(i)) * comps(i)%MolWt)
        End Do
        Mass_Comp = Mass_Comp / NA
	End Function

	! -------------------------------------------------------

	Pure Real Function Volume_Part(sp, comps, ncomps)
		! DESCRIPTION:
		!	Returns the volume of a stochastic particle.
		! RETURNS:
		!	Volume of particle (cm^3).

        Use SWPMECH_TYPES
		Implicit None

        ! ARGUMENTS.
		Type(StochParticle), Intent(IN) :: sp
        Type(Component), Intent(IN)     :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
		
        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        Volume_Part = 0.0E0
        Do i = 1, ncomps
            Volume_Part = Volume_Part + (Real(sp%Components(i)) * comps(i)%MolWt / comps(i)%Density)
        End Do
        Volume_Part = Volume_Part / NA
	End Function

	! -------------------------------------------------------

	Pure Real Function Volume_Comp(compValues, comps, ncomps)
		! DESCRIPTION:
		!	Returns the volume of a stochastic particle.
		! RETURNS:
		!	Volume of particle (cm^3).

        Use SWPMECH_TYPES
		Implicit None

        ! ARGUMENTS.
		Integer, Intent(IN)         :: compValues(:)
        Type(Component), Intent(IN) :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
		
        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        Volume_Comp = 0.0E0
        Do i = 1, ncomps
            Volume_Comp = Volume_Comp + (Real(compValues(i)) * comps(i)%MolWt / comps(i)%Density)
        End Do
        Volume_Comp = Volume_Comp / NA
	End Function

	! -------------------------------------------------------

	Pure Real Function SurfaceArea(sp)
		! DESCRIPTION:
		!	Returns the surface area of the given particle.
		! RETURNS:
		!	Surface area of particle (cm^2).
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
		SurfaceArea = sp%Surface
	End Function

	! -------------------------------------------------------

	Pure Real Function EquivSphereSurface_Part(sp, comps, ncomps)
		! DESCRIPTION:
		!	Calculates the equivalent sphere surface area
        !   of the given particle.
		! RETURNS:
		!	Equivalent sphere surface area (cm^2).
		Use SWPMECH_TYPES
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        Type(Component), Intent(IN)     :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
        EquivSphereSurface_Part = EquivSphereSurface_Vol(Volume(sp, comps, ncomps))
	End Function

	! -------------------------------------------------------

	Pure Real Function EquivSphereSurface_Comp(compValues, comps, ncomps)
		! DESCRIPTION:
		!	Calculates the equivalent sphere surface area
        !   of the given particle.
		! RETURNS:
		!	Equivalent sphere surface area (cm^2).
		Use SWPMECH_TYPES
		Implicit None
		Integer, Intent(IN)          :: compValues(:)
        Type(Component), Intent(IN)  :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
        EquivSphereSurface_Comp = EquivSphereSurface_Vol(Volume(compValues, comps, ncomps))
	End Function

	! -------------------------------------------------------

	Pure Real Function EquivSphereSurface_Vol(vol)
		! DESCRIPTION:
		!	Calculates the equivalent sphere surface area
        !   for the given volume.
		! RETURNS:
		!	Equivalent sphere surface area (cm^2).
		Use SWPPARAMS, only: PI, TWO_THIRDS
		Implicit None
		Real, Intent(IN) :: vol
		EquivSphereSurface_Vol = PI * (6.0*vol/PI)**TWO_THIRDS
	End Function

	! -------------------------------------------------------

	Pure Real Function EquivSphereDiameter_Part(sp, comps, ncomps)
		! DESCRIPTION:
		!	Calculates the equivalent sphere diameter
        !   of the given particle.
		! RETURNS:
		!	Equivalent sphere diameter (cm).
		Use SWPMECH_TYPES
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        Type(Component), Intent(IN)     :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
		EquivSphereDiameter_Part = (6.0 * Volume(sp, comps, ncomps) / PI) ** ONE_THIRD
	End Function

	! -------------------------------------------------------

	Pure Real Function EquivSphereDiameter_Vol(vol)
		! DESCRIPTION:
		!	Calculates the equivalent sphere diameter
        !   of the given particle.
		! RETURNS:
		!	Equivalent sphere diameter (cm).
		Use SWPPARAMS
		Implicit None
		Real, Intent(IN) :: vol ! cm^3.
		EquivSphereDiameter_Vol = (6.0 * vol / PI) ** ONE_THIRD
	End Function

	! -------------------------------------------------------

	Pure Real Function CollisionDiameter_Part(sp, mech)
		! DESCRIPTION:
		!	Calculates the collision diameter of the particle.
		! RETURNS:
		!	Collision diameter of particle.
		Use SWPMECH_TYPES
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        Type(Mechanism), Intent(IN) :: mech
        REAL vol

        CollisionDiameter_Part = (6.0 * Volume(sp, mech%Components, mech%ComponentCount) / PI) ** ONE_THIRD
        If (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
		!    CollisionDiameter_Part = (CollisionDiameter_Part + Sqrt(sp%Surface / PI)) * ONE_HALF
!********************************************************************************************
! DEFINITION DU DIAMETRE COLLISIONEL CONSISTENT AVEC LE MODEL SECTIONEL ==> PAS
! FORCEMENT GENERALE
!****************************************************************************************
                vol = Volume(sp, mech%Components, mech%ComponentCount)
                CollisionDiameter_Part = 6.0/(36.0*PI)**(1./1.8) &
                                  * vol**(1.0-2.0/1.8)   &
                                  * sp%Surface**(3.0/1.8 - 1.0)
        End If

	End Function

	! -------------------------------------------------------

	Pure Real Function CollisionDiameter_SurfVol(v, s, mech)
		! DESCRIPTION:
		!	Calculates the collision diameter of the particle.
		! RETURNS:
		!	Collision diameter of particle.
		Use SWPPARAMS
        Use SWPMECH_TYPES
		Implicit None
		Real, Intent(IN) :: v, s !v=cm^3, s=cm^2.
        Type(Mechanism), Intent(IN) :: mech

        CollisionDiameter_SurfVol = (6.0 * v / PI) ** ONE_THIRD
        If (mech%ParticleModel == SURFACE_VOLUME_MODEL) Then
!		    CollisionDiameter_SurfVol = (CollisionDiameter_SurfVol + Sqrt(s / PI)) * ONE_HALF
                CollisionDiameter_SurfVol = 6.0/(36.0*PI)**(1./1.8) &
                                  * v**(1.0-2.0/1.8)   &
                                  * s**(3.0/1.8 - 1.0)
        End If


	End Function

	! -------------------------------------------------------

	Pure Real Function GrowthRadius_Part(sp)
		! DESCRIPTION:
		!	Radius used to calculate change in surface area
        !   for a growth surface process in the surface-volume
        !   model.
		! RETURNS:
		!	Growth radius (cm).
		Use SWPPARAMS, only: PI
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        GrowthRadius_Part = Sqrt(sp%Surface / (4*PI))
	End Function

	! -------------------------------------------------------

	Pure Real Function GrowthRadius_Surf(surf)
		! DESCRIPTION:
		!	Radius used to calculate change in surface area
        !   for a growth surface process in the surface-volume
        !   model.
		! RETURNS:
		!	Growth radius (cm).
		Use SWPPARAMS, only: PI
		Implicit None
		Real, Intent(IN) :: surf
        GrowthRadius_Surf = Sqrt(surf / (4*PI))
	End Function

	! -------------------------------------------------------

	Pure Real Function OxidationRadius_Part(sp, comps, ncomps)
		! DESCRIPTION:
		!	Radius used to calculate change in surface area
        !   for an oxidationj surface process in the
        !   surface-volume model.
		! RETURNS:
		!	Oxidation radius.
		Use SWPMECH_TYPES
		Implicit None
		Type(StochParticle), Intent(IN) :: sp
        Type(Component), Intent(IN)     :: comps(:)
        Integer, Intent(IN) :: ncomps ! Number of components in mechanism.
        OxidationRadius_Part = 3.0 * Volume(sp, comps, ncomps) / sp%Surface
	End Function

	! -------------------------------------------------------

	Pure Real Function OxidationRadius_SurfVol(v, s)
		! DESCRIPTION:
		!	Radius used to calculate change in surface area
        !   for an oxidationj surface process in the
        !   surface-volume model.
		! RETURNS:
		!	Oxidation radius.
		Implicit None
		Real, Intent(IN) :: v, s ! v=cm^3, s=cm^2.
        OxidationRadius_SurfVol = 3.0 * v / s
	End Function

	! -------------------------------------------------------
	! PROCESS RATE ROUTINES.
	!
	!	These routines are used to calculate particle
    !   properties required for process rate calculations.
	!
	! -------------------------------------------------------

	Pure Real Function CalcParticleWeight(sp, masswt, surfwt, actsurfwt, diamwt)
		! DESCRIPTION:
		!	Calculates a custom particle weight given the
		!	exponents of particle properties.
		!
		! RETURNS:
		!	Custom weight.

		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(IN)	::	sp	! The particle to calculate the weight for.
        Real, Intent(IN) :: masswt, &    ! Mass/volume exponent.
                            surfwt, &    ! Surface area exponent.
                            actsurfwt, & ! Active surface exponent.
                            diamwt       ! Diameter exponent.

		! EXECUTABLE CODE.
		CalcParticleWeight = 1.0E0
		CalcParticleWeight = CalcParticleWeight * sp%Surface**surfwt
		CalcParticleWeight = CalcParticleWeight * sp%Properties(iAS)**actsurfwt
		CalcParticleWeight = CalcParticleWeight * sp%Properties(iD)**diamwt
		CalcParticleWeight = CalcParticleWeight * sp%Properties(iM)**masswt
	End Function

	! -------------------------------------------------------

	Function GetPropertyIndex(p) Result(ix)
		! DESCRIPTION:
        !   Takes a list of property exponents in the order
        !   Mass, Surface, Active Surface, and Diameter. Returns
        !   the index of a particle property which is described
        !   by this exponent list.  If no particle property
        !   exists returns -1.

        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
        Real             :: ix   ! Return value.
		Real, Intent(IN) :: p(4) ! List of property powers.

		! VARIABLES.
		Logical	::	test

		! EXECUTABLE CODE.

		! Assume reaction has a custom dependancy, then
		! check built-in weightings to see if we are wrong.
		ix = -1

		! ***********************************
		!  BUILT-IN VARIABLE COMBINATIONS.
		! ***********************************

		! M^-1/2 * D^2.
		test  =	(p(1) == -ONE_HALF) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 2.0E0)

		If (test) Then
			ix = iD2M_1_2
			Return
		End If

		! ***********************************
		! BUILT-IN SINGLE VARIABLE.
		! ***********************************

		! S^1.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 1.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 0.0E0)

		If (test) Then
			ix = iS
			Return
		End If

		! AS^1.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 1.0E0) .And. &
				(p(4) == 0.0E0)

		If (test) Then
			ix = iAS
			Return
		End If

		! D^1.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 1.0E0)

		If (test) Then
			ix = iD
			Return
		End If

		! D^-1.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == -1.0E0)

		If (test) Then
			ix = iD_1
			Return
		End If

		! D^2.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 2.0E0)

		If (test) Then
			ix = iD2
			Return
		End If

		! D^-2.
		test  =	(p(1) == 0.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == -2.0E0)

		If (test) Then
			ix = iD_2
			Return
		End If

		! M^1.
		test  =	(p(1) == 1.0E0) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 0.0E0)

		If (test) Then
			ix = iM
			Return
		End If

		! M^-1/2.
		test  =	(p(1) == -ONE_HALF) .And. &
				(p(2) == 0.0E0) .And. &
				(p(3) == 0.0E0) .And. &
				(p(4) == 0.0E0)

		If (test) Then
			ix = iM_1_2
			Return
		End If
	End Function

	! -------------------------------------------------------
	! STATISTICAL ITEMS ROUTINES.
	!
	!	The following routines calculate statistically
	!	useful properties for stochastic particles.
	!
	! -------------------------------------------------------

	Pure Function GetStats(sp) result (stats)
		! DESCRIPTION:
		!	Gets an array of statistically useful
		!	properties about a particle.
		! RETURNS:
		!	Particle statistics.

		Use SWPPART_STATS
		Implicit None

		! ARGUMENTS.
		Type(StochParticle), Intent(IN) :: sp
		Double Precision :: stats(Count)

		! EXECUTABLE CODE.

		! Stat items are:
		!	Particle count (M0), 
		!	M1 - M6,
		!	Volume (cm3), 
		!	Mass (g),
		!	Surface area (cm2),
		!	Active surface fraction,
		!	Diameter (cm).
		stats(iM0)	    = 1.0D0
		stats(iM1)	    = Dble(Sum(sp%Components))
		stats(iM2)	    = stats(iM1) * stats(iM1)
		stats(iM3)	    = stats(iM1) * stats(iM2)
		stats(iM4)	    = stats(iM1) * stats(iM3)
		stats(iM5)	    = stats(iM1) * stats(iM4)
		stats(iM6)	    = stats(iM1) * stats(iM5)
        stats(iMass)    = Dble(sp%Properties(iM))
        stats(iV)       = Dble(sp%Properties(iV))
        stats(iM2)      = stats(iV) * stats(iV) 
		stats(iSurf)    = Dble(sp%Surface)
!		stats(iActSurf) = Dble(sp%ActSurf)
		stats(iActSurf) = Dble(sp%Surface)**(3.0D0/2.0D0)
		stats(iDiam)    = Dble(sp%Properties(iD))
	End Function

End Module

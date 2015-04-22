! *****************************************************************************
!
! File:					swpmech_types.f90
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
!	Definition of types which describe the particle mechanism in Sweep.
! *****************************************************************************

Module SWPMECH_TYPES
    Use SWPPARAMS
	Implicit None
	Public

 	! -------------------------------------------------------
	! MECHANISM PARAMETERS.
	! -------------------------------------------------------
 
    ! Number of additional processes above Inceptions and Reactions
    ! solved by Sweep.  This is the number of coagulation terms +
    ! inflow and outflow.
    Integer, Parameter :: EXTRA_PROCESS_COUNT = NCOAGTERMS + 2
    Integer, Parameter :: EXTRA_GROUP_COUNT = 3

	! -------------------------------------------------------
	! PARTICLE COMPOSITION PARAMETERS.
	! -------------------------------------------------------

    Integer, Parameter :: MAX_COMP  = 1 ! Max. number of particle components.
	Integer, Parameter :: MAX_TRACK = 1 ! Max. tracking variable count.

	! -------------------------------------------------------
	! REACTION DEFINITION PARAMETERS.
	! -------------------------------------------------------

	! Max. number of species in a surface/inception reaction.
	Integer, Parameter, Private	::	NRXNSP = 6

	! Values of Reaction%ID other than particle property indices.
	Integer, Parameter	::	CUSTOM_ID = -2
	Integer, Parameter	::	INCEPT_ID = -1
	Integer, Parameter	::	UNIFORM_ID = 0

	! -------------------------------------------------------
	! INCEPTION RXN DEF'N.
	! -------------------------------------------------------

	Type Inception
		Integer	:: Species(NRXNSP)   = 0   ! Involved species (reactants + products).
		Integer	:: NREAC, NPROD	     = 0   ! Number of reactants and products.
		Integer	:: Stoich(NRXNSP)    = 0   ! Stoichiometric coefficients of involved species.
		Integer	:: PType             = 0   ! Type of the incepted soot particle.
		Real	:: A                 = 1.0 ! Inception rate constant.
		Real	:: kfm, ksf1, ksf2   = 0.0 ! Kernel parameters for free-mol and slip-flow.
        Real    :: dS                = 0.0 ! Surface area of incepted particles.
        Integer :: dComp(MAX_COMP)   = 0   ! Initial component values for incepted particle.
        Real    :: dTrack(MAX_TRACK) = 0   ! Initial tracker variable values.
	End Type

	! -------------------------------------------------------
	! SURFACE RXN DEF'N.
	! -------------------------------------------------------

	Type Reaction
		Integer	:: Species(NRXNSP)   = 0   ! Involved species (reactants + products).
		Integer	:: NREAC, NPROD	     = 0   ! Number of reactants and products.
		Integer	:: Stoich(NRXNSP)    = 0   ! Stoichiometric coefficients of involved species.
		Integer	:: PTypeIn		     = 0   ! Input and output types of the reacting soot particle.
		Integer	:: PTypeOut		     = 0
		Real	:: A    		     = 1.0 ! Arrhenius rate constant.
		Real	:: n, E 		     = 0.0 ! Arrhenius rate constants.
        Real    :: P(4)              = 0.0 ! Custom property powers.
		Integer	:: ID			     = 0   ! Particle weight ID (-2=Custom, -1=Inception, 0=Uniform).
        Integer :: dSID              = 0   ! Index of the radius to use for particle surface area change.
        Integer :: dComp(MAX_COMP)   = 0   ! Component value changes.
        Real    :: dTrack(MAX_TRACK) = 0   ! Tracker variable changes.
	End Type

	! -------------------------------------------------------
	! PARTICLE COMPONENT DEF'N.
	! -------------------------------------------------------

    Type Component
        Character(LEN=16) :: Symbol  = ""    ! Symbolic representation.
        Real              :: Density = 0.0E0 ! Density (mol/cm3).
        Real              :: MolWt   = 0.0E0 ! Molecular weight (g/mol).
    End Type

	! -------------------------------------------------------
	! COAGULATION RULES.
	! -------------------------------------------------------

    Type CoagRule
		Integer	:: In1=-1, In2=-1 ! Input particle types.
		Integer	:: Out=-1         ! Output particle type.
    End Type

	! -------------------------------------------------------
	! MECHANISM DEF'N.
	! -------------------------------------------------------

    Type Mechanism
        ! Total process count incl. coagulation.
        Integer :: ProcessCount = 0
        ! Inception processes.
        Integer                  :: InceptionCount = 0
        Type(Inception), Pointer :: Inceptions(:)
        ! Surface reactions incl. condensations.
        Integer                  :: ReactionCount = 0
        Type(Reaction), Pointer  :: Reactions(:)
        ! Mask of deferred processes.
        Logical, Pointer  :: DeferMask(:)
        Logical           :: AnyDeferred = .False.
        ! Reaction groupings.
        Integer                    :: GroupCount = 0
        Integer, Pointer           :: Groups(:)
        Character(LEN=50), Pointer :: GroupNames(:)
        ! Particle type and composition information.
        Integer           :: NParticleTypes = 0
        Character(LEN=16), Pointer :: ParticleTypeNames(:)
        Integer           :: ComponentCount
        Type(Component)   :: Components(MAX_COMP)
        Integer           :: TrackerCount = 0
        Character(LEN=16) :: Trackers(MAX_TRACK)
        ! Helpful indices in rate arrays.
        Integer :: ISR=0, ICG=0, IIN=0, IOUT=0
        ! Flags to switch on/off different processes.
        Logical :: SurfOn = .True.  ! Surface reactions.
        Logical :: CoagOn = .True.  ! Coagulation.
        Logical :: SintOn = .False. ! Sintering.
        ! Coagulation rules.
        Integer :: CoagRuleCount = 0
        Type(CoagRule), Pointer :: CoagRules(:)
        ! Sintering parameters.
        Real :: SintParams(2) = 0.0E0
        ! Active surface model.
        Logical :: UseHACA = .False.
        Integer :: ActSurfModel = ACT_SURF_CONST
        Real    :: ConstActSurf = 1.0E0
        ! Particle model.
        Integer :: ParticleModel = SPHERICAL_PARTICLE_MODEL
    End Type
End Module

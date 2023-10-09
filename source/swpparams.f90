! *****************************************************************************
!
! File:				swpparams.f90
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
!   A group of constants which define properties of soot particles.  Useful
!   for calculating ensemble properties and reaction rates.
!
!   File also contains control parameters for Sweep and definitions of model
!   constants.
!
! *****************************************************************************

Module SWPPARAMS
    Implicit None
    Public

    ! -------------------------------------------------------
    ! MECHANISM CONTROL.
    ! -------------------------------------------------------

    ! Maximum mechanism size.
    Integer, Parameter :: MAXICNS = 100 ! Max number of inceptions.
    Integer, Parameter :: MAXRXNS = 300 ! Max number of surface processes.
    Integer, Parameter :: MAXGRPS = 403 ! Max number of reaction groups.

    ! Flags to switch processes on and off.
    Logical, Parameter :: SURF_ON = .True.  ! Surface reactions.
    Logical, Parameter :: COAG_ON = .True.  ! Coagulation.
    Logical, Parameter :: SINT_ON = .False. ! Sintering.
    Logical, Parameter :: OBLT_ON = .True.  ! Obliteration.

    ! Majorant factor for surface reactions.
    Real, Parameter    :: SURF_MAJ = 2.0E0

    ! Scaling factor for sintering.  The minimum change
    ! in particle surface is SINT_SCALE * Ssph, where Ssph
    ! is the surface area of an equivalent volume sphere, i.e.
    ! the theoretical minimum surface area.
    Real, Parameter    :: SINT_SCALE = 0.1E0

    ! Turn on for an extra level of resolution when updating
    ! deferred particle events.
    Logical, Parameter :: LPDA_TWO_LEVELS = .False.

    ! -------------------------------------------------------
    ! COAGULATION PARAMETERS.
    ! -------------------------------------------------------

    ! Number of coagulation processes there are in the
    ! particle dynamics model.
    Integer, Parameter, Public :: NCOAGTERMS = 7
    Integer, Parameter :: COAG_TYPE = 3     ! (1) Free-molecular
                                            ! (2) Continuum      
                                            ! (3) Harmonic mean

    ! -------------------------------------------------------
    ! PARTICLE MODEL.
    ! -------------------------------------------------------

    ! Enumeration of different particle models.
    Integer, Parameter :: SPHERICAL_PARTICLE_MODEL = 1
    Integer, Parameter :: SURFACE_VOLUME_MODEL     = 2
    Integer, Parameter :: FRACTAL_MODEL            = 3

    ! Minimum particle size (number of carbon atoms).
    Integer, Parameter :: MIN_PARTICLE_SIZE = 2 !

    ! -------------------------------------------------------
    ! ACTIVE SURFACE MODEL.
    ! -------------------------------------------------------

    ! Enumeration of active surface models.
    Integer, Parameter :: ACT_SURF_CONST     = 1 ! Alpha kept constant at CONST_ALPHA.
    Integer, Parameter :: ACT_SURF_PARTICLE  = 2 ! Alpha given individually as a particle property.
    Integer, Parameter :: ACT_SURF_ABF       = 3 ! Alpha given by Frenklach correlation.
    Integer, Parameter :: ACT_SURF_PROFILE   = 4 ! Alpha given as a parameter in the chemistry profile.

    ! -------------------------------------------------------
    ! SURFACE AGING.
    ! -------------------------------------------------------

    Logical, Parameter :: AGE_SURFACE      = .False. ! On/off.
    Real, Parameter    :: BASE_ACT_SURF    = 0.01E0
    Real, Parameter    :: INITIAL_ACT_SURF = 1.0E0

    ! -------------------------------------------------------
    ! PARTICLE DIAMETERS ENUMERATION.
    ! -------------------------------------------------------

    ! Enumeration of different particle radii.
    Integer, Parameter :: EQUIV_SPHERE_DIAMETER = 0
    Integer, Parameter :: OXIDATION_RADIUS      = 1
    Integer, Parameter :: GROWTH_RADIUS         = 2
    Integer, Parameter :: COLLISION_DIAMETER    = 3

    ! -------------------------------------------------------
    ! SHARED CHEMISTRY PARAMETERS.
    ! -------------------------------------------------------

    ! Enurmation of types of chemistry used in Sweep.  Fixed
    ! chemistry is when the chemical conditions are not
    ! affected by the soot.  Variable chem is when the chemical
    ! conditions are affected by the soot.
    Integer, Parameter :: FIXED_CHEM=1, VARIABLE_CHEM=2

    ! -------------------------------------------------------
    ! OUTFLOW TREATMENT.
    ! -------------------------------------------------------

    ! Enumeration of different outflow methods.
    Integer, Parameter :: OUTFLOW_SCALE_METHOD    = 0
    Integer, Parameter :: OUTFLOW_PARTICLE_METHOD = 1
    Integer, Parameter :: OUTFLOW_METHOD = OUTFLOW_PARTICLE_METHOD

    ! -------------------------------------------------------
    ! STEP CONTROL.
    ! -------------------------------------------------------

    ! Set to true to reset counters every time Sweep is run.
    Logical, Parameter :: RESET_COUNTERS = .True.

    ! -------------------------------------------------------
    ! GENERAL CONTROL.
    ! -------------------------------------------------------

    ! Set to true to track computation time.
    Logical, Parameter :: USE_CT = .False.

    ! Flag to turn on/off console messages.
    Logical, Parameter :: MSGS = .True.
    
    ! -------------------------------------------------------
    ! PHYSICAL CONSTANTS.
    ! -------------------------------------------------------

    ! Molecular weight of Carbon.
    Real, Parameter :: CMW = 12.011E0 ! g/mol

    ! Average density of soot (density of graphite).
    Real, Parameter :: RHOS = 1.860E0 ! g/cm3

    ! Avogadro's number.
    Real, Parameter :: NA = 6.0221415E23

    ! PI.
    Real, Parameter :: PI = 3.141592653589793E0

    ! Boltzmann's constant.
    Real, Parameter :: KB = 1.3806488E-16 ! gcm2/s2K

    ! Site density on soot surface.
    Real, Parameter :: CHI = 2.32E15 ! 1/cm2

    ! Ideal gas constant.
    Real, Parameter :: R = 8.31434E0  ! J/molK
    Real, Parameter :: RCAL = 1.9872E-3  ! kcal/molK

    ! Constant term in calculation on Knudsen number.  Actual
    ! Kn = K * T [K] / (P [bar] * d [cm]).
    Real, Parameter :: KNUDSEN_K = 4.74151636E-8

    ! Fractal coefficient
    Real, Parameter :: Dfrac = 1.8E0
    
    ! -------------------------------------------------------
    ! USEFUL QUANTITIES.
    ! -------------------------------------------------------

    Real, Parameter     ::  ONE_THIRD       = 3.333333333333334E-01
    Real, Parameter     ::  TWO_THIRDS      = 6.666666666666667E-01
    Real, Parameter     ::  ONE_SIXTH       = 1.666666666666667E-01
    Real, Parameter     ::  ONE_HALF        = 5.000000000000000E-01
    Real, Parameter     ::  ROOT_TWO        = 1.414213562373095E+00

    ! -------------------------------------------------------
    ! PARTICLE SIZE CONVERSIONS.
    ! -------------------------------------------------------

    ! Conversion from number of carbon atoms, N.
    !  NB:  Conversion to D requires N^1/3 and conversion
    !       to S requires N^2/3.
    Real, Parameter :: NTOM    = 1.9945199E-23 ! g
    Real, Parameter :: NTOV    = 1.1080667E-23 ! cm3
    Real, Parameter :: NTOD    = 2.7660214E-08 ! cm
    Real, Parameter :: NTOS    = 2.4035932E-15 ! cm2

    ! Conversion from mass, M (g).
    !  NB:  Conversion to D requires M^1/3 and conversion
    !       to S requires M^2/3.
    Real, Parameter :: MTON    = 5.0137377E+22
    Real, Parameter :: MTOV    = 5.5555558E-01 ! cm3
    Real, Parameter :: MTOD    = 1.0610330E+00 ! cm
    Real, Parameter :: MTOS    = 3.5367770E+00 ! cm2

    ! Conversion from volume, V (cm3).
    !  NB:  Conversion to D requires V^1/3 and conversion
    !       to S requires V^2/3.
    Real, Parameter :: VTON    = 9.0247273E+22
    Real, Parameter :: VTOM    = 1.8000000E+00 ! g
    Real, Parameter :: VTOD    = 1.2407010E+00 ! cm
    Real, Parameter :: VTOS    = 4.8359756E+00 ! cm2

    ! Conversion from diameter, D (cm).
    !  NB:  Conversion to N, M, or V requires S^3/2 and conversion
    !       to D requires S^1/2.
    Real, Parameter :: DTON   = 4.7253366E+22
    Real, Parameter :: DTOM   = 9.4247776E-01 ! g
    Real, Parameter :: DTOV   = 5.2359879E-01 ! cm3
    Real, Parameter :: DTOS   = PI ! cm2

    ! Conversion from surface, S (cm2).
    !  NB:  Conversion to N, M, or V requires D^3 and conversion
    !       to S requires D^2.
    Real, Parameter :: STON    = 8.4861080E+21
    Real, Parameter :: STOM    = 1.5034483E-01 ! g
    Real, Parameter :: STOV    = 9.4031602E-02 ! cm3
    Real, Parameter :: STOD    = 5.6418955E-01 ! cm

    ! -------------------------------------------------------
    ! UNIT CONVERSIONS.
    ! -------------------------------------------------------

    ! Engineering Magnitudes (conversion to..)
    Real, Parameter :: MEGA    = 1.0E-6
    Real, Parameter :: KILO    = 1.0E-3
    Real, Parameter :: CENTI   = 1.0E2
    Real, Parameter :: MILLI   = 1.0E3
    Real, Parameter :: MICRO   = 1.0E6
    Real, Parameter :: NANO    = 1.0E9

    ! .. squares,
    Real, Parameter :: MEGA2   = 1.0E-12
    Real, Parameter :: KILO2   = 1.0E-06
    Real, Parameter :: CENTI2  = 1.0E+04
    Real, Parameter :: MILLI2  = 1.0E+06
    Real, Parameter :: MICRO2  = 1.0E+12
    Real, Parameter :: NANO2   = 1.0E+18
 
    ! .. cubes.
    Real, Parameter :: MEGA3   = 1.0E-18
    Real, Parameter :: KILO3   = 1.0E-09
    Real, Parameter :: CENTI3  = 1.0E+06
    Real, Parameter :: MILLI3  = 1.0E+09
    Real, Parameter :: MICRO3  = 1.0E+18
    Real, Parameter :: NANO3   = 1.0E+27

    ! Custom conversions.
    Real, Parameter :: CENTI_TO_NANO  = 1.0E+07
    Real, Parameter :: NANO_TO_CENTI  = 1.0E-07
End Module

! *****************************************************************************
!
! File:				swpstats.f90
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
!   Helper routines to get useful statistics from the particle ensemble in
!   Sweep2.
!
! Functions:
!   GetStats		-   Returns an array of statistical properties
!			    for the whole particle ensemble.
!   GetStatsType        -   Returns an array of statistical properties for the
!                           whole ensemble, grouped by particle type.
!   SumStats            -   Adds up all the particle stats for the ensemble.
!   SumStats            -   Adds up all the particle stats for the ensemble
!                           grouped by particle type.
!   CompilePSL          -   Generates a 4D particle size list for the ensemble.
!   CompilePCL          -   Generates particle composition list for the given
!                           ensemble.
! *****************************************************************************

Module SWPSTATS_INDICES
    Implicit None
    Public

    ! -------------------------------------------------------
    ! INDICES IN ARRAYS.
    ! -------------------------------------------------------

    Integer, Parameter :: Count = 29
    
    Integer, Parameter :: iSP=1,iNs=2,iFV=3,iSurf=4,iPrimPar=5
    Integer, Parameter :: iAvgFv=6,iFvSDv=7,iAvgSurf=8,iSurfSDv=9,iAvgPrimPar=10,iPrimParSDv=11
    Integer, Parameter :: iDiam=12,iDiamSDv=13,iDiamSk=14,iDiamKt=15,iSV=16,iScale=17
    Integer, Parameter :: iM0=18,iM1=19,iM2=20,iM3=21,iM4=22,iM5=23,iM6=24,iM7=25,iM8=26,iM9=27,iP1=28,iP2=29
 
    ! -------------------------------------------------------
    ! NAMES.
    ! -------------------------------------------------------
   
    Character(LEN=*), Parameter :: Names(Count) = (/"Stoch. Particle Count",&
        "Number Density (cm-3)","Volume Fraction", &
        "Surface Area (cm2/cm3)","Primary particle number (cm-3)",&
        "Average Volume (cm3)","Volume Stand. Deviation (cm3)",&
        "Average Surface Area (cm2)","Surface Area Stand. Deviation (cm2)",&
        "Average Primary Particle Number","Primary Particle Number Stand. Deviation",&
        "Average Primary Particle Diameter (nm)","Primary Particle Diameter Stand. Deviation (nm)",&
        "Primary Particle Diameter Skewness (nm^3)","Primary Particle Diameter kurtosis (nm^4)",&
        "Sample Volume (cm3)","Ensemble Scaling Factor",&
        "M0/3","M1/3","M2/3","M3/3","M4/3","M5/3","M6/3","M7/3","M8/3","M9/3","P1","P2"/)
End Module

! ==============================================================================

Module SWPSTATS
    Use SWPSTATS_INDICES, only: NSTATS=>Count, StatNames=>Names
    Implicit None
    Public
    Private :: SumStats, SumStatsTyped

    Integer, Public, Parameter :: PSL_VAR_COUNT = 5

    ! Upper and lower values of stat
    ! property bounds.  Only particles whose
    ! given property lies within these bounds
    ! shall be included in the statistics. 
    Type StatRange
        Double Precision :: Lower  = 0.0     ! Lower property bound.
        Double Precision :: Upper  = 1.0D300 ! Upper property bound.
        Integer          :: PropID = 2       ! Property used to calculated bounds.
    End Type

    Contains

    ! -------------------------------------------------------

    Function GetStats(soln, T, range, typeID) result (stats)
        ! DESCRIPTION:
        !   Returns an array of statistical properties of
        !   the current ensemble.
        ! RETURNS:
        !   Current ensemble statistics.

        Use SWPPARAMS
        Use SWPSTATS_INDICES
        Use SWPSOLN
        Use SWPPART_STATS, ParticleStatCount=>Count, iPM0=>iM0, iPM1=>iM1, &
                           iPM2=>iM2, iPM3=>iM3, iPM4=>iM4, iPM5=>iM5, &
                           iPM6=>iM6, iPM7=>iM7, iPM8=>iM8, iPM9=>iM9, iPP1=>iP1, iPP2=>iP2, &
                           iPV=>iVol, iPV2=>iVol2, iPSurf=>iSurf, iPSurf2=>iSurf2, &
                           iPDiam=>iDiam, iPDiam2=>iDiam2, iPDiam3=>iDiam3, iPDiam4=>iDiam4, & 
                           iPPrimPar=>iPrimPar, iPPrimPar2=>iPrimPar2
        Implicit None

        ! ARGUMENTS.
        Double Precision                :: stats(Count) ! Return value.
        Type(Solution), Intent(IN)      :: soln         ! Solution for which to calculate stats.
        Real, Intent(IN)                :: T            ! Temperature at which to calculate stats (K).
        Type(StatRange), Intent(IN)     :: range        ! Low and High cut-off values of the given particle
                                                        ! property.  particles outside this range will
                                                        ! not be included in the summation.
        Integer, Optional, Intent(IN)   :: typeID       ! Optional type of particle for which to 
                                                        ! get statistics.  If omitted gets statistics
                                                        ! for all particles.

        ! VARIABLES.
        Real :: v
        Double Precision :: s(ParticleStatCount)

        ! EXECUTABLE CODE.

        ! Pre-calc sample volume so we can scale the
        ! stats to a unit volume.
        v = SampleVolume(soln, T)

        ! Get statistics from the ensemble.
        If (Present(typeID)) Then
            s = SumStats(soln%Ensemble, range, typeID)
        Else
            s = SumStats(soln%Ensemble, range)
        End If

        ! Convert to statistics we want to export.
        stats(iSP)         = s(iPM0)
        stats(iM0)         = s(iPM0)
        stats(iM1)         = s(iPM1)
        stats(iM2)         = s(iPM2)
        stats(iM3)         = s(iPM3)
        stats(iM4)         = s(iPM4)
        stats(iM5)         = s(iPM5)
        stats(iM6)         = s(iPM6)
        stats(iM7)         = s(iPM7)
        stats(iM8)         = s(iPM8)
        stats(iM9)         = s(iPM9)
        stats(iP1)         = s(iPP1)
        stats(iP2)         = s(iPP2)
        stats(iNs)         = s(iPM0)
        stats(iFV)         = s(iPV)
        stats(iSurf)       = s(iPSurf)
        stats(iPrimPar)    = s(iPPrimPar)
        stats(iSV)         = Dble(v)
        stats(iScale)      = Dble(soln%LastKnownScale)

        ! Averaged properties.
        If (stats(iSP) > 0) Then
            stats(iAvgFV)      = stats(iFV) / stats(iSP)
            stats(iFvSDv)      = sqrt(max(0.0D0, s(iPV2)/stats(iM0) - stats(iAvgFV)**2.D0))
            stats(iAvgSurf)    = stats(iSurf) / stats(iSP)
            stats(iSurfSDv)    = sqrt(max(0.0D0, (s(iPsurf2) / stats(iSP)) - stats(iAvgSurf)**2.0D0 ))
            stats(iDiam)       = s(iPDiam) / stats(iSP)
            stats(iDiamSDv)    = sqrt(max(0.0D0, (s(iPDiam2) / stats(iSP)) -  stats(iDiam)**2.0D0 ))
            stats(iDiamSk)     = max(0.0D0, (s(iPDiam3) / stats(iSP)) - & 
                                    3.0D0*stats(iDiam)*stats(iDiamSDv)**2.0D0 - & 
                                          stats(iDiam)**3.0D0)
            stats(iDiamKt)     = max(0.0D0, (s(iPDiam4) / stats(iSP)) - & 
                                    4.0D0*stats(iDiam)*stats(iDiamSk) - & 
                                    6.0D0*stats(iDiam)**2.0D0*stats(iDiamSDv)**2.0D0 - & 
                                          stats(iDiam)**4.0D0)
            stats(iAvgPrimPar) = stats(iPrimPar) / stats(iSP)
            stats(iPrimParSDv) = sqrt(max(0.0D0, s(iPPrimPar2)/stats(iM0) - stats(iAvgPrimPar)**2.D0))
        Else
            stats(iAvgFV)      = 0.0D0
            stats(iFvSDv)      = 0.0D0
            stats(iAvgSurf)    = 0.0D0
            stats(iSurfSDv)    = 0.0D0
            stats(iDiam)       = 0.0D0
            stats(iDiamSDv)    = 0.0D0
            stats(iDiamSk)     = 0.0D0
            stats(iDiamKt)     = 0.0D0
            stats(iAvgPrimPar) = 1.0D0
            stats(iPrimParSDv) = 0.0D0
        End If

        ! Scale stats to unit volume.
        stats(iM0)   = stats(iM0)   / stats(iSV)
        stats(iM1)   = stats(iM1)   / stats(iSV)
        stats(iM2)   = stats(iM2)   / stats(iSV)
        stats(iM3)   = stats(iM3)   / stats(iSV)
        stats(iM4)   = stats(iM4)   / stats(iSV)
        stats(iM5)   = stats(iM5)   / stats(iSV)
        stats(iM6)   = stats(iM6)   / stats(iSV)
        stats(iM7)   = stats(iM7)   / stats(iSV)
        stats(iM8)   = stats(iM8)   / stats(iSV)
        stats(iM9)   = stats(iM9)   / stats(iSV)
        stats(iP1)   = stats(iP1)   / stats(iSV)
        stats(iP2)   = stats(iP2)   / stats(iSV)
        stats(iNs)   = stats(iNs)   / stats(iSV)
        stats(iFV)   = stats(iFV)   / stats(iSV)
        stats(iSurf) = stats(iSurf) / stats(iSV)
        stats(iPrimPar) = stats(iPrimPar) / stats(iSV)
    End Function

    ! -------------------------------------------------------

    Function GetStatsTyped(soln, T, range, N) result (stats)
        ! DESCRIPTION:
        !   Returns an array of statistical properties of
        !   the current ensemble for each particle type.
        ! RETURNS:
        !   Current ensemble statistics.

        Use SWPPARAMS
        Use SWPSTATS_INDICES
        Use SWPSOLN
        Use SWPPART_STATS, ParticleStatCount=>Count, iPM0=>iM0, iPM1=>iM1, &
                           iPM2=>iM2, iPM3=>iM3, iPM4=>iM4, iPM5=>iM5, &
                           iPM6=>iM6, iPM7=>iM7, iPM8=>iM8, iPM9=>iM9, iPP1=>iP1, iPP2=>iP2, &
                           iPV=>iVol, iPV2=>iVol2, iPSurf=>iSurf, iPSurf2=>iSurf2, &
                           iPDiam=>iDiam, iPDiam2=>iDiam2, iPDiam3=>iDiam3, iPDiam4=>iDiam4, & 
                           iPPrimPar=>iPrimPar, iPPrimPar2=>iPrimPar2
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)             :: N              ! Type count to get.
        Double Precision                :: stats(Count,N) ! Return value.
        Type(Solution), Intent(IN)      :: soln           ! Solution for which to calculate stats.
        Real, Intent(IN)                :: T              ! Temperature at which to calculate stats (K).
        Type(StatRange), Intent(IN)     :: range          ! Low and High cut-off values of the given particle
                                                          ! property.  particles outside this range will
                                                          ! not be included in the summation.

        ! VARIABLES.
        Real :: v
        Double Precision :: s(ParticleStatCount,N)

        ! EXECUTABLE CODE.

        ! Pre-calc sample volume so we can scale the
        ! stats to a unit volume.
        v = SampleVolume(soln, T)

        ! Get statistics from the ensemble.
        s = SumStatsTyped(soln%Ensemble, range, N)

        ! Convert to statistics we want to export.
        stats(iSP,:)         = s(iPM0,:)
        stats(iM0,:)         = s(iPM0,:)
        stats(iM1,:)         = s(iPM1,:)
        stats(iM2,:)         = s(iPM2,:)
        stats(iM3,:)         = s(iPM3,:)
        stats(iM4,:)         = s(iPM4,:)
        stats(iM5,:)         = s(iPM5,:)
        stats(iM6,:)         = s(iPM6,:)
        stats(iM7,:)         = s(iPM7,:)
        stats(iM8,:)         = s(iPM8,:)
        stats(iM9,:)         = s(iPM9,:)
        stats(iP1,:)         = s(iPP1,:)
        stats(iP2,:)         = s(iPP2,:)
        stats(iNs,:)         = s(iPM0,:)
        stats(iFV,:)         = s(iPV,:)
        stats(iSurf,:)       = s(iPSurf,:)
        stats(iPrimPar,:)    = s(iPPrimPar,:)
        stats(iSV,:)         = Dble(v)
        stats(iScale,:)      = Dble(soln%LastKnownScale)

        ! Averaged properties.
        Where (stats(iSP,:) > 0)
            stats(iAvgFV,:)      = stats(iFV,:)   / stats(iSP,:)
            stats(iFvSDv,:)      = sqrt(max(0.0D0, s(iPV2,:)/stats(iM0,:) - stats(iAvgFV,:)**2.D0))
            stats(iAvgSurf,:)    = stats(iSurf,:) / stats(iSP,:)
            stats(iSurfSDv,:)    = sqrt(max(0.0D0, (s(iPsurf2,:) / stats(iSP,:)) - stats(iAvgSurf,:)**2.0D0 ))
            stats(iDiam,:)       = s(iPDiam,:) / stats(iSP,:)
            stats(iDiamSDv,:)    = sqrt(max(0.0D0, (s(iPDiam2,:) / stats(iSP,:)) - (s(iPDiam,:) / stats(iSP,:))**2.0D0 ))
            stats(iDiamSk,:)     = max(0.D0, (s(iPDiam3,:) / stats(iSP,:)) - &
                                     3.0D0*stats(iDiam,:)*stats(iDiamSDv,:)**2.0D0 - & 
                                           stats(iDiam,:)**3.0D0)
            stats(iDiamKt,:)     = max(0.0D0, (s(iPDiam4,:) / stats(iSP,:)) - & 
                                     4.0D0*stats(iDiam,:)*stats(iDiamSk,:) - &
                                     6.0D0*stats(iDiam,:)**2.0D0*stats(iDiamSDv,:)**2.0D0 - &
                                           stats(iDiam,:)**4.0D0)
            stats(iAvgPrimPar,:) = stats(iPPrimPar,:) / stats(iSP,:)
            stats(iPrimParSDv,:) = sqrt(max(0.0D0, s(iPPrimPar2,:)/stats(iM0,:) - stats(iAvgPrimPar,:)**2.D0))
        ElseWhere
            stats(iAvgFV,:)      = 0.0D0
            stats(iFvSDv,:)      = 0.0D0
            stats(iAvgSurf,:)    = 0.0D0
            stats(iSurfSDv,:)    = 0.0D0
            stats(iDiam,:)       = 0.0D0
            stats(iDiamSDv,:)    = 0.0D0
            stats(iDiamSk,:)     = 0.0D0
            stats(iDiamKt,:)     = 0.0D0
            stats(iAvgPrimPar,:) = 1.0D0
            stats(iPrimParSDv,:) = 0.0D0
        End Where

        ! Scale stats to unit volume.
        stats(iM0,:)   = stats(iM0,:)   / stats(iSV,:)
        stats(iM1,:)   = stats(iM1,:)   / stats(iSV,:)
        stats(iM2,:)   = stats(iM2,:)   / stats(iSV,:)
        stats(iM3,:)   = stats(iM3,:)   / stats(iSV,:)
        stats(iM4,:)   = stats(iM4,:)   / stats(iSV,:)
        stats(iM5,:)   = stats(iM5,:)   / stats(iSV,:)
        stats(iM6,:)   = stats(iM6,:)   / stats(iSV,:)
        stats(iM7,:)   = stats(iM7,:)   / stats(iSV,:)
        stats(iM8,:)   = stats(iM8,:)   / stats(iSV,:)
        stats(iM9,:)   = stats(iM9,:)   / stats(iSV,:)
        stats(iP1,:)   = stats(iP1,:)   / stats(iSV,:)
        stats(iP2,:)   = stats(iP2,:)   / stats(iSV,:)
        stats(iNs,:)   = stats(iNs,:)   / stats(iSV,:)
        stats(iFV,:)   = stats(iFV,:)   / stats(iSV,:)
        stats(iSurf,:) = stats(iSurf,:) / stats(iSV,:)
        stats(iPrimPar,:) = stats(iPrimPar,:) / stats(iSV,:)
    End Function

    ! -------------------------------------------------------

    Function SumStats(ens, range, typeID) result (stats)
        ! DESCRIPTION:
        !   Gets the accumulative statistics of the particle
        !   ensemble.
        ! RETURNS:
        !   Ensemble statistics.

        Use SWPENSEMBLE
        Use SWPPART_STATS, StatCount=>Count
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens        ! Ensemble to count.
        Double Precision :: stats(StatCount)     ! Return value.
        Type(StatRange), Intent(IN)   :: range   ! Low and High cut-off values of the given particle
                                                 ! property.  particles outside this range will
                                                 ! not be included in the summation.
        Integer, Optional, Intent(IN) :: typeID  ! Optional type of particle for which to 
                                                 ! get statistics.  If omitted gets statistics
                                                 ! for all particles.

        ! VARIABLES.
        Integer :: i, N
        Double Precision :: sz

        ! EXECUTABLE CODE.
        stats = 0.0E0          ! Reset stats array!
        N = ParticleCount(ens) ! Get total particle count.

        If (N > 0) Then
            If (Present(typeID) .And. (typeID>=0)) Then
                ! Sum only particles of given type within cut-off range.
                Do i = 1, N
                    sz = ens%Particles(i)%Properties(range%PropID)
                    If ((ens%Particles(i)%TypeID == typeID) .And. &
                        (sz >= range%Lower) .And. (sz <= range%Upper)) Then
                        stats = stats + GetStatItems(ens%Particles(i))
                    End If
                End Do
            Else
                ! Sum all particles within cut-off range.
                Do i = 1, N
                    sz = ens%Particles(i)%Properties(range%PropID)
                    If ((sz >= range%Lower) .And. (sz <= range%Upper)) Then
                        stats = stats + GetStatItems(ens%Particles(i))
                    End If
                End Do
            End If
        End If
    End Function

    ! -------------------------------------------------------

    Function SumStatsTyped(ens, range, ntypes) result (stats)
        ! DESCRIPTION:
        !   Gets the accumulative statistics of the particle
        !   ensemble for all different particle types.
        ! RETURNS:
        !   Ensemble statistics.

        Use SWPENSEMBLE
        Use SWPPART_STATS, StatCount=>Count
        Implicit None

       ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens ! Ensemble to count.
        Integer, Intent(IN) :: ntypes
        Double Precision :: stats(StatCount, ntypes) ! Return value.
        Type(StatRange), Intent(IN)   :: range       ! Low and High cut-off values of the given particle
                                                     ! property.  particles outside this range will
                                                     ! not be included in the summation.

        ! VARIABLES.
        Integer :: i, N
        Double Precision :: sz, items(StatCount)

        ! EXECUTABLE CODE.
        stats = 0.0E0          ! Reset stats array!
        N = ParticleCount(ens) ! Get total particle count.

        If (N > 0) Then
            Do i = 1, N
                sz = ens%Particles(i)%Properties(range%PropID)
                If ((ens%Particles(i)%TypeID <= ntypes) .And. &
                    (sz >= range%Lower) .And. (sz <= range%Upper)) Then
                    items = GetStatItems(ens%Particles(i))
                    stats(:,ens%Particles(i)%TypeID+1) = &
                    stats(:,ens%Particles(i)%TypeID+1) + items
                End If
            End Do
        End If
    End Function

    ! -------------------------------------------------------

    Function CompilePSL(ens) result (psl)
        ! DESCRIPTION:
        !   Compiles a Particle Size List into a 2D array
        !   from the given ensemble.
        ! RETURNS:
        !   Ensemble PSL.

        Use SWPENSEMBLE
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens ! Ensemble to count.
        Double Precision :: psl(PSL_VAR_COUNT, ens%Capacity) ! Return PSL.

        ! VARIABLES.
        Integer :: i, N

        ! EXECUTABLE CODE.
        psl = 0.0D0
        N = ens%FirstSpace - 1

        If (N > 0) Then
           Do i = 1, N
                psl(1,i) = Dble(ens%Particles(i)%Properties(iV))
                psl(2,i) = Dble(ens%Particles(i)%Surface)
                psl(3,i) = Dble(ens%Particles(i)%Properties(iD))
                psl(4,i) = Dble(ens%Particles(i)%TypeID)
                psl(5,i) = Dble(ens%Particles(i)%CreateTime)
           End Do
        End If
    End Function

    ! -------------------------------------------------------

    Function CompilePCL(ens, mech) result (pcl)
        ! DESCRIPTION:
        !   Compiles a Particle Component List into a 2D array
        !   from the given ensemble.  Component list also contains
        !   tracker variables.
        ! RETURNS:
        !   Ensemble PCL.

        Use SWPENSEMBLE
        Use SWPMECH
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN)  :: ens  ! Ensemble to count.
        Type(Mechanism), Intent(IN) :: mech ! Mechanism with component info.
        Double Precision :: pcl(mech%ComponentCount+mech%TrackerCount, ens%Capacity) ! Return PCL.

        ! VARIABLES.
        Integer :: i, N

        ! EXECUTABLE CODE.
        pcl = 0.0D0
        N = ens%FirstSpace - 1

        If (N > 0) Then
            Do i = 1, N
                pcl(1:mech%ComponentCount,i) = Dble(ens%Particles(i)%Components(1:mech%ComponentCount))
                pcl(mech%ComponentCount+1:mech%ComponentCount+mech%TrackerCount,i) = Dble(ens%Particles(i)%Track(1:mech%TrackerCount))
            End Do
        End If
    End Function
End Module

! *****************************************************************************
!
! File:				swprng.f90
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
!   The random number generator for Sweep.
!
! Functions:
!   SeedRnd     -   Seeds the random number generator.
!   Rnd         -   Generates a uniform real random number in [0,1].
!   IRnd        -   Generates a random integer in [0,Top].
! *****************************************************************************


Module SWPRNG
    ! -------------------------------------------------------
    ! IMPORT PUBLIC MODULE FUNCTIONS.
    ! -------------------------------------------------------
    Use MersenneTwister, only: ignpoi
    ! -------------------------------------------------------

    Implicit none
    Public

    Contains

    ! -------------------------------------------------------

    Subroutine SeedRnd(seed)
        ! DESCRIPTION:
        !    Seeds the random number generator.
        Use MersenneTwister, only: sgrnd
        Implicit None
        Integer, Intent(IN) :: seed
!        Integer :: sd(1)

         Call sgrnd(seed)

!        sd(1) = seed
!        Call Random_Seed(PUT=sd)
    End Subroutine

    ! -------------------------------------------------------

    Real Function Rnd()
        ! DESCRIPTION:
        !  Returns a uniform real random number in [0,1].
        Use MersenneTwister, only: grnd
        Implicit None
        Rnd = Real(grnd())
!        Call Random_Number(Rnd)
    End Function

    ! -------------------------------------------------------

	Real Function IRnd(top)
        ! DESCRIPTION:
        !    Returns a uniform random integer in [1,top].
        Use MersenneTwister, only: igrnd
        Implicit None

        Integer, Intent(IN) :: top
!        Real :: num

        IRnd = IAnd(igrnd(), top-1) + 1

!        Call Random_Number(num)
!        IRnd = Int(num * Real(top) - 0.5E0) + 1
    End Function

End Module

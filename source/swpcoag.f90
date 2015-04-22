! *****************************************************************************
!
! File:					swpcoag.f90
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
!	Defintions of coagulation rates and processes.
!
! Functions:
!	CoagRateTerms		-	Calculates the rate of all coagulation processes
!							and returns them as an array.
! *****************************************************************************

Module SWPCOAG
	! -------------------------------------------------------
	! IMPORT REQUIRED MODULES.
	! -------------------------------------------------------
	Use SWPCOAG_MODEL ! Coagulation model.
	! -------------------------------------------------------

	Implicit None
	Public

	Contains

	! -------------------------------------------------------

	Function CoagRateTerms(T, P, N, sums, vol, flag)
		! DESCRIPTION:
		!	Calculates the rate of coagulation of soot
		!	particles.
		! RETURNS:
		!	The rates of all coagulation processes.

        Use SWPCHEM_SHARED, only: Visc
		Use SWPPART
        Use SWPPARAMS
		Implicit None

		! ARGUMENTS.
		Real ::	CoagRateTerms(NCOAGTERMS)   ! Return values.
		Real, Intent(IN)	 ::	T, P        ! Current chemical conditions.
        Integer, Intent(IN)  :: N           ! Number of particles in ensemble.
		Real, Intent(IN)	 ::	sums(PROPERTY_COUNT) ! Sums of properties of the particle ensemble.
		Real, Intent(IN)	 ::	vol  	    ! Sample volume of the ensemble.
		Integer, Intent(OUT) ::	flag	    ! Error flag.

		! VARIABLES.
		Real :: N_1, a, b, c, sf, fm

		! EXECUTABLE CODE.
		flag = 0

        If (N > 2.0E0) Then
		    ! Calculate some prerequisite values.
		    N_1	= N - 1.0E0
		    a	= C_sf * (T / Visc(T))
		    b	= a * KNUDSEN_K * 1.257E0 * T / P
		    c	= C_fm * Sqrt(T)

		    ! Slip-flow coagulation terms.
		    CoagRateTerms(1) = N * N_1 * a
		    CoagRateTerms(2) = (sums(iD) * sums(iD_1) - N) * a
		    CoagRateTerms(3) = sums(iD_1) * N_1 * b
		    CoagRateTerms(4) = (sums(iD) * sums(iD_2) - sums(iD_1)) * b

		    ! Free-molecular coagulation terms.
		    CoagRateTerms(5) = N_1 * sums(iD2M_1_2) * c
		    CoagRateTerms(6) = (sums(iM_1_2) * sums(iD2) - sums(iD2M_1_2)) * c
!    		CoagRateTerms(7) = 2.0E0 * (sums(i_1_6) * sums(i1_3) - sums(i1_6)) * c
		    CoagRateTerms(7) = 0.0E0

		    ! Majorant as used in DSA versions.
		    CoagRateTerms(5:7) = CoagRateTerms(5:7) * C_fmmaj

		    ! Scale coagulation rates to a unit voume.
		    CoagRateTerms = CoagRateTerms / vol

		    ! Avoid small negative numbers when two 'equal' 
		    ! numbers are subtracted, i.e. in terms 2 and 4.
		    Where(CoagRateTerms < 0.0E0)
			    CoagRateTerms = 0.0E0
		    End Where

		    ! Select which coagulation majorant to use.
		    sf = Sum(CoagRateTerms(1:4))
		    fm = Sum(CoagRateTerms(5:7))

		    If (sf > 0.0E0 .Or. fm > 0.0E0) Then
			    ! There is some coagulation.
			    If (sf > fm) Then
				    ! Use free molecular majorant, ignore
				    ! slip flow terms.
				    CoagRateTerms(1:4) = 0.0E0
			    Else
				    ! Use slip flow majorant, ignore
				    ! free molecular terms.
				    CoagRateTerms(5:7) = 0.0E0
			    End If
		    End If

!            ! These terms for HD soot mech.
!            CoagRateTerms(1:6) = 0.0E0
!            CoagRateTerms(7)   = 4.5E12 * Sqrt(T) * N * N / (NA * vol)
        Else
            CoagRateTerms = 0.0E0
        End If
	End Function

End Module

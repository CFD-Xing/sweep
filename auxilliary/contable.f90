! *****************************************************************************
!
! File:				contable.f90
! Project:			N/A.
! Author(s):			Matthew Celnik (msc37)
!
! Copyright (C) 2006  Matthew S Celnik
!
! Licence:
!   This library is free software; you can redistribute it and/or
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
!   This module contains routines which allow tables of numerical values
!   to be printed to the console.  Various standard formats are included.
!
! Functions:
!   PrintHeader			-	Prints a header row to the console with
!						the given column names.
!   PrintScientific		-	Prints a data row where the columns are numbers
!						given in scientific form.
! *****************************************************************************


Module ConsoleTable
    Implicit None
    Public

    ! -------------------------------------------------------
    ! CONSOLE TABLE PARAMETERS.
    ! -------------------------------------------------------

    ! Maximum number of columns allowed in the table (Note
    ! the leftmost column is also a row number, so this value
    ! is the number of columns not including the leftmost).
    Integer, Parameter :: MAX_COLS = 6

    ! -------------------------------------------------------
    ! STANDARD TABLE FORMATS.
    ! -------------------------------------------------------

    ! Header row format.
    Character(LEN=*), Parameter :: HEAD = "('| ',A4,6(' |  ',A8),' |')"

    ! Dividing lines.
    Character(LEN=*), Parameter :: DIV1 = "('|-',4('-'),6('-|-',9('-')),'-|')"

    ! Numerical data rows.
    Character(LEN=*), Parameter :: SCI          = "('| ',I4.4,6(' | ',ES9.2E2),' |')"
    Character(LEN=*), Parameter :: FLOAT        = "('| ',I4.4,6(' | ',F9.5),' |')"
    Character(LEN=*), Parameter :: INT          = "('| ',I4.4,6(' | ',I9),' |')"

    ! Column formats.
    Character(LEN=*), Parameter :: SCICOL       = "ES9.2E2"
    Character(LEN=*), Parameter :: FLOATCOL     = "F9.5"
    Character(LEN=*), Parameter :: FLOAT2COL    = "F9.2"
    Character(LEN=*), Parameter :: INTCOL       = "I9"

    Contains

    ! -------------------------------------------------------

    Subroutine PrintHeader(leftmost, cols, ncols)
        ! DESCRIPTION:
        !	Prints a header row to the console.

        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)    :: ncols ! Number of columns after the first.
        Character(LEN=*)       :: leftmost, cols(ncols) ! Column headings.

        ! EXECUTABLE CODE.
        Print DIV1
        Print HEAD, leftmost(1:4), cols(1:Min(ncols, MAX_COLS))
        Print DIV1
    End Subroutine

    ! -------------------------------------------------------

    Subroutine PrintScientific(row, vals, nvals)
        ! DESCRIPTION:
        !   Prints a data row using scientific notation for
        !   the numbers.  Requires the values to be given
        !   as reals (single precision).
     
        Implicit None

        ! ARGUMENTS.
        ! Row number (used to determine if header row should
        ! be reprinted).
        Integer, Intent(IN)     :: row
        ! Values in the remaining columns.
        Integer, Intent(IN)     :: nvals
        Real, Intent(IN)        :: vals(nvals)

        ! EXECUTABLE CODE.
        Print SCI, row, vals(1:Min(nvals, MAX_COLS))
    End Subroutine

    ! -------------------------------------------------------

    Subroutine PrintScientificD(row, vals, nvals)
        ! DESCRIPTION:
        !   Prints a data row using scientific notation for
        !   the numbers.  Requires the values to be given
        !   as reals (single precision).
        
        Implicit None

        ! ARGUMENTS.
        ! Row number (used to determine if header row should
        ! be reprinted).
        Integer, Intent(IN) :: row
        ! Values in the remaining columns.
        Integer, Intent(IN) ::nvals
        Double Precision, Intent(IN) :: vals(nvals)

        ! EXECUTABLE CODE.
        Print SCI, row, vals(1:Min(nvals, MAX_COLS))
    End Subroutine

    ! -------------------------------------------------------

    Subroutine PrintCustom(row, vals, nvals, fmt)
        ! DESCRIPTION:
        !   Prints a data row using custom formats for each
        !   column.  Requires the values to be given
        !   as reals (single precision).

        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)          :: row         ! Row number.
        Integer, Intent(IN)          :: nvals       ! Number of columns.
        Real, Intent(IN)             :: vals(nvals) ! Values in the remaining columns.
        Character(LEN=*), Intent(IN) :: fmt         ! Format for output.

        ! EXECUTABLE CODE.
        Print fmt, row, vals(1:Min(nvals, MAX_COLS))
    End Subroutine

    ! -------------------------------------------------------

    Subroutine CreateCustomFormat(fmt, colfmts, nvals)
        ! DESCRIPTION:
        !   Prints a data row using custom formats for each
        !   column.  Requires the values to be given
        !   as reals (single precision).

        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(OUT) :: fmt            ! Generated format
        Integer, Intent(IN)           :: nvals          ! Number of columns.
        Character(LEN=*), Intent(IN)  :: colfmts(nvals) ! Formats for each column.

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        fmt = "('| ',I4.4,"
        Do i = 1, nvals
            fmt = Trim(fmt) // "' | '," // Trim(colfmts(i))
        End Do
        fmt = Trim(fmt) // ",' |')"
    End Subroutine

End Module

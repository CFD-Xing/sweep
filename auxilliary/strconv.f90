! *****************************************************************************
!
! File:				strconv.f90
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
!   This module contains routines for converting numerical values into
!   strings.  This functionality is curiously complicated in F90.
!
!   This module also contains routines for processing file names, such as
!   returning the file extension, or the file name when given a path.
!
! Functions:
!   IntStr			-	Converts an integer to a string of 9 digits.  The string
!						starts at the beginning of the character array.
!   FStr			-	Converts a single precision number to a string of 9
!						characters.
!   ESStr			-	Converts a double precision number to a string of 9
!						characters.
!	---------------------------------------------------------------------------
!   GetFileNameNoExtension	-	Gets the part of the filename before the
!						extension.
!   GetFileExtension	        -	Gets the file extension for the given file name.
!	---------------------------------------------------------------------------
!   IndexOf			-	Returns the index of a string in an array of strings.
!   SplitString                 -       Splits a string into substrings, delimited with the given
!                                               string.
! *****************************************************************************


Module StrConv
    ! -------------------------------------------------------
    ! IMPORT REQUIRED MODULES.
    ! -------------------------------------------------------
    ! None.
    ! -------------------------------------------------------

    Implicit None
    Public

    Interface CStr
        Module Procedure IntStr
        Module Procedure ESStr
    End Interface

    Contains

    ! -------------------------------------------------------
    ! NUMERICAL CONVERSION ROUTINES.
    !
    !   These routines convert numerical values to 9
    !   character strings.  These functions are all wrapped
    !   in the CStr interface above.
    !
    ! -------------------------------------------------------

    Function IntStr(i) Result(sint)

        ! DESCRIPTION:
        !   Converts an integer into a string of 9 characters.
        !
        ! RETURNS:
        !   A 9 character long string representation of the
        !   integer.  String may include blank spaces after
        !   the number.

        ! ARGUMENTS.
        Integer, Intent(IN)     :: i

        ! VARIABLES.
        Character(LEN=6)        :: ft
        Character(LEN=9)        :: sint

        ! EXECUTABLE CODE.
        If (i > 99999999) Then
            ft = "(I9.9)"
        ElseIf (i > 9999999) Then
           ft = "(I8.8)"
        ElseIf (i > 999999) Then
           ft = "(I7.7)"
        ElseIf (i > 99999) Then
           ft = "(I6.6)"
        ElseIf (i > 9999) Then
           ft = "(I5.5)"
        ElseIf (i > 999) Then
           ft = "(I4.4)"
        ElseIf (i > 99) Then
           ft = "(I3.3)"
        ElseIf (i > 9) Then
           ft = "(I2.2)"
        Else
           ft = "(I1.1)"
        End If
        
        sint = ""
        Write (sint, FMT=ft) i
    End Function

    ! -------------------------------------------------------

    Character(LEN=9) Function FStr(float)
        ! DESCRIPTION:
        !   Converts an double into a string of 9 characters.
        !
        ! RETURNS:
        !   A 9 character long string representation of the
        !   integer.  String may include blank spaces after
        !   the number.

        ! ARGUMENTS.

        Real, Intent(IN) :: float

        ! VARIABLES.

        Character(LEN=*), Parameter :: fmt = "ES8.2E2"

        ! EXECUTABLE CODE.

        FStr = ""
        Write (FStr, "(" // fmt // ")") float
    End Function

    ! -------------------------------------------------------

    Character(LEN=9) Function ESStr(float, how)
        ! DESCRIPTION:
        !   Converts an double into a string of 9 characters
        !   in scientific format.
        !
        ! RETURNS:
        !   A 9 character long string representation of the
        !   integer.  String may include blank spaces after
        !   the number.

        ! ARGUMENTS.
        Double Precision, Intent(IN)  :: float
        Integer, Intent(IN), Optional :: how

        ! VARIABLES.

        Character(LEN=*), Parameter :: fmt1 = "ES8.2E2"
        Character(LEN=*), Parameter :: fmt2 = "F8.2"

        ! EXECUTABLE CODE.
        ESStr = ""

        If (Present(how)) Then
            If (how == 1) Then
                Write (ESStr, "(" // fmt2 // ")") float
            Else
                Write (ESStr, "(" // fmt1 // ")") float
            End If
        Else
            Write (ESStr, "(" // fmt1 // ")") float
        End If
    End Function

    ! -------------------------------------------------------

    Real Function CReal(str, def)
        ! DESCRIPTION:
        !   Converts an a string into a real number, returning
        !   a default value if there is an error.
        !
        ! RETURNS:
        !   A 9 character long string representation of the
        !   integer.  String may include blank spaces after
        !   the number.

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)  :: str ! String to convert.
		Real, Intent(IN), Optional :: def ! Default value.

        ! VARIABLES.
        Integer :: err, x
        Real    :: defval

        ! EXECUTABLE CODE.

        If (Present(def)) Then
            defval = def
        Else
            defval = 0.0E0
        End If

        If (IsInteger(str)) Then
            ! Read as an integer.
            err = 0
            Read(str, FMT=*, IOSTAT=err) x
            
            If (err /= 0) Then
                CReal = defval  ! Not a number.
            Else
                CReal = Real(x) ! It's an integer!
            End If
        Else
            ! Try to read value as a float.
            err = 0
            Read(str(1:Len_Trim(str)), FMT=*, IOSTAT=err) CReal
            If (err /= 0) CReal = defval
        End If
    End Function

    ! -------------------------------------------------------

    Logical Function IsNumeric(str)
        ! DESCRIPTION:
        !   Checks a string to see if it represents a number.
        ! RETURNS:
        !   TRUE  = String is a number.
        !   FALSE = String is not a number.

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: str ! String to test.

        ! VARIABLES.
        Integer :: i, j1, jdot, jexp, L

        ! EXECUTABLE CODE.
        
        If (Len_Trim(str) < 1) Then
            IsNumeric = .False.
            Return
        End If

        ! First check to see if string is an integer.
        IsNumeric = IsInteger(str)

        If (.Not. IsNumeric) Then
            ! String is not an integer, but might be a float.  First
            ! check valid starting characters.
            If (Scan(str(1:1), ".-+1234567890") == 0) Then
                ! Invalid starting character.
                IsNumeric = .False.
            Else
                L = Len_Trim(str) ! Length of string.

                ! Get location of decimal point.
                jdot = Index(str, ".")
                
                ! Get location of exponent marker.
                jexp = Scan(str,"dDeE")

                If (jdot <= L) Then
                    ! There is a decimal point in the string.
                    If ((jexp > 0) .And. (jexp < L)) Then
                        ! There is an exponent marker in the string.
                        If (jexp < jdot) Then
                            ! Decimal point is in the exponent, this is
                            ! wrong.
                            IsNumeric = .False.
                        Else
                            IsNumeric = IsInteger(str(1:jdot-1))
                            If (IsNumeric .And. (jexp > jdot + 1)) Then
                                IsNumeric = IsInteger(str(jdot+1:jexp-1))
                            End If
                            If (IsNumeric .And. (Scan(str(jexp+1:jexp+1), "-+1234567890") /= 0)) Then
                                If (jexp+2 <= L) IsNumeric = IsInteger(str(jexp+2:))
                            Else
                                IsNumeric = .False.
                            End If
                        End If
                    ElseIf (jexp == L) Then
                        ! Exponent marker is at end of string.  This
                        ! is wrong.
                        IsNumeric = .False.
                    Else
                        ! There is no exponent marker in the string.
                        IsNumeric = IsInteger(str(1:jdot-1))
                        If (IsNumeric .And. (jdot+1 <= L)) IsNumeric = IsInteger(str(jdot+1:))
                    End If
                Else
                    ! There is no decimal point in the string (or it as it the beginning).
                    If ((jexp > 0) .And. (jexp < L)) Then
                        ! There is an exponent marker in the string.
                        IsNumeric = IsInteger(str(1:jexp-1))
                        If (IsNumeric .And. (Scan(str(jexp+1:jexp+1), "-+1234567890") /= 0)) Then
                            If (jexp+2 <= L) IsNumeric = IsInteger(str(jexp+2:))
                        Else
                            IsNumeric = .False.
                        End If
                    Else
                        ! There is no exponent marker in the string
                        ! or it is at the end of the string.  Both wrong.
                        IsNumeric = .False.
                    End If
                End If
            End If
        End If
    End Function

    ! -------------------------------------------------------

    Logical Function IsInteger(str)
        ! DESCRIPTION:
        !   Checks a string to see if it represents an integer.
        ! RETURNS:
        !   TRUE  = String is an integer.
        !   FALSE = String is not an integer.

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: str ! String to test.

        ! VARIABLES.
        Integer :: i, L
        Character :: c

        ! EXECUTABLE CODE.
        IsInteger = .True.

        ! Check first character.
        If (Scan(str, "-+1234567890") /= 1) Then
            IsInteger = .False.
        Else
            L = Len_Trim(str)
            Do i = 2, L
                If (Scan(str(i:), "1234567890") /= 1) Then
                    IsInteger = .False.
                    Exit
                End If
            End Do
        End If
    End Function

    ! -------------------------------------------------------
    ! FILE NAME PROCESSING ROUTINES
    !
    !   These routines process file names.  For example
    !   returning the file extension given the file name.
    !
    ! -------------------------------------------------------

    Subroutine GetFileNameNoExtension(file, name)

        ! DESCRIPTION:
        !   Takes a file name and returns the part before
        !   the extension (not including the full-stop).
        !
        ! RETURNS:
        !   The file name with the file extension removed.

        ! ARGUMENTS.

        Character*(*), Intent(IN)  :: file
        Character*(*), Intent(OUT) :: name

        ! VARIABLES.

        Integer :: idot

        ! EXECUTABLE CODE.

        idot = Index(file, ".", .True.)

        If ((idot > 0) .And. (idot <= Len_Trim(file))) Then
           name = file(1:idot-1)
        End If

    End Subroutine

    ! -------------------------------------------------------

    Subroutine GetFileExtension(file, ext)

        ! DESCRIPTION:
        !   Takes a file name and returns the extension, that
        !   is the part after the last full-stop.
        !
        ! RETURNS:
        !   The file extension.

        ! ARGUMENTS.

        Character*(*), Intent(IN)  :: file
        Character*(*), Intent(OUT) :: ext

        ! VARIABLES.

        Integer :: idot

        ! EXECUTABLE CODE.

        idot = Index(file, ".", .True.)

        If ((idot > 0) .And. (idot <= Len_Trim(file))) Then
            ext = file(idot+1:Len_Trim(file))
        End If

    End Subroutine

    ! -------------------------------------------------------

    Integer Function IndexOf(s, list)

        ! DESCRIPTION:
        !   Locates a string in a list of strings, returns
        !   0 if the string is not present.
        !
        ! RETURNS:
        !   Index of string in list.

        Implicit None

        ! ARGUMENTS.

        ! String to locate.
        Character*(*), Intent(IN) :: s
        ! List to find it in.
        Character*(*), Intent(IN) :: list(:)

        ! VARIABLES.

        Integer :: L, LTS, LTL, N, i, j
        Logical :: found

        ! EXECUTABLE CODE.

        N = Size(list)

        If (N > 0) Then
            IndexOf = 0

            ! Get the trimmed length of the
            ! string to find.
            LTS = Len_Trim(s)

            ! Get the minimum string length.
            L = Min(LTS, Len(list(1)))

            ! Loop through string in list comparing
            ! each character.
            Do i = 1, N
                LTL = Len_Trim(list(i))

                ! Check that the string length
                ! in the lsit is the same as the search
                ! string length.
                If (LTL == L) Then
                    found = .True.

                    ! Loop through all the characters
                    ! in the strings.
                        Do j = 1, L
                            If (s(j:j) /= list(i)(j:j)) Then
                                found = .False.
                                Exit
                            End If
                        End Do

                        If (found) Then
                           IndexOf = i
                        Return
                    End If
                End If
            End Do
        Else
           ! No items in list.
           IndexOf = 0
        End If
    End Function

    ! -------------------------------------------------------

    Subroutine SplitString(str, delim, parts, lparts, nparts, mergeDelims, blanks)
        ! DESCRIPTION:
        !   Splits a string into substrings using the given
        !   delimiter.

        Implicit None

        ! ARGUMENTS.
        Character*(*), Intent(IN)     :: str, delim     ! The string to split and the delimiter.
        Integer, Intent(IN)           :: lparts, nparts ! Length of a substring and max. number
                                                        ! of substrings.
        Character(LEN=lparts), Intent(OUT) :: parts(nparts)  ! Array of substrings to return.
        Logical, Intent(IN)           :: mergeDelims    ! Treat consecutive delimeters as one.
        Logical, Optional, Intent(IN) :: blanks         ! FALSE=Remove leading/trailing blanks from
                                                        ! substrings.
                                                        ! TRUE=Leave blanks in substrings.

        ! VARIABLES.
        Integer :: N, L, LD, i, pos0, pos, pos1
        Logical :: trm

        ! EXECUTABLE CODE.
        N  = 0
        L  = Len_Trim(str)
        LD = Len(delim)
        parts = ""

        If (Present(blanks)) Then
            trm = .Not. blanks
        Else
            trm = .True.
        End If
        
        ! Start index of search.
        pos0 = 1

        Do i = 1, nparts
            ! Get position of next delimiter.
            pos = Index(str(pos0:), delim)
            
            ! Calculate end of substring.
            pos1 = pos0 + Min(pos, lparts) - 2

            ! Save substring.
            If (trm) Then
                parts(i) = AdjustL(str(pos0:pos1))
            Else
                parts(i) = str(pos0:pos1)
            End If

            ! Adjust starting index.
            pos0 = pos1 + LD + 1

            ! Check for end of delimiters.
            If (pos0 > L) Exit

            ! Check for consecutive delimeters.
            If (mergeDelims) Then
                pos = pos0
                Do While (str(pos:pos) == delim)
                    pos = pos + 1
                End Do
                pos0 = pos
            End If
        End Do
    End Subroutine
End Module

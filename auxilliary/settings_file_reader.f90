! *****************************************************************************
!
! File:				settings_file_reader.f90
! Project:			Coupled Chemistry-Soot Solver.
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
!   A standard text file parser for files of the form:
!   
!   KEYWORD param1 param2 ...
!
!
! Functions:
! *****************************************************************************

Module SettingsFileReader
    Implicit None
    Private
    Public :: LoadSettingsFile, GetParams

    ! -------------------------------------------------------
    ! PARAMETERS.
    ! -------------------------------------------------------

    Integer, Private, Parameter :: LKEY = 20  ! Length of a keyword string.
    Integer, Private, Parameter :: LSTR = 200 ! Length of a descriptive string.
    Integer, Private, Parameter :: MAX_NPARAMS = 200 ! Max. number of parameters after a keyword.
    Integer, Private, Parameter :: MAX_NLINES  = 200 ! Max. number of lines in a file.

    ! -------------------------------------------------------
    ! DERIVED TYPES.
    ! -------------------------------------------------------

    Type, Public :: LineData
        Character(LEN=LKEY) :: Key=""
        Character(LEN=LSTR) :: Params(MAX_NPARAMS)=""
    End Type

    Type, Public :: FileData
        Character(LEN=LSTR) :: Name=""           ! File name.
        Integer             :: NLines=0          ! Line count in file.
        Type(LineData)      :: Lines(MAX_NLINES) ! Lines in file.
    End Type

    ! -------------------------------------------------------
    ! ROUTINE OVERLOADS.
    ! -------------------------------------------------------

    Interface GetParams
        Module Procedure GetParams_Str
        Module Procedure GetParams_StrArray
        Module Procedure GetParams_Real
        Module Procedure GetParams_RealArray
        Module Procedure GetParams_Int
        Module Procedure GetParams_IntArray
        Module Procedure GetParams_MixedSRI
        Module Procedure GetParams_MixedSR
        Module Procedure GetParams_MixedSI
        Module Procedure GetParams_MixedRI
        Module Procedure GetParams_MixedRIArray
    End Interface

    Contains

    ! -------------------------------------------------------

    Subroutine LoadSettingsFile(file, unit, fdata, flag)
        ! DESCRIPTION:
        !   Parses a settings file and stores the data
        !   in a FileData type, supplied by the calling code.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character*(*), Intent(IN)   :: file  ! File name to read.
        Integer, Intent(IN)         :: unit  ! Unit number on which to open file.
        Type(FileData), Intent(OUT) :: fdata ! Settings data from file.
        Integer, Intent(OUT)        :: flag  ! Error flag.

        ! VARIABLES.
        Integer :: i, err
        Character(LEN=5000) :: line
        Character(LEN=LSTR) :: subs(MAX_NPARAMS+1)

        ! EXECUTABLE CODE.
        flag = 0
        i = 0

        ! Open file for reading.
        Open(UNIT=unit, FILE=file, FORM="FORMATTED", &
             STATUS="OLD", ACTION="READ", IOSTAT=err)

        If (err /= 0) Then
            ! Failed to open file.
            flag = -1
            Return
        End If


        Do
            ! Read line from text file.
            line = ""
            Read(UNIT=unit, FMT="(A)", IOSTAT=err) line

            If (err == -2) Then
                ! EOR: not applicable here.
            ElseIf (err == -1) Then
                ! EOF: we have read all the data rows so
                ! exit loop.
                Exit
            ElseIf (err == 0) Then
                ! Normal read, no problems.
            ElseIf (err > 0) Then
                ! Error occured when reading.
                flag = -2
                Exit
            End If


            ! Split the line into substrings.
            Call SplitString(line, " ", subs, LSTR, MAX_NPARAMS+1, .True.)

            ! Save line in FileData.
            i = i + 1
            fdata%Lines(i)%Key    = subs(1)
            fdata%Lines(i)%Params = subs(2:)
        End Do

        ! Save file name and line count.
        fdata%Name   = file
        fdata%NLines = i

        ! Close input file.
        Close (UNIT=unit, IOSTAT=err)
    End Subroutine

    ! -------------------------------------------------------

    Integer Function GetParams_Str(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameter as a string.
        !   If the keyword was not found then returns a
        !   blank strings and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)  :: key       ! Keyword to find.
        Character*(*), Intent(OUT) :: params    ! Parameters returned.
        Integer, Intent(OUT)       :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, i1
        Character(LEN=LSTR) :: subs(MAX_NPARAMS+1)

        ! EXECUTABLE CODE.
        GetParams_Str = 0
        flag = -1
        params = ""

        If (Present(startLine)) Then
            i1 = Max(startLine, 1)
        Else
            i1 = 1
        End If

        Do i = i1, fdata%NLines
            If (fdata%Lines(i)%Key == key) Then
                flag = 0
                GetParams_Str = i
                params = fdata%Lines(i)%Params(1)
                Exit
            End If
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_StrArray(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameters as a list of strings.
        !   If the keyword was not found then returns a list
        !   of blank strings and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)  :: key       ! Keyword to find.
        Character*(*), Intent(OUT) :: params(:) ! Parameters returned.
        Integer, Intent(OUT)       :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, i1

        ! EXECUTABLE CODE.
        GetParams_StrArray = 0
        flag = -1
        params = ""

        If (Present(startLine)) Then
            i1 = Max(startLine, 1)
        Else
            i1 = 1
        End If

        Do i = i1, fdata%NLines
            If (fdata%Lines(i)%Key == key) Then
                flag = 0
                GetParams_StrArray = i
                params = fdata%Lines(i)%Params
                Exit
            End If
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_Real(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameters as a list of double reals.
        !   If the keyword was not found then returns a list
        !   of 0.0 and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)  :: key       ! Keyword to find.
        Double Precision, Intent(OUT) :: params ! Parameters returned.
        Integer, Intent(OUT)       :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i
        Character(LEN=LSTR) :: sparams

        ! EXECUTABLE CODE.
        flag = -1
        params = 0.0E0
        sparams = ""

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_Real = GetParams_Str(fdata, key, sparams, flag, startLine)
        Else
            GetParams_Real = GetParams_Str(fdata, key, sparams, flag)
        End If
        If (flag < 0) Return

        ! Convert parameter to real.
        If (sparams == "") Return

        If (IsNumeric(sparams)) Then
            params = CReal(sparams)
        Else
            flag = -2
            params = 0.0E0
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_RealArray(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameters as a list of double reals.
        !   If the keyword was not found then returns a list
        !   of 0.0 and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN)    :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)     :: key       ! Keyword to find.
        Double Precision, Intent(OUT) :: params(:) ! Parameters returned.
        Integer, Intent(OUT)          :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i
        Character(LEN=LSTR) :: sparams(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        params = 0.0E0
        sparams = ""

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_RealArray = GetParams_StrArray(fdata, key, sparams, flag, startLine)
        Else
            GetParams_RealArray = GetParams_StrArray(fdata, key, sparams, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to reals.
        Do i = 1, MAX_NPARAMS
            If (sparams(i) == "") Exit

            If (IsNumeric(sparams(i))) Then
                params(i) = CReal(sparams(i))
            Else
                flag = -2
                params(i) = 0.0E0
            End If
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_Int(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameters as a list of integers.
        !   If the keyword was not found then returns a list
        !   of 0 and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)  :: key       ! Keyword to find.
        Integer, Intent(OUT)       :: params    ! Parameters returned.
        Integer, Intent(OUT)       :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, err
        Character(LEN=LSTR) :: sparams

        ! EXECUTABLE CODE.
        flag = -1
        params = 0
        sparams = ""

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_Int = GetParams_Str(fdata, key, sparams, flag, startLine)
        Else
            GetParams_Int = GetParams_Str(fdata, key, sparams, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to integers.
        If (sparams == "") Return

        If (IsInteger(sparams)) Then
            Read(sparams, FMT=*, IOSTAT=err) params
        Else
            flag = -2
            params = 0
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_IntArray(fdata, key, params, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   and returns the parameters as a list of integers.
        !   If the keyword was not found then returns a list
        !   of 0 and sets flag=-1.
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)  :: key       ! Keyword to find.
        Integer, Intent(OUT)       :: params(:) ! Parameters returned.
        Integer, Intent(OUT)       :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, err
        Character(LEN=LSTR) :: sparams(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        params = 0
        sparams = ""

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_IntArray = GetParams_StrArray(fdata, key, sparams, flag, startLine)
        Else
            GetParams_IntArray = GetParams_StrArray(fdata, key, sparams, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to integers.
        Do i = 1, MAX_NPARAMS
            If (sparams(i) == "") Exit

            If (IsInteger(sparams(i))) Then
                Read(sparams(i), FMT=*, IOSTAT=err) params(i)
            Else
                flag = -2
                params(i) = 0
            End If
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function GetParams_MixedSRI(fdata, key, fmt, sparams, rparams, iparams, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   andd splits the parameters up into strings, reals
        !   or integers based on the fmt string.  The fmt string
        !   contains one character for each expected parameter
        !   where S=String, R=Real and I=Integer.  Hence the 
        !   following string reads two reals followed by an integer
        !   followed by three strings:
        !
        !   fmt = "RRISSS"
        !
        !   Easy!
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata      ! Settings data from file.
        Character*(*), Intent(IN)  :: key        ! Keyword to find.
        Character*(*), Intent(IN)  :: fmt        ! Formats of parameters to read.
        Character*(*), Intent(OUT) :: sparams(:) ! Parameters returned.
        Double Precision, Intent(OUT) :: rparams(:) ! ''
        Integer, Intent(OUT)       :: iparams(:) ! ''
        Integer, Intent(OUT)       :: flag       ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, is, ir, ii, err, nparams
        Character(LEN=LSTR) :: strs(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        sparams = ""
        rparams = 0.0E0
        iparams = 0
        strs = ""
        is = 0
        ir = 0
        ii = 0
        nparams = Len_Trim(fmt)

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_MixedSRI = GetParams_StrArray(fdata, key, strs, flag, startLine)
        Else
            GetParams_MixedSRI = GetParams_StrArray(fdata, key, strs, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to required formats.
        Do i = 1, Min(nparams, MAX_NPARAMS)
            If (strs(i) == "") Exit

            Select Case (fmt(i:i))
                Case ("S","s")
                    is = is + 1
                    sparams(is) = sparams(i)
                Case ("R","r")
                    ir = ir + 1

                    If (IsNumeric(strs(i))) Then
                        rparams(ir) = CReal(strs(i))
                    Else
                        flag = -2
                        rparams(ir) = 0.0E0
                    End If
                Case ("I","i")
                    ii = ii + 1
                    If (IsInteger(strs(i))) Then
                        Read(strs(i), FMT=*, IOSTAT=err) iparams(ii)
                    Else
                        flag = -2
                        iparams(ii) = 0
                    End If
                Case Default
                    ! A bad fmt character found!
                    flag = -3
                    Exit
            End Select
        End Do
    End Function      

    ! -------------------------------------------------------

    Integer Function GetParams_MixedSR(fdata, key, fmt, sparams, rparams, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   andd splits the parameters up into strings or reals
        !   based on the fmt string.  The fmt string
        !   contains one character for each expected parameter
        !   where S=String, R=Real.  Hence the 
        !   following string reads two reals followed by three strings:
        !
        !   fmt = "RRSSS"
        !
        !   Easy!
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata      ! Settings data from file.
        Character*(*), Intent(IN)  :: key        ! Keyword to find.
        Character*(*), Intent(IN)  :: fmt        ! Formats of parameters to read.
        Character*(*), Intent(OUT) :: sparams(:) ! Parameters returned.
        Double Precision, Intent(OUT) :: rparams(:) ! ''
        Integer, Intent(OUT)       :: flag       ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, is, ir, err, nparams
        Character(LEN=LSTR) :: strs(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        sparams = ""
        rparams = 0.0E0
        strs = ""
        is = 0
        ir = 0
        nparams = Len_Trim(fmt)

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_MixedSR = GetParams_StrArray(fdata, key, strs, flag, startLine)
        Else
            GetParams_MixedSR = GetParams_StrArray(fdata, key, strs, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to required formats.
        Do i = 1, Min(nparams, MAX_NPARAMS)
            If (strs(i) == "") Exit

            Select Case (fmt(i:i))
                Case ("S","s")
                    is = is + 1
                    sparams(is) = sparams(i)
                Case ("R","r")
                    ir = ir + 1

                    If (IsNumeric(strs(i))) Then
                        rparams(ir) = CReal(strs(i))
                    Else
                        flag = -2
                        rparams(ir) = 0.0E0
                    End If
                Case Default
                    ! A bad fmt character found!
                    flag = -3
                    Exit
            End Select
        End Do
    End Function      

    ! -------------------------------------------------------

    Integer Function GetParams_MixedSI(fdata, key, fmt, sparams, iparams, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   andd splits the parameters up into strings
        !   or integers based on the fmt string.  The fmt string
        !   contains one character for each expected parameter
        !   where S=String, I=Integer.  Hence the 
        !   following string reads an integer
        !   followed by three strings:
        !
        !   fmt = "ISSS"
        !
        !   Easy!
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata      ! Settings data from file.
        Character*(*), Intent(IN)  :: key        ! Keyword to find.
        Character*(*), Intent(IN)  :: fmt        ! Formats of parameters to read.
        Character*(*), Intent(OUT) :: sparams(:) ! Parameters returned.
        Integer, Intent(OUT)       :: iparams(:) ! ''
        Integer, Intent(OUT)       :: flag       ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, is, ii, err, nparams
        Character(LEN=LSTR) :: strs(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        sparams = ""
        iparams = 0
        strs = ""
        is = 0
        ii = 0
        nparams = Len_Trim(fmt)

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_MixedSI = GetParams_StrArray(fdata, key, strs, flag, startLine)
        Else
            GetParams_MixedSI = GetParams_StrArray(fdata, key, strs, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to required formats.
        Do i = 1, Min(nparams, MAX_NPARAMS)
            If (strs(i) == "") Exit

            Select Case (fmt(i:i))
                Case ("S","s")
                    is = is + 1
                    sparams(is) = sparams(i)
                Case ("I","i")
                    ii = ii + 1
                    If (IsInteger(strs(i))) Then
                        Read(strs(i), FMT=*, IOSTAT=err) iparams(ii)
                    Else
                        flag = -2
                        iparams(ii) = 0
                    End If
                Case Default
                    ! A bad fmt character found!
                    flag = -3
                    Exit
            End Select
        End Do
    End Function      

    ! -------------------------------------------------------

    Integer Function GetParams_MixedRI(fdata, key, fmt, rparams, iparams, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   andd splits the parameters up into reals
        !   or integers based on the fmt string.  The fmt string
        !   contains one character for each expected parameter
        !   where R=Real and I=Integer.  Hence the 
        !   following string reads two reals followed by an integer:
        !
        !   fmt = "RRI"
        !
        !   Easy!
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN)    :: fdata     ! Settings data from file.
        Character*(*), Intent(IN)     :: key       ! Keyword to find.
        Character*(*), Intent(IN)     :: fmt       ! Formats of parameters to read.
        Double Precision, Intent(OUT) :: rparams   ! Parameters returned.
        Integer, Intent(OUT)          :: iparams   ! ''
        Integer, Intent(OUT)          :: flag      ! Error flag.
        Integer, Optional, Intent(IN) :: startLine ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, ir, ii, err, nparams
        Character(LEN=LSTR) :: strs(2)

        ! EXECUTABLE CODE.
        flag = -1
        rparams = 0.0E0
        iparams = 0
        strs = ""
        ir = 0
        ii = 0

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_MixedRI = GetParams_StrArray(fdata, key, strs, flag, startLine)
        Else
            GetParams_MixedRI = GetParams_StrArray(fdata, key, strs, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to required formats.
        Do i = 1, 2
            If (strs(i) == "") Exit

            Select Case (fmt(i:i))
                Case ("R","r")
                     If (IsNumeric(strs(i))) Then
                        rparams = CReal(strs(i))
                    Else
                        flag = -2
                        rparams = 0.0E0
                    End If
                Case ("I","i")
                    If (IsInteger(strs(i))) Then
                        Read(strs(i), FMT=*, IOSTAT=err) iparams
                    Else
                        flag = -2
                        iparams = 0
                    End If
                Case Default
                    ! A bad fmt character found!
                    flag = -3
                    Exit
            End Select
        End Do
    End Function      

    ! -------------------------------------------------------

    Integer Function GetParams_MixedRIArray(fdata, key, fmt, rparams, iparams, flag, startLine)
        ! DESCRIPTION:
        !   Searches a FileData type for the given keyword,
        !   andd splits the parameters up into reals
        !   or integers based on the fmt string.  The fmt string
        !   contains one character for each expected parameter
        !   where R=Real and I=Integer.  Hence the 
        !   following string reads two reals followed by an integer:
        !
        !   fmt = "RRI"
        !
        !   Easy!
        ! RETURNS:
        !   Line number at which keyword was found.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Type(FileData), Intent(IN) :: fdata         ! Settings data from file.
        Character*(*), Intent(IN)  :: key           ! Keyword to find.
        Character*(*), Intent(IN)  :: fmt           ! Formats of parameters to read.
        Double Precision, Intent(OUT) :: rparams(:) ! Parameters returned.
        Integer, Intent(OUT)       :: iparams(:)    ! ''
        Integer, Intent(OUT)       :: flag          ! Error flag.
        Integer, Optional, Intent(IN) :: startLine  ! Line number to start search at.

        ! VARIABLES.
        Integer :: i, ir, ii, err, nparams
        Character(LEN=LSTR) :: strs(MAX_NPARAMS)

        ! EXECUTABLE CODE.
        flag = -1
        rparams = 0.0E0
        iparams = 0
        strs = ""
        ir = 0
        ii = 0
        nparams = Len_Trim(fmt)

        ! Get the parameter list as strings.
        If (Present(startLine)) Then
            GetParams_MixedRIArray = GetParams_StrArray(fdata, key, strs, flag, startLine)
        Else
            GetParams_MixedRIArray = GetParams_StrArray(fdata, key, strs, flag)
        End If
        If (flag < 0) Return

        ! Convert parameters to required formats.
        Do i = 1, Min(nparams, MAX_NPARAMS)
            If (strs(i) == "") Exit

            Select Case (fmt(i:i))
                Case ("R","r")
                    ir = ir + 1

                    If (IsNumeric(strs(i))) Then
                        rparams(ir) = CReal(strs(i))
                    Else
                        flag = -2
                        rparams(ir) = 0.0E0
                    End If
                Case ("I","i")
                    ii = ii + 1
                    If (IsInteger(strs(i))) Then
                        Read(strs(i), FMT=*, IOSTAT=err) iparams(ii)
                    Else
                        flag = -2
                        iparams(ii) = 0
                    End If
                Case Default
                    ! A bad fmt character found!
                    flag = -3
                    Exit
            End Select
        End Do
    End Function     

End Module

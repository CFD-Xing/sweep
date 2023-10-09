! *****************************************************************************
!
! File:				dict_t.f90
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
!   An implementation of a variable length dictionary of different types.  Dictionary
!   entries have a key (string) and a data value. Currently
!   included data value types are:
!
!   > Double Precision (DICT_DP).
!   > Integer (DICT_INT).
!
! Functions:
!   Destroy     -   Clears memory associated with a dict type.
!   SetP        -   Sets the pointer of ones dict's items to another dict's items.
!   Copy        -   Copies the items from one dict to another.
!   Add         -   Adds two dicts together and returns a new dict (this is safe
!                   for returning one of the dicts being added).
!   AddItem     -   Adds a single item to a dict and returns a new dict (this is safe
!                   for returning the RH dict).
!   AddItems    -   Adds a multiple items to a dict and returns a new dict (this is
!                   safe for returning the RH dict).
! *****************************************************************************

Module DICT_DP_T
    Implicit None
    Public

    Type DICT_DP_ENTRY
        Character(LEN=20) :: Key
        Double Precision  :: Value
    End Type

    Type DICT_DP
        Integer :: Count = 0
        Type(DICT_DP_ENTRY), Pointer :: Items(:)
    End Type

    Interface Assignment(=)
        Module Procedure SetP
    End Interface

    Interface Operator(+)
        Module Procedure Add
        Module Procedure AddItem
        Module Procedure AddEntry
        Module Procedure AddItems
    End Interface

    Contains

    ! -------------------------------------------------------

    Subroutine Destroy(dict)
        ! DESCRIPTION
        !   Clears memory associated with a dictionary type.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: dict

        ! VARIABLES.
        Integer :: err = 0

        ! EXECUTABLE CODE.
        dict%Count = 0
        Deallocate(dict%Items, STAT=err)
        Nullify(dict%Items)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine SetP(new, old)
        ! DESCRIPTION
        !   Assigns the dict of one dict type to another
        !   by passing the dict pointer.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: new
        Type(DICT_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items => old%Items
    End Subroutine

    ! -------------------------------------------------------

    Subroutine Copy(new, old)
        ! DESCRIPTION
        !   Copies the items of one dict type to another.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: new
        Type(DICT_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items = old%Items
    End Subroutine

    ! -------------------------------------------------------

    Function Add(dict1, dict2) Result (sumdict)
        ! DESCRIPTION
        !   Adds together two dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN) :: dict1
        Type(DICT_DP), Intent(IN) :: dict2
        Type(DICT_DP) :: sumdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict1%Count+dict2%Count)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict1%Count > 0) items(1:dict1%Count) = dict1%Items
        If (dict2%Count > 0) items(dict1%Count+1:) = dict2%Items

        ! Resize dict.        
        Deallocate(sumdict%Items, STAT=err)
        sumdict%Count = dict1%Count + dict2%Count
        Allocate(sumdict%Items(sumdict%Count), STAT=err)

        ! Assign back the dict.
        sumdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItem(dict, key, value) Result (newdict)
        ! DESCRIPTION
        !   Adds a number to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)     :: dict
        Character(LEN=20), Intent(IN) :: key
        Double Precision, Intent(IN)  :: value
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) items(1:dict%Count) = dict%Items
        items(dict%Count+1)%Key   = key
        items(dict%Count+1)%Value = value

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + 1
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddEntry(dict, ent) Result (newdict)
        ! DESCRIPTION
        !   Adds a number to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)       :: dict
        Type(DICT_DP_ENTRY), Intent(IN) :: ent
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) items(1:dict%Count) = dict%Items
        items(dict%Count+1) = ent

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + 1
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItems(dict, keys, values) Result (newdict)
        ! DESCRIPTION
        !   Adds some numbers to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)     :: dict
        Character(LEN=20), Intent(IN) :: keys(:)
        Double Precision, Intent(IN)  :: values(:)
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: i, err = 0
        Integer :: N = Min(Size(keys), Size(values))
        Type(DICT_DP_ENTRY) :: allitems(dict%Count+N)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) allitems(1:dict%Count) = dict%Items

        Do i = 1, N
            allitems(dict%Count+i)%Key   = keys(i)
            allitems(dict%Count+i)%Value = values(i)
        End Do

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + N
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = allitems
    End Function

    ! -------------------------------------------------------

    Function AddEntries(dict, entries) Result (newdict)
        ! DESCRIPTION
        !   Adds some numbers to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)       :: dict
        Type(DICT_DP_ENTRY), Intent(IN) :: entries(:)
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: i, err = 0
        Integer :: N = Size(entries)
        Type(DICT_DP_ENTRY) :: allitems(dict%Count+N)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) allitems(1:dict%Count) = dict%Items
        allitems(dict%Count+1:) = entries

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + N
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = allitems
    End Function

End Module

! =======================================================

Module DICT_INT_T
    Implicit None
    Public

    Type DICT_INT_ENTRY
        Character(LEN=20) :: Key
        Integer  :: Value
    End Type

    Type DICT_INT
        Integer :: Count = 0
        Type(DICT_INT_ENTRY), Pointer :: Items(:)
    End Type

    Interface Assignment(=)
        Module Procedure SetP
    End Interface

    Interface Operator(+)
        Module Procedure Add
        Module Procedure AddItem
        Module Procedure AddItems
    End Interface

    Contains

    ! -------------------------------------------------------

    Subroutine Destroy(dict)
        ! DESCRIPTION
        !   Clears memory associated with a dictionary type.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: dict

        ! VARIABLES.
        Integer :: err = 0

        ! EXECUTABLE CODE.
        dict%Count = 0
        Deallocate(dict%Items, STAT=err)
        Nullify(dict%Items)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine SetP(new, old)
        ! DESCRIPTION
        !   Assigns the dict of one dict type to another
        !   by passing the dict pointer.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: new
        Type(DICT_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items => old%Items
    End Subroutine

    ! -------------------------------------------------------

    Subroutine Copy(new, old)
        ! DESCRIPTION
        !   Copies the items of one dict type to another.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(INOUT) :: new
        Type(DICT_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items = old%Items
    End Subroutine

    ! -------------------------------------------------------

    Function Add(dict1, dict2) Result (sumdict)
        ! DESCRIPTION
        !   Adds together two dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN) :: dict1
        Type(DICT_DP), Intent(IN) :: dict2
        Type(DICT_DP) :: sumdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict1%Count+dict2%Count)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict1%Count > 0) items(1:dict1%Count) = dict1%Items
        If (dict2%Count > 0) items(dict1%Count+1:) = dict2%Items

        ! Resize dict.        
        Deallocate(sumdict%Items, STAT=err)
        sumdict%Count = dict1%Count + dict2%Count
        Allocate(sumdict%Items(sumdict%Count), STAT=err)

        ! Assign back the dict.
        sumdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItem(dict, key, value) Result (newdict)
        ! DESCRIPTION
        !   Adds a number to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)     :: dict
        Character(LEN=20), Intent(IN) :: key
        Integer, Intent(IN)  :: value
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) items(1:dict%Count) = dict%Items
        items(dict%Count+1)%Key   = key
        items(dict%Count+1)%Value = value

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + 1
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddEntry(dict, ent) Result (newdict)
        ! DESCRIPTION
        !   Adds a number to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)       :: dict
        Type(DICT_DP_ENTRY), Intent(IN) :: ent
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: err = 0
        Type(DICT_DP_ENTRY) :: items(dict%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) items(1:dict%Count) = dict%Items
        items(dict%Count+1) = ent

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + 1
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItems(dict, keys, values) Result (newdict)
        ! DESCRIPTION
        !   Adds some numbers to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)     :: dict
        Character(LEN=20), Intent(IN) :: keys(:)
        Integer, Intent(IN)  :: values(:)
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: i, err = 0
        Integer :: N = Min(Size(keys), Size(values))
        Type(DICT_DP_ENTRY) :: allitems(dict%Count+N)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) allitems(1:dict%Count) = dict%Items

        Do i = 1, N
            allitems(dict%Count+i)%Key   = keys(i)
            allitems(dict%Count+i)%Value = values(i)
        End Do

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + N
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = allitems
    End Function

    ! -------------------------------------------------------

    Function AddEntries(dict, entries) Result (newdict)
        ! DESCRIPTION
        !   Adds some numbers to the end of the dict and
        !   returns a new dict.

        ! ARGUMENTS.
        Type(DICT_DP), Intent(IN)       :: dict
        Type(DICT_DP_ENTRY), Intent(IN) :: entries(:)
        Type(DICT_DP) :: newdict

        ! VARIABLES.
        Integer :: i, err = 0
        Integer :: N = Size(entries)
        Type(DICT_DP_ENTRY) :: allitems(dict%Count+N)
        
        ! EXECUTABLE CODE.

        ! Save current dict.
        If (dict%Count > 0) allitems(1:dict%Count) = dict%Items
        allitems(dict%Count+1:) = entries

        ! Resize dict.        
        Deallocate(newdict%Items, STAT=err)
        newdict%Count = dict%Count + N
        Allocate(newdict%Items(newdict%Count), STAT=err)

        ! Assign back the dict.
        newdict%Items = allitems
    End Function

End Module

! =======================================================

Module DICT_T
    Use DICT_DP_T
    Use DICT_INT_T
    Implicit None
    Public
End Module


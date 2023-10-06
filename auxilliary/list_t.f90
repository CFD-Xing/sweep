! *****************************************************************************
!
! File:				list_t.f90
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
!   An implementation of a variable length list of different types.  Currently
!   included list types are:
!
!   > Double Precision (LIST_DP).
!   > Integer (LIST_INT).
!
! Functions:
!   Destroy     -   Clears memory associated with a list type.
!   SetP        -   Sets the pointer of ones list's items to another list's items.
!   Copy        -   Copies the items from one list to another.
!   Add         -   Adds two lists together and returns a new list (this is safe
!                   for returning one of the lists being added).
!   AddItem     -   Adds a single item to a list and returns a new list (this is safe
!                   for returning the RH list).
!   AddItems    -   Adds a multiple items to a list and returns a new list (this is
!                   safe for returning the RH list).
! *****************************************************************************

Module LIST_DP_T
    Implicit None
    Public

    Type LIST_DP
        Integer :: Count = 0
        Double Precision, Pointer :: Items(:)
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

    Subroutine Destroy(list)
        ! DESCRIPTION
        !   Clears memory associated with a list type.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(INOUT) :: list

        ! VARIABLES.
        Integer :: err = 0

        ! EXECUTABLE CODE.
        list%Count = 0
        Deallocate(list%Items, STAT=err)
        Nullify(list%Items)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine SetP(new, old)
        ! DESCRIPTION
        !   Assigns the list of one list type to another
        !   by passing the list pointer.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(INOUT) :: new
        Type(LIST_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items => old%Items
    End Subroutine

    ! -------------------------------------------------------

    Subroutine Copy(new, old)
        ! DESCRIPTION
        !   Copies the list of one list type to another.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(INOUT) :: new
        Type(LIST_DP), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items = old%Items
    End Subroutine

    ! -------------------------------------------------------

    Function Add(list1, list2) Result (sumlist)
        ! DESCRIPTION
        !   Adds together two list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(IN) :: list1
        Type(LIST_DP), Intent(IN) :: list2
        Type(LIST_DP) :: sumlist

        ! VARIABLES.
        Integer :: err = 0
        Double Precision :: items(list1%Count+list2%Count)
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list1%Count > 0) items(1:list1%Count) = list1%Items
        If (list2%Count > 0) items(list1%Count+1:) = list2%Items

        ! Resize list.        
        Deallocate(sumlist%Items, STAT=err)
        sumlist%Count = list1%Count + list2%Count
        Allocate(sumlist%Items(sumlist%Count), STAT=err)

        ! Assign back the list.
        sumlist%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItem(list, item) Result (newlist)
        ! DESCRIPTION
        !   Adds a number to the end of the list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(IN)    :: list
        Double Precision, Intent(IN) :: item
        Type(LIST_DP) :: newlist

        ! VARIABLES.
        Integer :: err = 0
        Double Precision :: items(list%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list%Count > 0) items(1:list%Count) = list%Items
        items(list%Count+1) = item

        ! Resize list.        
        Deallocate(newlist%Items, STAT=err)
        newlist%Count = list%Count + 1
        Allocate(newlist%Items(newlist%Count), STAT=err)

        ! Assign back the list.
        newlist%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItems(list, items) Result (newlist)
        ! DESCRIPTION
        !   Adds some numbers to the end of the list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_DP), Intent(IN)    :: list
        Double Precision, Intent(IN) :: items(:)
        Type(LIST_DP) :: newlist

        ! VARIABLES.
        Integer :: err = 0
        Double Precision :: allitems(list%Count+Size(items))
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list%Count > 0) allitems(1:list%Count) = list%Items
        allitems(list%Count+1:) = items

        ! Resize list.        
        Deallocate(newlist%Items, STAT=err)
        newlist%Count = list%Count + Size(items)
        Allocate(newlist%Items(newlist%Count), STAT=err)

        ! Assign back the list.
        newlist%Items = allitems
    End Function

    ! -------------------------------------------------------

    Double Precision Function SumList(list)
        ! DESCRIPTION
        !   Returns the sum of the list items.
        Type(LIST_DP), Intent(IN) :: list
        SumList = Sum(list%Items)
    End Function

End Module

! =======================================================

Module LIST_INT_T
    Implicit None
    Public

    Type LIST_INT
        Integer :: Count = 0
        Integer, Pointer :: Items(:)
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

    Subroutine Destroy(list)
        ! DESCRIPTION
        !   Clears memory associated with a list type.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(INOUT) :: list

        ! VARIABLES.
        Integer :: err = 0

        ! EXECUTABLE CODE.
        list%Count = 0
        Deallocate(list%Items, STAT=err)
        Nullify(list%Items)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine SetP(new, old)
        ! DESCRIPTION
        !   Assigns the list of one list type to another
        !   by passing the list pointer.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(INOUT) :: new
        Type(LIST_INT), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items => old%Items
    End Subroutine

    ! -------------------------------------------------------

    Subroutine Copy(new, old)
        ! DESCRIPTION
        !   Copies the list of one list type to another.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(INOUT) :: new
        Type(LIST_INT), Intent(IN)    :: old

        ! EXECUTABLE CODE.
        new%Count = old%Count
        new%Items = old%Items
    End Subroutine

    ! -------------------------------------------------------

    Function Add(list1, list2) Result (sumlist)
        ! DESCRIPTION
        !   Adds together two list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(IN) :: list1
        Type(LIST_INT), Intent(IN) :: list2
        Type(LIST_INT) :: sumlist

        ! VARIABLES.
        Integer :: err = 0
        Integer :: items(list1%Count+list2%Count)
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list1%Count > 0) items(1:list1%Count) = list1%Items
        If (list2%Count > 0) items(list1%Count+1:) = list2%Items

        ! Resize list.        
        Deallocate(sumlist%Items, STAT=err)
        sumlist%Count = list1%Count + list2%Count
        Allocate(sumlist%Items(sumlist%Count), STAT=err)

        ! Assign back the list.
        sumlist%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItem(list, item) Result (newlist)
        ! DESCRIPTION
        !   Adds a number to the end of the list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(IN)    :: list
        Integer, Intent(IN) :: item
        Type(LIST_INT) :: newlist

        ! VARIABLES.
        Integer :: err = 0
        Integer :: items(list%Count+1)
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list%Count > 0) items(1:list%Count) = list%Items
        items(list%Count+1) = item

        ! Resize list.        
        Deallocate(newlist%Items, STAT=err)
        newlist%Count = list%Count + 1
        Allocate(newlist%Items(newlist%Count), STAT=err)

        ! Assign back the list.
        newlist%Items = items
    End Function

    ! -------------------------------------------------------

    Function AddItems(list, items) Result (newlist)
        ! DESCRIPTION
        !   Adds some numbers to the end of the list and
        !   returns a new list.

        ! ARGUMENTS.
        Type(LIST_INT), Intent(IN)    :: list
        Integer, Intent(IN) :: items(:)
        Type(LIST_INT) :: newlist

        ! VARIABLES.
        Integer :: err = 0
        Integer :: allitems(list%Count+Size(items))
        
        ! EXECUTABLE CODE.

        ! Save current list.
        If (list%Count > 0) allitems(1:list%Count) = list%Items
        allitems(list%Count+1:) = items

        ! Resize list.        
        Deallocate(newlist%Items, STAT=err)
        newlist%Count = list%Count + Size(items)
        Allocate(newlist%Items(newlist%Count), STAT=err)

        ! Assign back the list.
        newlist%Items = allitems
    End Function

    ! -------------------------------------------------------

    Integer Function SumList(list)
        ! DESCRIPTION
        !   Returns the sum of the list items.
        Type(LIST_INT), Intent(IN) :: list
        SumList = Sum(list%Items)
    End Function

End Module

! =======================================================

Module LIST_T
    Use LIST_DP_T
    Use LIST_INT_T
    Implicit None
    Public
End Module


! *****************************************************************************
!
! File:				swpensemble.f90
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
!   The binary tree is an efficient method of storing structured data, and
!   is used in Sweep as a way of viewing the particle ensemble in a way that
!   is convenient for particle selection by different properties.
!
!   The binary tree is designed to allow the random selection of particles
!   according to arbitrary distributions speicified only by the weights 
!   assigned to individual particles.  The FORTRAN representation of the 
!   particles is given in the module Particles.  In this module the most 
!   important thing is that the type representing a stochastic particle 
!   contains a field (which may also be a derived type) called "properties", 
!   this will be some array of numbers which are the weights mentioned above.
!
!   This implementation uses a binary tree which implememnts particle doubling 
!   and maintains a sample volume - ie a scale factor.  Arguably in this 
!   respect I have extended the functionality to cover an issue that should
!   be handled at a higher level in the program.
!
!   The storage capcity of the tree is 2**levels particles but for convenience
!   (in the implementation of add_particle) we never store a particle is the 
!   last position so in fact the largest number of particles ever in the tree
!   is 2**levels - 1
!
!   Integer types are used as indexes into the tree and list but callers are 
!   not supposed to use this except to make appropriate variable declarations.
!   In particular callers should only get index values from choose_particle().
!
!   Multiple trees for different particle ensembles are not supported.
!
!   Usage must begin with a call to bintree_prepare(levels) to create a 
!   tree with the desired number of levels.  To release the memory call
!   bintree_destroy.
!
! Functions:
!   Init			-	Initialises the binary tree to hold at least the
!						given number of particles.  This routine must be
!						called before the tree is used.
!   Destroy			-	Clears all dynamic memory associated with the tree
!						and resets all variables.
!   AscendingRecalc		-	Recalculates the values in the tree based on particle
!						properties, starting at the given tree node.
!   AddParticle			-	Adds a particle to the ensemble, contracting the
!						population if required.
!   ReplaceParticle		-	Replaces a particle currently in the tree with the
!						given particle.
!   ChooseParticle		-	Takes advantage of the binary tree architecture to
!						select a particle using the given weighting method.
!   CustomChooseParticle	-	Selects a particle from the list using a weight
!						which is not stored in the tree.
!   RemoveParticle		-	Removes a particle from the tree.  Also runs the
!						particle doubling algorithm if required.
!   RemoveParticleND	        -	Internal function only.  Removes a particle from the
!						list without running the doubling algorithm.
!   FreeSpace			-	Returns the number of free spaces in the particle list.
!   GetParticleSums		-	Returns the sums of soot properties stored in the
!						tree.  As these values are calculated in order to
!						build the tree, they do not to be summed using the
!						actual particle list.
!   CustomParticleSum	        -	Returns the sum of a custom particle property over
!						all particles.  This routine cannot take advantage
!						of the binary tree architecture.
!   GetParticle			-	Gets the particle at the given index in the list.
!   ParticleCount		-	Returns the number of stochastic particles in the
!						list.
!   TypedParticleCount	        -	Returns the number of stochastic particles of a 
!						given type in the list.
!   MaxParticleCount	        -	Returns the capacity of the binary tree.
!   ScalingFactor		-	Returns the current scaling factor of the binary
!						tree due to tree operations.
!   ClearEnsemble		-	Removes all particles from the list and resets the
!						binary tree.
!   GetEnsembleArray	        -	Returns a 2D array of all particle properties.  Only
!						the unique properties are given (count given by DIM
!						parameter in SWPPART).
!   GetTypedList		-	Returns the list of particles using the StochParticle
!						derived type.
!   LoadEnsembleArray	        -	Builds the particle list given a 2D of particle
!						properties.
!   SetTypedList		-	Sets the particle list given a list of StochParticle
!						derived types.
!   GetStats			-	Gets the statistical items summed over all particles.
!   TreeIndex			-	Returns the index of a particle in the tree given the
!						index of the particle in the particle list.
!   GetScalingCounters	        -	Returns the scaling counters (doubling, expansions, 
!						contractions).
!   GetCompositionMatrix	-	Returns the composition of all particles in a
!						VERY LARGE matrix.
!   GetAverageComposition	-	Returns the average composition of all particles
!						in the ensmeble.
! *****************************************************************************

Module SWPENSEMBLE
    ! -------------------------------------------------------
    ! IMPORT PUBLIC MODULE FUNCTIONS.
    ! -------------------------------------------------------
    Use SWPPART, GetStatItems=>GetStats ! Definition of a stochastic particle.
    ! -------------------------------------------------------

    Implicit None
    Public
    Private :: AscendingRecalc, RemoveParticleND

    ! -------------------------------------------------------
    ! BINARY TREE PARAMETERS.
    ! -------------------------------------------------------

    ! Default number of tree levels, if it is not specified
    ! by the calling program.
    Integer, Parameter, Private :: DEFAULT_LEVELS = 10

    ! Flag to turn on console messages from the binary tree.
    Logical, Parameter, Private :: MSGS = .True.

    ! Switch for the doubling mechanism.
    Logical, Parameter, Private :: USE_DOUBLING = .True.

    ! -------------------------------------------------------
    ! DEF'N OF A BINARY TREE NODE.
    ! -------------------------------------------------------

    Type, Private :: TreeNode
        Real :: Left(PROPERTY_COUNT)
        Real :: Right(PROPERTY_COUNT)
    End Type

    Type, Private :: BinaryTree
        Integer :: Capacity, Levels ! Number of levels in tree.
        Type(TreeNode), Pointer :: Nodes(:) ! Binary tree nodes.  
    End Type

    Type :: Ensemble
        Type(StochParticle), Pointer :: Particles(:) ! List of particles.
        Type(BinaryTree) :: Tree         ! Binary tree structure.
        Integer :: Capacity              ! Max. number of particles which may be stored. 
        Real    :: Scaling               ! Scaling factor based on ensemble operations.
        Real    :: contFactor, expFactor ! Constant factors for contraction and expansion.
        Integer :: NEXP, NCONT, NDBLE    ! Expansion, contraction and doubling counts.
        Logical :: CONTWARN              ! Has a contraction warning been issued.
        Integer :: FirstSpace            ! First free space in particle list.
        Integer :: HalfCapacity          ! Half the capacity, precalculated for convenience.
        ! Doubling algorithm variables.
        Logical :: ActiveDoubling ! Is doubling active?
        Integer :: DoublingCutOff ! Minimum particle count at which doubling is activated.
        Integer :: DoublingLimit  ! Particle count at which ensemble is doubled.
        Integer :: DoublingSlack  ! Slack space in particle list after doubling.
    End Type

    Contains

    ! -------------------------------------------------------

    Subroutine Init(ens, N, flag)
        ! DESCRIPTION:
        !   Initialises a particle ensemble.  This function
        !   must be called before the ensemble can be used.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(OUT) :: ens  ! Ensemble to initialise.
        Integer, Intent(IN)         :: N    ! Number of particles to enter into tree.  Must
                                            ! be a power of two.  If not a power of two, then
                                            ! the first power of two larger than N is used.
        Integer, Intent(OUT)        :: flag ! Error flag.

        ! VARIABLES.
        Integer :: err
	Real    :: RL

        ! EXECUTABLE CODE.
        flag = 0

        ! Determine the number of tree levels
        ! required to hold at least N particles.
        RL = Log(Real(N)) / Log(2.0E0)
        ens%Tree%Levels  = NInt(RL)
        If (Real(ens%Tree%Levels) < RL) Then
            ens%Tree%Levels = ens%Tree%Levels + 1
            flag = 1 ! Tell program ensemble size has been changed.
        End If

        ! Save ensemble capacity.
        ens%Capacity = 2**ens%Tree%Levels
        ens%Tree%Capacity = ens%Capacity
        ens%HalfCapacity = ens%Capacity / 2

        ! Initialise ensemble scaling variables.
        ens%Scaling        = 1.0
        ens%contFactor     = Real(ens%Capacity-1) / Real(ens%Capacity)
        ens%expFactor      = Real(ens%Capacity - 1) / Real(ens%Capacity - 2)
        ens%NEXP           = 0
        ens%NCONT          = 0
        ens%NDBLE          = 0
        ens%CONTWARN       = .False.
        ens%ActiveDoubling = .False.
        ens%DoublingLimit  = ens%HalfCapacity - 2**Max(ens%Tree%Levels-5,0)
        ens%DoublingSlack  = 2**Max(ens%Tree%Levels-4,0)
        ens%DoublingCutOff = 3 * ens%Capacity / 4

        ! Allocate memory for particle list and binary tree.
        Allocate(ens%Particles(ens%Capacity), STAT=err)
        Allocate(ens%Tree%Nodes(ens%Capacity-1), STAT=err)

        ! Reset list and tree memory.
        Call ClearEnsemble(ens)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine CopyEnsemble(e1, e2)
        ! DESCRIPTION:
        !   Copies the data from one ensemble to another.

        Implicit None
        
        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: e1, e2

        ! VARIABLES.
        Integer :: err

        ! EXECUTABLE CODE.
        e2%Tree%Capacity = e1%Tree%Capacity
        e2%Tree%Levels = e1%Tree%Levels
        e2%Tree%Nodes = e1%Tree%Nodes
        e2%Particles = e1%Particles
        e2%Scaling = e1%Scaling
        e2%contFactor = e1%ContFactor
        e2%expFactor = e1%expFactor
        e2%NEXP = e1%NEXP
        e2%NCONT = e1%NCONT
        e2%NDBLE = e1%NDBLE
        e2%CONTWARN = e1%CONTWARN
        e2%FirstSpace = e1%FirstSpace
        e2%ActiveDoubling = e1%ActiveDoubling
    End Subroutine

    ! -------------------------------------------------------

    Subroutine DeleteEnsemble(ens)
        ! DESCRIPTION:
        !   Deallocates memory used for an ensemble.
        Implicit None
        
        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: ens ! Ensemble to destroy.

        ! VARIABLES.
        Integer :: err

        ! EXECUTABLE CODE.
        Deallocate(ens%Particles, ens%Tree%Nodes, STAT=err)
        ens%Tree%Levels    = 0
        ens%Capacity       = 0
        ens%halfCapacity   = 0
        ens%Scaling        = 0.0
        ens%contFactor     = 0.0
        ens%expFactor      = 0.0
        ens%NEXP           = 0
        ens%NCONT          = 0
        ens%NDBLE          = 0
        ens%CONTWARN       = .False.
        ens%ActiveDoubling = .False.
        ens%DoublingLimit  = 0
        ens%DoublingSlack  = 0
        ens%DoublingCutOff = 0
    End Subroutine

    ! -------------------------------------------------------

    Subroutine AscendingRecalc(tree, start, flag)
        ! DESCRIPTION:
        !   Works its way up the tree from a given node
        !   updating the accumulated values, definitely not
        !   for users to call.
        !
        !   The first argument is the node number to start from, 
        !   NB root is node 1, its two children are 2 and 3,
        !   the children of node n are 2n and 2n+1.

        Implicit None

        ! ARGUMENTS.
        Type(BinaryTree), Intent(INOUT) :: tree     ! Tree to recalc.
        Integer, Intent(IN)             :: start    ! Particle index from which to start.
        Integer, Intent(OUT)            :: flag     ! Error flag.

        ! VARIABLES.
        Integer :: j, k

        ! EXECUTABLE CODE.
        flag = 0

        ! Check start has a sensible value.
        If (1 < start .And. start < tree%Capacity) Then
            j = start

            Do While (j > 1)
                ! These methods exploit the truncation
                ! properties of integer division.

                k = j / 2

                ! Work out whether we are to the left
                ! or right of the parent.
                If (2*k == j) Then
                    tree%Nodes(k)%Left = tree%Nodes(j)%Left + tree%Nodes(j)%Right
                Else
                    tree%Nodes(k)%Right = tree%Nodes(j)%Left + tree%Nodes(j)%Right
                End If

                ! Move up tree to parent.
                j = j / 2
           End Do
        Else
           ! Invalid value of start.
           flag = -1
        End If
    End Subroutine

    ! -------------------------------------------------------

    Integer Function AddParticle(ens, sp)
        ! DESCRIPTION:
        !   Stores sp in a leaf of the tree, updates the tree,
        !   returns a positive integer indicating the position of
        !   the particle in the list of particles on success or
        !   a negative integer on failure.
        !
        !   The function returns a positive integer on success
        !   (this should be the position of the new	entry in the
        !   internal list but the particle may be moved so callers
        !   should not use the precise value).  A negative value
        !   indicates an error.
        !
        ! RETURNS:
        !   Index of particle in list is successful, otherwise
        !   < 0.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT)   :: ens ! Ensemble in which to put particle.
        Type(StochParticle), Intent(IN) :: sp  ! Particle to add to ensemble.

        ! VARIABLES.
        Integer :: i, dest, err

        ! EXECUTABLE CODE.
        AddParticle = 0

        If (ens%FirstSpace < ens%Capacity) Then
            ! There is space to put new particle
            ! at end of list.
            dest = ens%FirstSpace
            ens%FirstSpace = ens%FirstSpace + 1
        Else
            ! Uniformly select a particle to remove
            ! and rescale the ensemble.
            dest = IRnd(ens%Capacity) ! Random integer in [1, 2**levels].

            ! NOTE:
            ! If dest == 2**levels then the new particle
            ! is selected for removal.

            ! Record the rescaling of the sample volume (removal of a
            ! particle is a contraction).
            Call NoteContraction(ens, dest, .True.)
        End If

        If (dest < ens%Capacity) Then
            ! We are have not selected the
            ! new particle for deletion.

            ! Put the particle into the list.
            ens%Particles(dest) = sp

            ! NOTE:
            ! This section exploits integer division:

            ! First work out which tree node this
            ! particle is on.
            i = TreeIndex(dest, ens)

            ! Now work out whether this entry is to the left
            ! or the right of the node i and make an entry
            ! in the node accordingly.
            If (Mod(dest, 2) == 1) Then
               ens%Tree%Nodes(i)%Left = GetPreCalcs(sp)
            Else
               ens%Tree%Nodes(i)%Right = GetPreCalcs(sp)
            End If

            ! Update the cumulative values in the tree.
            Call AscendingRecalc(ens%Tree, i, err)

            If (err == 0) Then
                ! Success!
                AddParticle = dest
            Else
                ! Failed to update tree.
                AddParticle = -1
            End If
        Else
            ! No need to do anything since the new
            ! particle is immediately forgotten.
            AddParticle = ens%Capacity + 1 ! (not a normal return value).
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function ReplaceParticle(ens, dest, sp)
        ! DESCRIPTION:
        !    Replaces the particles at the position given
        !    in the first argument with the particle passed
        !    as the second argument.  Returns true on success,
        !    false if a problem is encountered.  Replacing a
        !    non-existant entry is an error.
        !
        ! RETURNS:
        !    < 0 if not successful.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT)   :: ens  ! Ensemble in which to replace particle.
        Integer, Intent(IN)             :: dest ! The index of the particle to be replaced.
        Type(StochParticle), Intent(IN) :: sp   ! Particle to add to ensemble.

        ! VARIABLES.
        Integer :: i, err

        ! EXECUTABLE CODE.
        ReplaceParticle = 0

        If (0 < dest .And. dest < ens%FirstSpace) Then
            ! This is a valid index to replace.

            ! Replace the particle.
            ens%Particles(dest) = sp

            ! Calculate index of tree leaf referring
            ! to this particle.
            i = TreeIndex(dest, ens)

            ! and update this leaf
            If (Mod(dest, 2) == 1) Then
                ens%Tree%Nodes(i)%left = GetPreCalcs(sp)
            Else
                ens%Tree%Nodes(i)%right = GetPreCalcs(sp)
            End If

            !propagate the change up the tree
            Call AscendingRecalc(ens%Tree, i, err)

            If (err < 0) Then
                ReplaceParticle = -1
            End If
        Else
            ! Invalid particle index.
            ReplaceParticle = -2
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function ChooseParticle(ens, wt)
        ! DESCRIPTION:
        !   Randomly selects a particle according the weights
        !   given in the wt position in the rate arrays
        !   and returns the position of the selected particle
        !   in the list.
        !
        !   This value is guaranteed to be valid until the next
        !   time a particle is removed.  A negative value
        !   indicates some kind of error.  Never returns 0.
        !
        ! RETURNS:
        !   Index of a particle in the particle list, <0 if
        !   not successful.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens ! Ensemble from which to choose particle.
        Integer, Intent(IN)        :: wt  ! Index of the particle property to use as a 
                                          ! weight for particle selection.

        ! VARIABLES.
        Double Precision :: rand, rand2
        Integer          :: j, k

        ! EXECUTABLE CODE.
        ChooseParticle = 0

        If (wt >= 1 .And. wt <= PROPERTY_COUNT) Then
            ! We can choose using the tree.

            ! Scale random number up according to total
            ! weight of all particles performing this
            ! procedure in two stages is convenient during
            ! debugging, in a release build I expect the
            ! compiler to optimise the temporary rand2 away.
            rand2 = Rnd()
            rand  = ens%Tree%Nodes(1)%Left(wt) + ens%Tree%Nodes(1)%Right(wt)
            rand  = rand * rand2

            j = 1

            Do k = 1, ens%Tree%Levels
                If (rand <= ens%Tree%Nodes(j)%Left(wt)) Then
                    ! Go left
                    j = 2*j
                ElseIf (ens%Tree%Nodes(j)%Right(wt) > 0) Then
                    ! Go right (subtract weight of left
                    ! subtree now it has been rejected).
                    rand = rand - ens%Tree%Nodes(j)%Left(wt) 
                    j = 2*j + 1
                Else 
                    ! The accumulated error in rand that is inevitable
                    ! with finite precision arithmetic is trying to
                    ! send us down the right subtree which has nothing
                    ! in it.  The machine has the following equality
                    ! as only approximate rand == tree(j)%left(rate_id) so 
                    ! the obvious solution is to make it as exact as
                    ! possible and go left.
                    rand = ens%Tree%Nodes(j)%Left(wt)
                    j = 2*j
                End If
            End Do

            ! This DO loop updates j after reaching a leaf
            ! of the tree so j is no longer a valid index
            ! for tree.
            ChooseParticle = j - ens%Capacity + 1
        Else
            If (wt == 0) Then
                ! Select uniformly.

                ChooseParticle = ens%FirstSpace

                ! Get a random integer, assuming all bits
                ! are indep taking value 1 with prob 1/2 then
                ! we use a mask to reduce to an integer uniformly
                ! distributed in [0, 2**levels - 1] and then add
                ! 1 to get a uniformly distributed list index the
                ! selected position may be empty so repeat the
                ! procedure until we get a choice < first_space.

                Do While (ChooseParticle >= ens%FirstSpace) ! Potential infinite loop
                    ChooseParticle = IRnd(ens%Capacity)
                End Do
            Else
                ! Invalid rate index specified.
                ChooseParticle = -1
            End If
        End If

        If (ChooseParticle >= ens%FirstSpace) Then
            ! Chosen a non existant particle.
            ChooseParticle = -2
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function ChooseTypedParticle(ens, wt, typeid)
        ! DESCRIPTION:
        !   Randomly selects a particle according the weights
        !   given in the wt position in the rate arrays
        !   and returns the position of the selected particle
        !   in the list.  The particle will have the given
        !   type.
        !
        !   This value is guaranteed to be valid until the next
        !   time a particle is removed.  A negative value
        !   indicates some kind of error.  Never returns 0.
        ! RETURNS:
        !   Index of a particle in the particle list, <0 if
        !   not successful.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens    ! Ensemble from which to choose particle.
        Integer, Intent(IN)        :: wt     ! Index of the particle property to use as a 
                                             ! weight for particle selection.
        Integer, Intent(IN)        :: typeid ! The type of particle to return.

        ! VARIABLES.
        Integer :: N = 0

        ! EXECUTABLE CODE.
        N = TypedParticleCount(ens, typeid)
        If (N > 0) Then
            ChooseTypedParticle = ChooseParticle(ens, wt)
            Do While(ens%Particles(ChooseTypedParticle)%TypeID /= typeid)
                ChooseTypedParticle = ChooseParticle(ens, wt)
            End Do
        Else
            ChooseTypedParticle = -1
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function CustomChooseParticle(ens, masswt, surfwt, actsurfwt, diamwt)
        ! DESCRIPTION:
        !   Randomly selects a particle according to a custom
        !   particle weight given by the exponents of
        !   standard particle properties.
        !
        !   This value is guaranteed to be valid until the next
        !   time a particle is removed.  A negative value
        !   indicates some kind of error.  Never returns 0.
        !
        ! RETURNS:
        !   Index of a particle in the particle list, <0 if
        !   not successful.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens ! Ensemble from which to choose particle.
        Real, Intent(IN) :: masswt, &     ! Mass/volume exponent.
                            surfwt, &     ! Surface area exponent.
                            actsurfwt, &  ! Active surface exponent.
                            diamwt        ! Diameter exponent.

        ! VARIABLES.
        Integer :: i, N
        Real    :: rand, wts(ens%FirstSpace - 1)

        ! EXECUTABLE CODE.

        ! Find the number of particles.
        N = ens%FirstSpace - 1

        ! Calculate the custom weight of each particle.
        Do i = 1, N
           wts(i) = CalcParticleWeight(ens%Particles(i), masswt, surfwt, actsurfwt, diamwt)
        End Do

        ! Generate a uniform random number scaled by the
        ! sum of all weights.
        rand = Rnd() * Sum(wts)

        ! Select a particle using DIV.
        CustomChooseParticle = 1
        rand = rand - wts(1)
        Do While (rand > 0.0E0 .And. CustomChooseParticle < N)
            CustomChooseParticle = CustomChooseParticle + 1
            rand = rand - wts(CustomChooseParticle)
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function RemoveParticle(ens, id)
        ! DESCRIPTION:
        !   Removes a particle from the list and updates the
        !   tree, this may invalidate any pre-existing
        !   index values for the tree and list.
        !
        ! RETURNS:
        !   <0 if not successful.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: ens ! Ensmeble from which to remove particle.
        Integer, Intent(IN)           :: id  ! Index of the particle to remove.

        ! VARIABLES.
        Integer :: i, j, k, err

        ! EXECUTABLE CODE.
        RemoveParticle = 0

        If (USE_DOUBLING) Then
            If ((.Not. ens%ActiveDoubling) .And. (ens%FirstSpace >= ens%DoublingCutOff)) Then
                ens%ActiveDoubling = .True.
                If (MSGS) Print *, "Sweep:\> Particle doubling activated."
            End If
        End If

        If (1 <= id .And. id < ens%Capacity .And. id < ens%FirstSpace) Then !argument in range
           If (ens%FirstSpace == ens%DoublingLimit + 1 .And. ens%ActiveDoubling) Then 
                ! Particle ensemble should be doubled after the removal.

                ! First remove the particle.
                err = RemoveParticleND(ens, id)

                If (err == 0) Then
                ! Initial removal successful.

                    ! Now perform the doubling:
                    If (MSGS) Print *, "Sweep:\> Doubling particle ensemble."

                    ! Copy the particles.
                    ens%Particles(ens%DoublingLimit:ens%Capacity-ens%DoublingSlack-2) = &
                    ens%Particles(1:ens%DoublingLimit-1)

                    ! Load particles into the bottom tree level.
                    Do j = 0, ens%DoublingLimit-2
                        ens%Tree%Nodes(ens%HalfCapacity+j)%Left  = GetPreCalcs(ens%Particles(2*j+1))
                        ens%Tree%Nodes(ens%HalfCapacity+j)%Right = GetPreCalcs(ens%Particles(2*j+2))
                    End Do

                    ! Update all tree levels.
                    Do i = ens%Tree%Levels - 1, 1, -1
                        Do j=0, 2**(i-1) - 1
                            k = 2**(i-1) + j
                            ens%Tree%Nodes(k)%Left  = ens%Tree%Nodes(2*k)%Left   + ens%Tree%Nodes(2*k)%Right
                            ens%Tree%Nodes(k)%Right = ens%Tree%Nodes(2*k+1)%left + ens%Tree%Nodes(2*k+1)%Right
                        End Do
                    End Do

                    ! Update scaling.
                    ens%FirstSpace = ens%Capacity - ens%DoublingSlack - 1
                    ens%NDBLE = ens%NDBLE + 1

                Else 
                    ! Initial particle removal failed.
                    ! Could not remove particle prior to doubling the ensemble.
                    RemoveParticle = -1
                End If

            Else
               ! Do not need to double.
               err = RemoveParticleND(ens, id)
               If (err < 0) RemoveParticle = -2
            End If
        Else
           ! Tried to remove from an invalid position.
           RemoveParticle = -3
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function RemoveParticleND(ens, id)
        ! DESCRIPTION:
        !   This is a private helper function, it removes a
        !   particle from the list and updates the tree
        !   without any particle doubling, this may invalidate
        !   any pre-existing index values for the tree and list.
        !
        !   The key idea is that we overwrite the particle to
        !   be removed with the one at the end of list so that
        !   we do not create a gap in the list.
        !
        ! RETURNS:
        !    <0 if not successful.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: ens ! Ensmeble from which to remove particle.
        Integer, Intent(IN)           :: id  ! Index of the particle to remove.

        ! VARIABLES.
        Integer :: j, k, err

        ! EXECUTABLE CODE.
        RemoveParticleND = 0

        If (1 <= id .And. id < ens%FirstSpace) Then
            ! We have a sensible index to remove.

            ! NOTE:
            ! This code works even if id == FirstSpace - 1.

            ! Overwrite removed particle and wipe particle at
            ! end of list.
            ens%Particles(id) = ens%Particles(ens%FirstSpace - 1) 
            ens%Particles(ens%FirstSpace - 1)%Properties = 0

            ! Calculate the tree leaf above FirstSpace - 1.
            j = TreeIndex(ens%FirstSpace - 1, ens)

            ! Put 0 entries into the tree leaf above the
            ! newly empty position (FirstSpace - 1 in list).
            If (Mod(ens%FirstSpace, 2) == 0) Then
                ! FirstSpace - 1 is to the left of leaf j.
                ens%Tree%Nodes(j)%Left  = 0
            Else
                ! to the right.
                ens%Tree%Nodes(j)%Right = 0
            End If

            ! Calculate the tree leaf above id (where the entry has changed)
            k = TreeIndex(id, ens)

            ! Put the new property entries into the
            ! tree leaf above id.
            If(Mod(id, 2) == 1) Then
                ! id is to the left of leaf k.
                ens%Tree%Nodes(k)%Left  = GetPreCalcs(ens%Particles(id))
            Else
                ! to the right.
                ens%Tree%Nodes(k)%Right = GetPreCalcs(ens%Particles(id))
            End If

            ! Now propagate the changes to the tree leaves upwards.
            Call AscendingRecalc(ens%Tree, j, err)
            If (err < 0) RemoveParticleND = -1
            Call AscendingRecalc(ens%Tree, k, err)
            If (err < 0) RemoveParticleND = -1
            ! The first space moves back by 1 now a
            ! particle has been removed.
            ens%FirstSpace = ens%FirstSpace - 1
        Else
            ! Non-existant particle.
            RemoveParticleND = -2
        End If
    End Function

    ! -------------------------------------------------------

    Integer Function FreeSpace(ens)
        ! DESCRIPTION:
        !   Returns the space available in the ensemble for
        !   additional particles.
        ! RETURNS:
        !   Free space in ensemble.
        Implicit None
        Type(Ensemble), Intent(IN) :: ens
        FreeSpace = ens%Capacity - ens%FirstSpace + 1
    End Function

    ! -------------------------------------------------------

    Function GetParticleSums(ens)
        ! DESCRIPTION:
        !   Returns an array of particle property totals.
        ! RETURNS:
        !   Particle sums.
        Implicit None
        Real                       :: GetParticleSums(PROPERTY_COUNT)
        Type(Ensemble), Intent(IN) :: ens
        GetParticleSums = ens%Tree%Nodes(1)%Left + ens%Tree%Nodes(1)%Right
    End Function

    ! -------------------------------------------------------

    Real Function CustomParticleSum(ens, masswt, surfwt, actsurfwt, diamwt)
        ! DESCRIPTION:
        !   Returns an array of particle property totals.
        ! RETURNS:
        !   Particle sums.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens ! Ensemble to calculate for.
        Real, Intent(IN) :: masswt, &     ! Mass/volume exponent.
                            surfwt, &     ! Surface area exponent.
                            actsurfwt, &  ! Active surface exponent.
                            diamwt        ! Diameter exponent.

        ! VARIABLES.
        Integer          :: i, N

        ! EXECUTABLE CODE.

        ! Find the number of particles.
        N = ens%FirstSpace - 1

        ! Calculate the custom weight of each particle.
        Do i = 1, N
            CustomParticleSum = CustomParticleSum + &
                                CalcParticleWeight(ens%Particles(i), masswt, surfwt, actsurfwt, diamwt)
        End Do
    End Function

    ! -------------------------------------------------------

    Integer Function GetParticle(ens, i, sp)
        ! DESCRIPTION:
        !   Put the particle at position i in the list into
        !   the second argument.
        ! RETURNS:
        !   <0 if not successful.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens  ! Ensemble from which to retrieve particle.
        Integer, Intent(IN)        :: i    ! Index of the particle to get.
        Type(StochParticle), Intent(OUT) :: sp ! Particle.

        ! EXECUTABLE CODE.

        If (0 < i .And. i < ens%FirstSpace) Then
            sp = ens%Particles(i)
            GetParticle = 0
        Else
            GetParticle = -1
        End If
    End Function

    ! -------------------------------------------------------

    Pure Integer Function ParticleCount(ens)
        ! DESCRIPTION:
        !   Returns the current number of particles in the
        !   ensemble.
        !
        ! RETURNS:
        !   Particle count.
        Implicit None
        Type(Ensemble), Intent(IN) :: ens
        ParticleCount = ens%FirstSpace - 1
    End Function

    ! -------------------------------------------------------

    Pure Integer Function TypedParticleCount(ens, typeID)
        ! DESCRIPTION:
        !   Returns the current number of particles in the
        !   ensemble of the given type.
        ! RETURNS:
        !   Count of particles of given type.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens    ! Ensemble to look in.
        Integer, Intent(IN)        :: typeID ! Type to count.

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        TypedParticleCount = 0
        Do i = 1, ens%FirstSpace - 1
            If (ens%Particles(i)%TypeID == typeID) TypedParticleCount = TypedParticleCount + 1
        End Do
    End Function

    ! -------------------------------------------------------

    Pure Real function ScalingFactor(ens)
        ! DESCRIPTION:
        !    Returns the current scaling factor for the
        !    ensemble.
        ! RETURNS:
        !    Scaling factor.
        Implicit None
        Type(Ensemble), Intent(IN) :: ens
        ScalingFactor = ens%Scaling * &
                            ens%contFactor**ens%NCONT * &
                            Real(2**ens%NDBLE)
    End function

    ! -------------------------------------------------------

    Subroutine NoteContraction(ens, i, disp)
        ! DESCRIPTION:
        !   Updates the ensemble NCONT as a contraction is
        !   performed.
        Implicit None
        Type(Ensemble), Intent(INOUT) :: ens  ! Ensemble to check.
        Integer, Intent(IN)           :: i    ! Particle removed.
        Logical, Intent(IN)           :: disp ! Print contraction to screen.

        If ((Real(ens%NCONT)/Real(ens%Capacity) > 0.01E0) .And. .Not. ens%CONTWARN .And. disp) Then
           Print *, "Sweep:\> Ensemble contracting too often: Possible stiffness issue."
           ens%CONTWARN = .True.
        End If
        ens%NCONT = ens%NCONT + 1
    End Subroutine

    ! -------------------------------------------------------

    Subroutine ClearEnsemble(ens)
        ! DESCRIPTION:
        !   Clears the current ensemble by removing all
        !   particles and setting counters to default
        !   values.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: ens

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        ! Set all nodes and particles to 0.
        Do i = 1, ens%Capacity - 1
            ens%Tree%Nodes(i)%Left  = 0.0E0
            ens%Tree%Nodes(i)%Right = 0.0E0
            ens%Particles(i)%Properties = 0.0E0
        End Do

        ! Tree is completely empty now
        ens%FirstSpace = 1

        ! Nothing has happened yet.
        ens%Scaling        = 1.0
        ens%NEXP           = 0
        ens%NCONT          = 0
        ens%NDBLE          = 0
        ens%CONTWARN       = .False.
        ens%ActiveDoubling = .False.
    End Subroutine

    ! -------------------------------------------------------

    Subroutine GetEnsembleArray(ens, N, p)
        ! DESCRIPTION:
        !   Put the particles in the tree into a two
        !   dimensional array.  The argument should be an
        !   unassociated array pointer to have the appropriate
        !   space allocated to it - the caller is responsible
        !   for later deallocation.  The first dimension of 
        !   the array (the fastest varying subscript) will have
        !   dimension to hold all the real numbers needed to
        !   describe a particle.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens    ! Ensemble to get.
        Integer, Intent(OUT)       :: N      ! Return the number of particles in the list.
        Real, Intent(OUT)          :: p(:,:) ! Returns the particle list

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        p = 0.0E0
        N = ParticleCount(ens)

        ! Loop over all the particles, putting
        ! them into the array.
        Do i = 1, Min(Size(p,2), N)
           p(:,i) = ens%Particles(i)
        End Do
    End Subroutine

    ! -------------------------------------------------------

    Subroutine GetTypedList(ens, N, p)
        ! DESCRIPTION:
        !   Copies the particles currently in the
        !   ensemble into a second list.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN)       :: ens  ! Ensemble to get.
        Integer, Intent(OUT)             :: N    ! Return the number of particles in the list
        Type(StochParticle), Intent(OUT) :: p(:) ! Returns the particle list.

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.
        N = ParticleCount(ens)

        ! Loop over all the particles, putting
        ! them into the array.
        Do i = 1, Min(Size(p), N)
           p(i) = ens%Particles(i)
        End Do
    End Subroutine

    ! -------------------------------------------------------

    Subroutine LoadEnsembleArray(N, p, ens, flag)
        ! DESCRIPTION:
        !   Read particles into the tree from a two
        !   dimensional array.  The first dimension of the
        !   array (the fastest varying subscript) will have
        !   dimension to hold all the real numbers needed to
        !   describe a particle.  vol is used to return the
        !   sample volume which the tree/list describes.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)           :: N      ! Number of particles to be added
	Real, Intent(IN)              :: p(:,:) ! Particles to be added to the ensemble.
        Type(Ensemble), Intent(INOUT) :: ens    ! Ensemble to get.
        Integer, Intent(OUT)          :: flag   ! Error flag.

        ! VARIABLES.
        Integer :: i, m, err
        Logical :: save_dble

        ! EXECUTABLE CODE.

        ! Empty the tree, but keep the ActiveDoubling flag
        ! so that doubling can continue.  This shows my 
        ! data flow and control is a bit mixed up since we
        ! treat ActiveDoubling differently to the
        ! sample volume / scaling.
        save_dble = ens%ActiveDoubling
        Call ClearEnsemble(ens)
        ens%ActiveDoubling = save_dble

        If (N > ens%Capacity - 1) Then
            ! Tree is not big enough.
            flag = -1
        ElseIf (N <= Size(p, 2) .And. N >= 0) Then 
            ! N is a reasonable value.

            ! Number of times we can fit the list
            ! into the tree.
            m = (ens%Capacity - ens%DoublingSlack) / N
            If (.Not. ens%ActiveDoubling) Then
                m = 1
            End If
            If (N > ens%DoublingCutOff) Then
                ens%ActiveDoubling = USE_DOUBLING
            End If

            ! Put the particles into the list of
            ! tree contents.
            Do i = 1, N
                ens%Particles((i-1)*m+1:i*m) = p(:,i)
            End Do

            ens%FirstSpace = N * m + 1

            If (ens%FirstSpace == ens%Capacity + 1) Then
                ! We need to create a space at the end of the tree
                i = IRnd(ens%Capacity)
                Call NoteContraction(ens, i, .False.)
                ens%Particles(i) = ens%Particles(ens%Capacity)
                ens%FirstSpace = ens%Capacity
            End If

            ! This ensures the next loop can go one entry
            ! past the particles just added.
            ens%Particles(ens%FirstSpace)%Properties = 0

            ! Set the leaves of the tree
            Do i=1, ens%FirstSpace - 1, 2
                ens%Tree%Nodes(i/2 + ens%HalfCapacity)%Left  = GetPreCalcs(ens%Particles(i))
                ens%Tree%Nodes(i/2 + ens%HalfCapacity)%Right = GetPreCalcs(ens%Particles(i+1))
            End Do

            ! Now calculate the remainder of the tree
            Do i = ens%HalfCapacity - 1, 1, -1
                ens%Tree%Nodes(i)%Left  = ens%Tree%Nodes(2*i)%Left    + ens%Tree%Nodes(2*i)%Right
                ens%Tree%Nodes(i)%Right = ens%Tree%Nodes(2*i +1)%Left + ens%Tree%Nodes(2*i +1)%Right
            End Do

            flag = N
            ens%Scaling = m

        Else
            ! Problem with arguments.
            flag = -2
        End If
    End Subroutine

    ! -------------------------------------------------------

    Subroutine SetTypedList(N, p, ens, flag)
        ! DESCRIPTION:
        !   Sets the ensemble particle list given a
        !   separate list of particles.

        Use SWPRNG
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN)             :: N    ! Number of particles to be added.
        Type(StochParticle), Intent(IN) :: p(N) ! Particles to be added to the ensemble.
        Type(Ensemble), Intent(INOUT)   :: ens  ! Ensemble to get.
        Integer, Intent(OUT)            :: flag ! Error flag.

        ! VARIABLES.
        Integer :: i, m, err
	Logical :: save_dble

        ! EXECUTABLE CODE.
        flag = 0

        ! Empty the tree, but keep the ActiveDoubling flag
        ! so that doubling can continue.  This shows my
        ! data flow and control is a bit mixed up since we
        ! treat ActiveDoubling differently to the
        ! sample volume / scaling.
        save_dble  = ens%ActiveDoubling
        Call ClearEnsemble(ens)
        ens%ActiveDoubling = save_dble

        If (N > ens%Capacity - 1) Then
           ! Tree is not big enough.
           flag = -1
        ElseIf(N <= Size(p) .And. N >= 0) Then
           ! N is a reasonable value.

           ! Number of times we can fit the
           ! list into the tree.
           m = 1
           If (ens%ActiveDoubling) Then
               m = Max((ens%Capacity - ens%DoublingSlack) / N, 1)
           End If
           If (N > ens%DoublingCutOff) Then
               ens%ActiveDoubling = USE_DOUBLING
           End If

           Do i = 1, N
               ! Put the particles into the list of tree contents.
               ens%Particles((i-1)*m+1:i*m) = p(i)
           End Do

           ens%FirstSpace = N * m + 1

           If (ens%FirstSpace == ens%Capacity + 1) Then
               ! We need to create a space at the end of the tree.
               i = IRnd(ens%Capacity)
               Call NoteContraction(ens, i, .False.)
               ens%Particles(i) = ens%Particles(ens%Capacity)
               ens%FirstSpace = ens%Capacity
           End If

           ! This ensures the next loop can go one entry
           ! past the particles just added to particle_list.
           ens%Particles(ens%FirstSpace)%Properties = 0

           ! Set the leaves of the tree.
           Do i=1, ens%FirstSpace - 1, 2
               ens%Tree%Nodes(i/2 + ens%HalfCapacity)%Left  = GetPreCalcs(ens%Particles(i))
               ens%Tree%Nodes(i/2 + ens%HalfCapacity)%Right = GetPreCalcs(ens%Particles(i+1))
           End Do

           ! Now calculate the remainder of the tree.
           Do i = ens%HalfCapacity - 1, 1, -1
               ens%Tree%Nodes(i)%Left  = ens%Tree%Nodes(2*i)%Left    + ens%Tree%Nodes(2*i)%Right
               ens%Tree%Nodes(i)%Right = ens%Tree%Nodes(2*i +1)%Left + ens%Tree%Nodes(2*i +1)%Right
           End Do

           flag = N
           ens%Scaling = m
        Else
           ! Problem with arguments.
           flag = -2
        End If
    End Subroutine

    ! -------------------------------------------------------

    Integer Function TreeIndex(i, ens)
        ! DESCRIPTION:
        !   Gets the index of the tree leaf referring to
        !   a particle.  The particle is given by its index
        !   in the particle list.
        !
        !   This function takes advantage of integer
        !   division.
        !
        ! RETURNS:
        !   Tree leaf index.
        Implicit None
        Integer, Intent(IN)        :: i
        Type(Ensemble), Intent(IN) :: ens
        TreeIndex = ens%HalfCapacity + ((i - 1)/2)
    End Function

    ! -------------------------------------------------------

    Subroutine GetScalingCounters(ens, D, E, C)
        ! DESCRIPTION:
        !   Returns the scaling counters for the binary tree.
        !   The counters are D=Doubling, E=Expansions,
        !   C=Contractions.
        Implicit None
        Type(Ensemble), Intent(IN) :: ens
        Integer, Intent(OUT)       :: D, E, C
        D = ens%NDBLE
        E = ens%NEXP
        C = ens%NCONT
    End Subroutine

    ! -------------------------------------------------------

    Subroutine GetCompositionMatrix(ens, comp)
        ! DESCRIPTION:
        !   Returns a matrix of particle compositions, that is
        !   the number of each type of species that have gone
        !   into making each particle.  Note:  This matrix
        !   is very large!

        Implicit None

        ! ARGUMENTS.

        ! The index of the particle for which to get
        ! the tree leaf index.
        Type(Ensemble), Intent(IN) :: ens
        Integer, Intent(OUT)       :: comp(:,:)

        ! VARIABLES.
        Integer :: i, N

        ! EXECUTABLE CODE.

        ! Initialise compostion matrix.
        comp = 0
        N    = ens%FirstSpace - 1

        ! Get the compositions of all known particles.
        If (Size(comp,1) >= MAX_COMP+1) Then
            Do i = 1, Min(Size(comp,2), N)
                comp(1,i)  = ens%Particles(i)%TypeID
                comp(2:,i) = ens%Particles(i)%Components
            End Do
        End If
    End Subroutine

    ! -------------------------------------------------------

    Subroutine GetAverageComposition(ens, comp, typeID)
        ! DESCRIPTION:
        !   Returns the average particle composition, that is
        !   the average number of each species which has gone
        !   unto each particle.

        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens       ! Ensemble to get from.
        Integer, Intent(OUT) :: comp(:)         ! Output composition array.
        Integer, Intent(IN), Optional :: typeID ! Optional type of particle for which to 
                                                ! get statistics.  If omitted gets statistics
                                                ! for all particles.

        ! VARIABLES.
        Integer :: i, N, M

        ! EXECUTABLE CODE.

        ! Initialise compostion matrix.
        comp = 0
        M = 0
        N = ens%FirstSpace - 1

        If ((N > 0) .And. (Size(comp) >= MAX_COMP)) Then
            If (Present(typeID)) Then
                ! Get the compositions of all known particles of the
                ! given type.
                Do i = 1, N
                    If (ens%Particles(i)%TypeID == typeID) Then
                        M = M + 1
                        comp = comp + ens%Particles(i)%Components
                    End If
                End Do                                                                              
                If (M > 0) comp = comp / M
            Else
                ! Get the compositions of all known particles.
                Do i = 1, N
                    comp = comp + ens%Particles(i)%Components
                End Do                                                                              
                comp = comp / N
           End If
        End If
    End Subroutine

    ! -------------------------------------------------------

    Subroutine WriteEnsemble(ens, fnum)
        ! DESCRIPTION:
        !   Writes a solution to a binary file.
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(IN) :: ens  ! Ensemble to output.
        Integer, Intent(IN)        :: fnum ! UNIT number on which to output.

        ! EXECUTABLE CODE.
        Write(fnum) ens%Particles, ens%Tree%Capacity, ens%Tree%Levels, ens%Tree%Nodes, &
                    ens%Capacity, ens%Scaling, ens%contFactor, ens%expFactor, ens%NEXP, &
                    ens%NCONT, ens%NDBLE, ens%CONTWARN, ens%FirstSpace, ens%HalfCapacity, &
                    ens%ActiveDoubling, ens%DoublingCutOff, ens%DoublingLimit, &
                    ens%DoublingSlack
    End Subroutine

    ! -------------------------------------------------------

    Subroutine ReadEnsemble(ens, fnum)
        ! DESCRIPTION:
        !   Reads a solution from a binary file.
        Implicit None

        ! ARGUMENTS.
        Type(Ensemble), Intent(INOUT) :: ens  ! Ensemble to read.
        Integer, Intent(IN)           :: fnum ! UNIT number on which to output.

        ! EXECUTABLE CODE.
        Read(fnum) ens%Particles, ens%Tree%Capacity, ens%Tree%Levels, ens%Tree%Nodes, &
                   ens%Capacity, ens%Scaling, ens%contFactor, ens%expFactor, ens%NEXP, &
                   ens%NCONT, ens%NDBLE, ens%CONTWARN, ens%FirstSpace, ens%HalfCapacity, &
                   ens%ActiveDoubling, ens%DoublingCutOff, ens%DoublingLimit, &
                   ens%DoublingSlack
    End Subroutine

End Module

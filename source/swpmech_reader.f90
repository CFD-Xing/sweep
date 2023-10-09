! *****************************************************************************
!
! File:				swpmech_reader.f90
! Project:			Sweep 2
! Author(s):			Matthew Celnik (msc37)
!
! Copyright (C) 2006  Matthew S Celnik
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
!   A module for reading a soot mechanism from an XML file.  This module requires
!   that the xmlf90 library is included with the build.
!
! Functions:
!   ReadMechXML         - Reads a soot mechanism from an XML file.
!   -- (xmlf90 handler routines) ----------------------------------------------
!   BeginElement        - Called when the parser begins to read an XML element.
!   EndElement          - Called when the parser stops reading an XML element.
!   HandleData          - Called when the parser reads a chunk of data.
!   -- (Reaction completion) --------------------------------------------------
!   CompleteInception    - Completes an inception reaction with info that is not
!                          contained in the XML file.
!   CompleteReaction     - Completes a surface reaction with info that is not
!                          contained in the XML file.
!   CompleteCondensation - Completes a condensation reaction with info that is not
!                          contained in the XML file.
! *****************************************************************************

Module SWPMECH_READER
    Use SWPMECH_TYPES

    Implicit None
    Public
    Private :: BeginElement, EndElement, HandleData
    Private :: CompleteInception, CompleteReaction, CompleteCondensation

    Integer, Private, Parameter :: LSYM = 20  ! Length of a symbolic string.
    Integer, Private, Parameter :: LSTR = 100 ! Length of a descriptive string.

    ! Flags indicating what type of element is currently being read.
    Logical, Private :: parseMech, parseModel, parseCoag !, parseTransform
    Logical, Private :: parseParticle, parseEnergyUnit, parseTracker
    Logical, Private :: parseComponent, parseDensity, parseMolWt
    Logical, Private :: parseInception, parseReaction, parseCondensation
    Logical, Private :: parseFormula, parseReactant, parseProduct
    Logical, Private :: parseParticleIn, parseParticleOut
    Logical, Private :: parseA, parseN, parseE, parseDX, parseParticleTerm

    Integer, Private :: error ! Error flag.

    ! Temporary mechanism into which to read from file.
    Type(Mechanism), Private :: m_mech

    ! Counters.
    Integer, Private :: jIcn, jRxn, jGrp

    Character(LEN=LSTR), Private :: m_mechName

    ! Array of species names.
    Character(LEN=LSYM), Private, Allocatable :: m_spnames(:)

    ! Particle symbols.
    Integer, Parameter, Private :: MAXSYMS=20 ! Max. particle symbol count.

    ! Particle components.
    Character(LEN=LSYM), Allocatable :: m_compNames(:) ! Component names.

    ! Energy units.
    Integer, Private, Parameter :: ENERGY_CAL = 1 ! Enerumation of different energy units.
    Integer, Private, Parameter :: ENERGY_SI  = 2
    Integer, Private :: m_units

    ! Inception and condensation reading.
    Real, Private :: m_m(2), m_d(2)

    Contains

    ! -------------------------------------------------------

    Subroutine ReadMechXML(file, mech, speciesNames, flag)
        ! DESCRIPTION.
        !   Reads a soot mechanism from an XML file.

        Use flib_sax, only: open_xmlfile, close_xmlfile, xml_t, xml_parse
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: file            ! File name.
        Type(Mechanism), Intent(OUT) :: mech            ! Mechanism to return.
        Character(LEN=*), Intent(IN) :: speciesNames(:) ! Names of gas-phase chemistry species.
        Integer, Intent(OUT)         :: flag            ! Error flag.

        ! VARIABLES.
        Type(xml_t) :: fxml
        Integer     :: iostat, i, err

        ! EXECUTABLE CODE.
        error = 0

        ! Allocate temporary memory.
        flag = AllocateTempMemory(Size(speciesNames))
        If (flag /= 0) Return
        m_spnames(1:Size(speciesNames)) = speciesNames(1:Size(speciesNames))

        ! Set all element flags to false.
        parseMech         = .False.
        parseModel        = .False.
        parseParticle     = .False.
        parseEnergyUnit   = .False.
        parseTracker      = .False.
        parseComponent    = .False.
        parseDensity      = .False.
        parseMolWt        = .False.
        parseInception    = .False.
        parseReaction     = .False.
        parseCondensation = .False.
        parseCoag         = .False.
        parseFormula      = .False.
        parseReactant     = .False.
        parseProduct      = .False.
        parseParticleIn   = .False.
        parseParticleOut  = .False.
        parseA            = .False.
        parseN            = .False.
        parseE            = .False.
        parseDX           = .False.
        parseParticleTerm = .False.

        ! Set counters to zero.
        jIcn = 0
        jRxn = 0
        jGrp = 0

        ! Parse the XML file.  This will fill the data structures.
        Call open_xmlfile(file, fxml, iostat)
        If (iostat == 0) Then
            Call xml_parse(fxml, begin_element_handler=BeginElement, &
                                 end_element_handler=EndElement, &
                                 pcdata_chunk_handler=HandleData)
        End If

        If (error == 0) Then
            ! Save particle type information.
            mech%NParticleTypes = m_mech%NParticleTypes
            If (mech%NParticleTypes > 0) Then
                Allocate(mech%ParticleTypeNames(mech%NParticleTypes), STAT=err)
                mech%ParticleTypeNames = m_mech%ParticleTypeNames(1:mech%NParticleTypes)
            End If

            ! Save inceptions.
            If (jIcn > 0) Then
                Allocate (mech%Inceptions(jIcn), STAT=err)
                mech%Inceptions = m_mech%Inceptions(1:jIcn)
                mech%InceptionCount = jIcn
            End If

            ! Save surface reactions.
            If (jRxn > 0) Then
                Allocate(mech%Reactions(jRxn), STAT=err)
                mech%Reactions = m_mech%Reactions(1:jRxn)
                mech%ReactionCount = jRxn
            End If

            ! Save reaction groups.
            If (jGrp > 0) Then
                Allocate(mech%Groups(jGrp+3), STAT=err)
                Allocate(mech%GroupNames(jGrp+3), STAT=err)
                mech%Groups(1:jGrp) = m_mech%Groups(1:jGrp)
                mech%GroupNames(1:jGrp) = m_mech%GroupNames(1:jGrp)
                mech%GroupCount = jGrp + EXTRA_GROUP_COUNT
            End If

            ! Save deferred process information.
            If (jIcn + jRxn > 0) Then
                Allocate(mech%DeferMask(jIcn+jRxn+EXTRA_PROCESS_COUNT), STAT=err)
                mech%DeferMask(1:jIcn) = .False.
                mech%DeferMask(jIcn+1:) = m_mech%DeferMask
                mech%AnyDeferred = Any(mech%DeferMask)
            End If

            ! Save coagulation information.
            If (m_mech%CoagRuleCount > 0) Then
                Allocate(mech%CoagRules(m_mech%CoagRuleCount), STAT=err)
                mech%CoagRules = m_mech%CoagRules(1:m_mech%CoagRuleCount)
                mech%CoagRuleCount = m_mech%CoagRuleCount
            End If

            ! Save component and tracking information.
            mech%Components     = m_mech%Components
            mech%ComponentCount = m_mech%ComponentCount
            mech%Trackers       = m_mech%Trackers
            mech%TrackerCount   = m_mech%TrackerCount

            ! Save model information.
            mech%UseHACA       = m_mech%UseHACA
            mech%ActSurfModel  = m_mech%ActSurfModel
            mech%ParticleModel = m_mech%ParticleModel

            ! Save other information.
            mech%ProcessCount = jIcn + jRxn + EXTRA_PROCESS_COUNT

            Print *, "Sweep:\> Successfully read """ // Trim(m_mechName) // """ mechanism from " // Trim(file) // "."
        Else
            ! An error occured somewhere!
            Print *, "Sweep:\> Failed to read mechanism!"
            flag = error
        End If

        ! Clear temporary memory.
        Call close_xmlfile(fxml)
        Call DeallocateTempMemory()
    End Subroutine

    ! -------------------------------------------------------
    ! XMLF90 HANDLER ROUTINES.
    !
    !   These routines are called by xmlf90 when parsing the
    !   mechanism file.
    !
    ! -------------------------------------------------------

    Subroutine BeginElement(name, attributes)
        ! DESCRIPTION.
        !   Called when the parser begins reading an element.

        Use flib_sax, only: has_key, dictionary_t, get_value
        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN)   :: name
        Type(dictionary_t), Intent(IN) :: attributes

        ! VARIABLES.
        Integer             :: stat, i
        Character(LEN=LSTR) :: str
        Character(LEN=LSYM) :: id
        Logical             :: flag

        ! EXECUTABLE CODE.
        Select Case (name)
            Case ("mechanism")

                parseMech = .True.

                ! Mechanism name.
                If (has_key(attributes, "name")) Then
                    Call get_value(attributes, "name", m_mechName, stat)
                Else
                    m_mechName = "no name"
                End If

            Case ("model")

                parseModel = .True.

                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", id, stat)

                    Select Case (id)
                        Case ("haca")
                            ! Turn on HACA model (requires certain species
                            ! to be present in the chemistry).
                            m_mech%UseHACA = .True.
                        Case ("actsurf")
                            ! An active surface model.
                            If (has_key(attributes, "type")) Then
                                Call get_value(attributes, "type", str, stat)
                                
                                Select Case (str)
                                    Case ("abf")
                                        m_mech%ActSurfModel = ACT_SURF_ABF
                                End Select
                            Else
                                ! Active surface model should have a type!
                                error = -1
                            End If
                    End Select
                Else
                    ! Model should have an ID!
                    error = -1
                End If

            Case ("particle")

                parseParticle = .True.
                m_mech%NParticleTypes = Min(m_mech%NParticleTypes+1, MAXSYMS)

                ! Particle symbol.
                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", m_mech%ParticleTypeNames(m_mech%NParticleTypes), stat)
                End If

            Case ("component")

                parseComponent = .True.

                ! Get component ID.
                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", id, stat)
                Else
                    id = ""
                End If

                If (parseInception .Or. parseReaction .Or. parseCondensation) Then
                    ! Get component index.
                    i = IndexOf(id, m_compNames)

                    ! Get component change.
                    If (has_key(attributes, "dx")) Then
                        Call get_value(attributes, "dx", str, stat)
                    Else
                        str = ""
                    End If

                    ! Note component change.
                    If (parseInception) Then
                        m_mech%Inceptions(jIcn)%dComp(i) = CReal(str, 0.0E0)
                    ElseIf (parseReaction) Then
                        m_mech%Reactions(jRxn)%dComp(i) = CReal(str, 0.0E0)
                    ElseIf (parseCondensation) Then
                        m_mech%Reactions(jRxn)%dComp(i)   = CReal(str, 0.0E0)
                        m_mech%Reactions(jRxn-1)%dComp(i) = m_mech%Reactions(jRxn)%dComp(i)
                        m_mech%Reactions(jRxn-2)%dComp(i) = m_mech%Reactions(jRxn)%dComp(i)
                    End If
                Else
                    ! Add a new component symbol to mechanism.
                    m_mech%ComponentCount = Min(m_mech%ComponentCount+1, MAX_COMP)
                    m_mech%Components(m_mech%ComponentCount)%Symbol = AdjustL(id)
                    m_compNames(m_mech%ComponentCount) = AdjustL(id)
                End If

            Case ("density")

                parseDensity = .True.

            Case ("molwt")

                parseMolWt = .True.

            Case ("track")

                parseTracker = .True.

                ! Get tracker ID.
                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", id, stat)
                Else
                    id = ""
                End If

                If (parseInception .Or. parseReaction .Or. parseCondensation) Then
                    ! Get tracker index.
                    i = IndexOf(id, m_mech%Trackers)

                    ! Get tracker change.
                    If (has_key(attributes, "dx")) Then
                        Call get_value(attributes, "dx", str, stat)
                    Else
                        str = ""
                    End If

                    ! Note tracker change.
                    If (parseInception) Then
                        m_mech%Inceptions(jIcn)%dTrack(i) = CReal(str, 0.0E0)
                    ElseIf (parseReaction) Then
                        m_mech%Reactions(jRxn)%dTrack(i) = CReal(str, 0.0E0)
                    ElseIf (parseCondensation) Then
                        m_mech%Reactions(jRxn)%dTrack(i)   = CReal(str, 0.0E0)
                        m_mech%Reactions(jRxn-1)%dTrack(i) = m_mech%Reactions(jRxn)%dTrack(i)
                        m_mech%Reactions(jRxn-2)%dTrack(i) = m_mech%Reactions(jRxn)%dTrack(i)
                    End If
                Else
                    ! We are noting a new tracker variable.
                    m_mech%TrackerCount = Min(m_mech%TrackerCount+1, MAX_TRACK)
                    m_mech%Trackers(m_mech%TrackerCount) = AdjustL(id)
                End If

            Case ("inception")

                parseInception = .True.
                
                jIcn = jIcn + 1
                jGrp = jGrp + 1
                m_mech%Groups(jGrp) = 1

                ! Inception name.
                If (has_key(attributes, "name")) Then
                    Call get_value(attributes, "name", m_mech%GroupNames(jGrp), stat)
                End If

                m_m = 0.0
                m_d = 0.0

            Case ("reaction")

                ! Reaction name.
                jGrp = jGrp + 1
                If (has_key(attributes, "name")) Then
                    Call get_value(attributes, "name", m_mech%GroupNames(jGrp), stat)
                End If

                ! Deferred process, default is not deferred.
                If (has_key(attributes, "defer")) Then
                    Call get_value(attributes, "defer", str, stat)
                    If (str == "true") Then
                        flag = .True.
                    Else
                        flag = .False.
                    End If
                Else
                    flag = .False.
                End If

                ! Check reaction type, default "surface".
                If (has_key(attributes, "type")) Then
                    Call get_value(attributes, "type", str, stat)
                Else 
                    str = "surface"
                End If

                ! Increment counters based on reaction type.
                If (str == "surface") Then
                    jRxn = jRxn + 1
                    m_mech%Groups(jGrp) = 1
                    m_mech%DeferMask(jRxn) = flag
                    parseReaction = .True.
                ElseIf (str == "condensation") Then
                    jRxn = jRxn + 3
                    m_mech%Groups(jGrp) = 3
                    m_mech%DeferMask(jRxn-2:jRxn) = flag
                    parseCondensation = .True.
                    m_m = 0.0
                    m_d = 0.0
                ElseIf (str == "transform") Then
                    jRxn = jRxn + 1
                    m_mech%Groups(jGrp) = 1
                    m_mech%DeferMask(jRxn) = .False.
                    parseReaction = .True.
                End If

            Case ("coagulation")
            
                parseCoag = .True.
                m_mech%CoagRuleCount = m_mech%CoagRuleCount + 1

            Case ("formula")

                parseFormula = .True.

            Case ("reactant")

                parseReactant = .True.

                If (parseInception) Then
                    ! Add another reactant to this inception.
                    m_mech%Inceptions(jIcn)%NREAC = m_mech%Inceptions(jIcn)%NREAC + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        i = IndexOf(id, m_spnames)
                        m_mech%Inceptions(jIcn)%Species(m_mech%Inceptions(jIcn)%NREAC) = IndexOf(id, m_spnames)
                    Else
                        m_mech%Inceptions(jIcn)%Species(m_mech%Inceptions(jIcn)%NREAC) = 0
                    End If
                     
                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Inceptions(jIcn)%Stoich(m_mech%Inceptions(jIcn)%NREAC) = CReal(str, 1.0E0)
                    Else
                        m_mech%Inceptions(jIcn)%Stoich(m_mech%Inceptions(jIcn)%NREAC) = 1.0E0
                    End If

                    ! Read species size data.
                    If (has_key(attributes, "m")) Then
                        Call get_value(attributes, "m", str, stat)
                        If (m_m(1) == 0.0) Then
                            m_m(1) = CReal(str, 0.0E0)
                            m_m(2) = m_m(1) ! Assume dimer.
                        Else
                            m_m(2) = CReal(str, 0.0E0)
                        End If
                    End If
                    If (has_key(attributes, "d")) Then
                        Call get_value(attributes, "d", str, stat)
                        If (m_d(1) == 0.0) Then
                            m_d(1) = CReal(str, 0.0E0)
                            m_d(2) = m_d(1) ! Assume dimer.
                        Else
                            m_d(2) = CReal(str, 0.0E0)
                        End If
                    End If
                ElseIf (parseReaction) Then
                    ! Add another reactant to this reaction.
                    m_mech%Reactions(jRxn)%NREAC = m_mech%Reactions(jRxn)%NREAC + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC) = IndexOf(id, m_spnames)
                    Else
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC) = 0
                    End If

                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC) = CReal(str, 1.0E0)
                    Else
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC) = 1.0E0
                    End If
                ElseIf (parseCondensation) Then
                    ! Add another reactant to this condensation.
                    m_mech%Reactions(jRxn-2)%NREAC = m_mech%Reactions(jRxn-2)%NREAC + 1
                    m_mech%Reactions(jRxn-1)%NREAC = m_mech%Reactions(jRxn-1)%NREAC + 1
                    m_mech%Reactions(jRxn)%NREAC   = m_mech%Reactions(jRxn)%NREAC + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC)     = IndexOf(id, m_spnames)
                        m_mech%Reactions(jRxn-1)%Species(m_mech%Reactions(jRxn-1)%NREAC) = m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC)
                        m_mech%Reactions(jRxn-2)%Species(m_mech%Reactions(jRxn-2)%NREAC) = m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC)
                    Else
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC)     = 0
                        m_mech%Reactions(jRxn-1)%Species(m_mech%Reactions(jRxn-1)%NREAC) = 0
                        m_mech%Reactions(jRxn-2)%Species(m_mech%Reactions(jRxn-2)%NREAC) = 0
                    End If

                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC)     = CReal(str, 1.0E0)
                        m_mech%Reactions(jRxn-1)%Stoich(m_mech%Reactions(jRxn-1)%NREAC) = m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC)
                        m_mech%Reactions(jRxn-2)%Stoich(m_mech%Reactions(jRxn-2)%NREAC) = m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC)
                    Else
                        m_mech%Reactions(jRxn-2)%Stoich(m_mech%Reactions(jRxn-2)%NREAC) = 1.0E0
                        m_mech%Reactions(jRxn-1)%Stoich(m_mech%Reactions(jRxn-1)%NREAC) = 1.0E0
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC)     = 1.0E0
                    End If

                    ! Read species size data.
                    If (has_key(attributes, "m")) Then
                        Call get_value(attributes, "m", str, stat)
                        m_m(1) = CReal(str, 0.0E0)
                    End If
                    If (has_key(attributes, "d")) Then
                        Call get_value(attributes, "d", str, stat)
                        m_d(1) = CReal(str, 0.0E0)
                    End If
                Else
                    ! A misplaced element!
                    error = -1
                End If
                 
            Case ("product")

                parseProduct = .True.

                If (parseInception) Then
                    ! Add another product to this inception.
                    m_mech%Inceptions(jIcn)%NPROD = m_mech%Inceptions(jIcn)%NPROD + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        m_mech%Inceptions(jIcn)%Species(m_mech%Inceptions(jIcn)%NREAC+ &
                                                        m_mech%Inceptions(jIcn)%NPROD) = IndexOf(id, m_spnames)
                    Else
                        m_mech%Inceptions(jIcn)%Species(m_mech%Inceptions(jIcn)%NREAC+ &
                                                        m_mech%Inceptions(jIcn)%NPROD) = 0
                    End If

                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Inceptions(jIcn)%Stoich(m_mech%Inceptions(jIcn)%NREAC+ &
                                                       m_mech%Inceptions(jIcn)%NPROD) = CReal(str, 1.0E0)
                    Else
                        m_mech%Inceptions(jIcn)%Stoich(m_mech%Inceptions(jIcn)%NREAC+ &
                                                       m_mech%Inceptions(jIcn)%NPROD) = 1.0E0
                    End If
                ElseIf (parseReaction) Then
                    ! Add another product to this reaction.
                    m_mech%Reactions(jRxn)%NPROD = m_mech%Reactions(jRxn)%NPROD + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC+ &
                                                       m_mech%Reactions(jRxn)%NPROD) = IndexOf(id, m_spnames)
                    Else
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC+ &
                                                       m_mech%Reactions(jRxn)%NPROD) = 0
                    End If

                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC+ &
                                                      m_mech%Reactions(jRxn)%NPROD) = CReal(str, 1.0E0)
                    Else
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC+ &
                                                      m_mech%Reactions(jRxn)%NPROD) = 1.0E0
                    End If
                ElseIf (parseCondensation) Then
                    ! Add another product to this condensation.
                    m_mech%Reactions(jRxn-2)%NPROD = m_mech%Reactions(jRxn-2)%NPROD + 1
                    m_mech%Reactions(jRxn-1)%NPROD = m_mech%Reactions(jRxn-1)%NPROD + 1
                    m_mech%Reactions(jRxn)%NPROD   = m_mech%Reactions(jRxn)%NPROD + 1

                    ! Read species ID.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC+ &
                                                       m_mech%Reactions(jRxn)%NPROD) = IndexOf(id, m_spnames)
                        m_mech%Reactions(jRxn-1)%Species(m_mech%Reactions(jRxn-1)%NREAC + &
                                                         m_mech%Reactions(jRxn-1)%NPROD) &
                                                       = m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC + &
                                                                                        m_mech%Reactions(jRxn)%NPROD)
                        m_mech%Reactions(jRxn-2)%Species(m_mech%Reactions(jRxn-2)%NREAC + &
                                                         m_mech%Reactions(jRxn-2)%NPROD) &
                                                       = m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC+ &
                                                                                        m_mech%Reactions(jRxn)%NPROD)
                    Else
                        m_mech%Reactions(jRxn)%Species(m_mech%Reactions(jRxn)%NREAC)     = 0
                        m_mech%Reactions(jRxn-1)%Species(m_mech%Reactions(jRxn-1)%NREAC) = 0
                        m_mech%Reactions(jRxn-2)%Species(m_mech%Reactions(jRxn-2)%NREAC) = 0
                    End If

                    ! Read the stoichiometry, default 1.0E0.
                    If (has_key(attributes, "stoich")) Then
                        Call get_value(attributes, "stoich", str, stat)
                        m_mech%Reactions(jRxn-2)%Stoich(m_mech%Reactions(jRxn-2)%NREAC+ &
                                                        m_mech%Reactions(jRxn-2)%NPROD) = CReal(str, 1.0E0)
                        m_mech%Reactions(jRxn-1)%Stoich(m_mech%Reactions(jRxn-1)%NREAC+ &
                                                        m_mech%Reactions(jRxn-1)%NPROD) = CReal(str, 1.0E0)
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC+ &
                                                      m_mech%Reactions(jRxn)%NPROD) = CReal(str, 1.0E0)
                    Else
                        m_mech%Reactions(jRxn-2)%Stoich(m_mech%Reactions(jRxn-2)%NREAC+ &
                                                        m_mech%Reactions(jRxn-2)%NPROD) = 1.0E0
                        m_mech%Reactions(jRxn-1)%Stoich(m_mech%Reactions(jRxn-1)%NREAC+ &
                                                        m_mech%Reactions(jRxn-1)%NPROD) = 1.0E0
                        m_mech%Reactions(jRxn)%Stoich(m_mech%Reactions(jRxn)%NREAC+ &
                                                      m_mech%Reactions(jRxn)%NPROD) = 1.0E0
                    End If
                Else
                    ! A misplaced element!
                    error = -1
                End If
 
            Case ("particlein")

                parseParticleIn = .True.

                ! Get particle ID.
                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", id, stat)
                    i = IndexOf(id, m_mech%ParticleTypeNames)
                Else
                    i = 1
                End If

                If (parseReaction) Then
                    m_mech%Reactions(jRxn)%PTypeIn = Max(i-1,0)
                ElseIf (parseCondensation) Then
                    m_mech%Reactions(jRxn-2)%PTypeIn = Max(i-1,0)
                    m_mech%Reactions(jRxn-1)%PTypeIn = Max(i-1,0)
                    m_mech%Reactions(jRxn)%PTypeIn   = Max(i-1,0)
                ElseIf (parseCoag) Then
                    If (m_mech%CoagRules(m_mech%CoagRuleCount)%In1 < 0) Then
                        m_mech%CoagRules(m_mech%CoagRuleCount)%In1 = Max(i-1,0)
                    Else
                        m_mech%CoagRules(m_mech%CoagRuleCount)%In2 = Max(i-1,0)
                    End If
                End If            

            Case ("particleout")

                parseParticleOut = .True.

                ! Get particle ID.
                If (has_key(attributes, "id")) Then
                    Call get_value(attributes, "id", id, stat)
                    i = IndexOf(id, m_mech%ParticleTypeNames)
                Else
                    i = 1
                End If

                If (parseReaction) Then
                    m_mech%Reactions(jRxn)%PTypeOut = Max(i-1,0)
                ElseIf (parseCondensation) Then
                    m_mech%Reactions(jRxn-2)%PTypeOut = Max(i-1,0)
                    m_mech%Reactions(jRxn-1)%PTypeOut = Max(i-1,0)
                    m_mech%Reactions(jRxn)%PTypeOut   = Max(i-1,0)
                ElseIf (parseCoag) Then
                    m_mech%CoagRules(m_mech%CoagRuleCount)%Out = Max(i-1,0)
                End If            

            Case ("A")

                parseA = .True.

            Case ("n")

                parseN = .True.

            Case ("E")

                parseE = .True.
                
                ! Read units.
                If (has_key(attributes, "units")) Then
                    Call get_value(attributes, "units", str, stat)
                    Select Case (str)
                        Case ("cal","CAL","Calories","CALORIES","calories")
                            m_units = ENERGY_CAL
                        Case ("SI", "si", "Si")
                            m_units = ENERGY_SI
                        Case Default
                            m_units = ENERGY_SI
                    End Select
                Else
                    m_units = ENERGY_SI
                End If

            Case ("particleterm")

                parseParticleTerm = .True.

                If (parseReaction) Then
                    ! Get ID of soot term, default is "v"=volume.
                    If (has_key(attributes, "id")) Then
                        Call get_value(attributes, "id", id, stat)
                    Else
                        id = "v"
                    End If

                    ! Note soot term powers, default 1.0E0.
                    If (has_key(attributes, "power")) Then
                        Call get_value(attributes, "power", str, stat)
                    Else
                        str = ""
                    End If

                    Select Case (id)
                        Case ("v")
                            m_mech%Reactions(jRxn)%P(1) = CReal(str, 1.0E0)
                        Case ("s")
                            m_mech%Reactions(jRxn)%P(2) = CReal(str, 1.0E0)
                        Case ("as")
                            m_mech%Reactions(jRxn)%P(3) = CReal(str, 1.0E0)
                        Case ("d")
                            m_mech%Reactions(jRxn)%P(4) = CReal(str, 1.0E0)
                        Case Default
                            ! An unknown particle term!
                            error = -1
                    End Select
                Else
                    ! A misplaced element!
                    error = -1
                End If

        End Select
    End Subroutine

    ! -------------------------------------------------------

    Subroutine EndElement(name)
        ! DESCRIPTION.
        !   Called when the parser stops reading an element.

        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: name

        ! VARIABLES.

        ! EXECUTABLE CODE.
        Select Case (name)
            Case ("mechanism")
                parseMech = .False.
            Case ("particle")
                parseParticle = .False.
            Case ("component")
                parseComponent = .False.
            Case ("density")
                parseDensity = .False.
            Case ("molwt")
                parseMolWt = .False.
            Case ("energyunits")
                parseEnergyUnit = .False.
            Case ("track")
                parseTracker = .False.
            Case ("inception")
                Call CompleteInception(jIcn)
                parseInception = .False.
            Case ("reaction")
                If (parseReaction) Then
                    Call CompleteReaction(jRxn)
                    parseReaction = .False.
                ElseIf (parseCondensation) Then
                    Call CompleteCondensation(jRxn-2)
                    parseCondensation = .False.
                End If
            Case ("coagulation")
                parseCoag = .False.
            Case ("formula")
                parseFormula = .False.
            Case ("reactant")
                parseReactant = .False.
            Case ("product")
                parseProduct = .False.
            Case ("particlein")
                parseParticleIn = .False.
            Case ("particleout")
                parseParticleOut = .False.
            Case ("A")
                parseA = .False.
            Case ("n")
                parseN = .False.
            Case ("E")
                parseE = .False.
            Case ("particleterm")
                parseParticleTerm = .False.
        End Select
    End Subroutine

    ! -------------------------------------------------------

    Subroutine HandleData(chunk)
        ! DESCRIPTION.
        !   Called when the parser handles a chunk of data.

        Use StrConv
        Implicit None

        ! ARGUMENTS.
        Character(LEN=*), Intent(IN) :: chunk

        ! VARIABLES.
        Integer :: i

        ! EXECUTABLE CODE.

        If (parseDensity) Then
            
            m_mech%Components(m_mech%ComponentCount)%Density = CReal(chunk, 0.0E0)

        ElseIf (parseMolWt) Then

             m_mech%Components(m_mech%ComponentCount)%MolWt = CReal(chunk, 0.0E0)

        ElseIf (parseA) Then

            If (parseReaction) Then
                m_mech%Reactions(jRxn)%A = CReal(chunk, 0.0E0)
            ElseIf (parseCondensation) Then
                m_mech%Reactions(jRxn)%A = CReal(chunk, 1.0E0)
                m_mech%Reactions(jRxn-1)%A = m_mech%Reactions(jRxn)%A
                m_mech%Reactions(jRxn-2)%A = m_mech%Reactions(jRxn)%A
            ElseIf (parseInception) Then
                m_mech%Inceptions(jIcn)%A = CReal(chunk, 1.0E0)
            End If

        ElseIf (parseN) Then

            If (parseReaction) Then
                m_mech%Reactions(jRxn)%n = CReal(chunk, 0.0E0)
            End If

        ElseIf (parseE) Then

            If (parseReaction) Then
                m_mech%Reactions(jRxn)%E = CReal(chunk, 0.0E0)

                If (m_units == ENERGY_CAL) Then
                    m_mech%Reactions(jRxn)%E = m_mech%Reactions(jRxn)%E * R / RCAL
                End If
            End If

        End If
    End Subroutine

    ! -------------------------------------------------------
    ! PROCESS COMPLETION ROUTINES.
    !
    !   These routines are called by when an inception or
    !   reaction element is completely read from the XML.
    !   They are used to complete any missing parts from the
    !   reaction types.
    !
    ! -------------------------------------------------------

    Subroutine CompleteInception(i)
        ! DESCRIPTION.
        !   Completes an inception reaction with all information
        !   the is not explicitly in the XML file.

        Use SWPCOAG_MODEL
        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN) :: i ! Index of inception reaction.

        ! EXECUTABLE CODE.
        Call NucSlipFlowK(m_d(1), m_d(2), m_mech%Inceptions(i)%ksf1, m_mech%Inceptions(i)%ksf2)
        m_mech%Inceptions(i)%kfm = NucFreeMolK(m_d(1), m_d(2), m_m(1), m_m(2))
        m_mech%Inceptions(i)%A   = m_mech%Inceptions(i)%A * 0.5E0
        m_mech%Inceptions(i)%dS  = EquivSphereSurface(m_mech%Inceptions(i)%dComp, &
                                                      m_mech%Components, m_mech%ComponentCount)
    End Subroutine
 
    ! -------------------------------------------------------

    Subroutine CompleteReaction(i)
        ! DESCRIPTION.
        !   Completes a surface reaction with all information
        !   that is not explicitly in the XML file.

        Use SWPPART
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN) :: i ! Index of surface reaction.

        ! EXECUTABLE CODE.

        ! Note correct soot property dependence.
        If (All(m_mech%Reactions(i)%P == 0.0E0)) Then
            m_mech%Reactions(i)%ID = UNIFORM_ID
        Else
            m_mech%Reactions(i)%ID = GetPropertyIndex(m_mech%Reactions(i)%P)
        End If

        ! If active surface dependence then include number of
        ! surface sites in rate constant.
        If (m_mech%Reactions(i)%P(3) == 1.0E0) m_mech%Reactions(i)%A = m_mech%Reactions(i)%A * CHI

        ! Note how particle surface changes.
        If (Sum(m_mech%Reactions(i)%dComp) > 0.0E0) Then
            m_mech%Reactions(i)%dSID = GROWTH_RADIUS
        Else
            m_mech%Reactions(i)%dSID = OXIDATION_RADIUS
        End If

        ! Cannot defer a reaction which changes a particle type.
        m_mech%DeferMask(i) = m_mech%DeferMask(i) .And. &
                             (m_mech%Reactions(i)%PTypeIn == m_mech%Reactions(i)%PTypeOut)
    End Subroutine

    ! -------------------------------------------------------

    Subroutine CompleteCondensation(i)
        ! DESCRIPTION.
        !   Completes a condensation reaction with all information
        !   that is not explicitly in the XML file.

        Use SWPCOAG
        USe SWPPART
        Implicit None

        ! ARGUMENTS.
        Integer, Intent(IN) :: i ! Index of first (of 3) condensation reaction.

        ! VARIABLES.
	Real :: k1, k2, k3, d, m

        ! EXECUTABLE CODE.
        
        ! Get free-molecular condensation terms.
        Call FreeMolCondTerms(m_d(1), m_m(1), k1, k2, k3)

        ! Reaction 1.
        m_mech%Reactions(i)%A  = m_mech%Reactions(i)%A * NA * k1
        m_mech%Reactions(i)%n  = 0.5E0
        m_mech%Reactions(i)%E  = 0.0E0
        m_mech%Reactions(i)%ID = UNIFORM_ID

        ! Reaction 2.
        m_mech%Reactions(i+1)%A  = m_mech%Reactions(i+1)%A * NA * k2
        m_mech%Reactions(i+1)%n  = 0.5E0
        m_mech%Reactions(i+1)%E  = 0.0E0
        m_mech%Reactions(i+1)%ID = iD

        ! Reaction 3.
        m_mech%Reactions(i+2)%A  = m_mech%Reactions(i+2)%A * NA * k3
        m_mech%Reactions(i+2)%n  = 0.5E0
        m_mech%Reactions(i+2)%E  = 0.0E0
        m_mech%Reactions(i+2)%ID = iD2

        ! All reactions.
        If (Sum(m_mech%Reactions(i)%dComp) > 0.0E0) Then
            m_mech%Reactions(i)%dSID   = GROWTH_RADIUS
            m_mech%Reactions(i+1)%dSID = GROWTH_RADIUS
            m_mech%Reactions(i+2)%dSID = GROWTH_RADIUS
        Else
            m_mech%Reactions(i)%dSID   = OXIDATION_RADIUS
            m_mech%Reactions(i+1)%dSID = OXIDATION_RADIUS
            m_mech%Reactions(i+2)%dSID = OXIDATION_RADIUS
        End If
    End Subroutine
       
    ! -------------------------------------------------------

    Integer Function AllocateTempMemory(nsp)
        ! DESCRIPTION.
        !   Allocates temporary memory required to read a
        !   mechanism from a XML file.
        Implicit None
        Integer, Intent(IN) :: nsp
        Integer :: err

        Allocate(m_mech%Inceptions(MAXICNS), m_mech%Reactions(MAXRXNS), &
                 m_mech%Groups(MAXGRPS), m_mech%GroupNames(MAXGRPS), &
                 m_mech%DeferMask(MAXRXNS), m_mech%ParticleTypeNames(MAXSYMS), &
                 m_compNames(MAX_COMP), m_spnames(nsp), &
                 m_mech%CoagRules(MAXSYMS*(MAXSYMS+1)/2), STAT=err)

        If (err /= 0) Then
            AllocateTempMemory = -1
        Else
            AllocateTempMemory = 0
            m_mech%Groups = 0
            m_mech%GroupNames = ""
            m_mech%DeferMask = .False.
            m_mech%Trackers = ""
            m_mech%NParticleTypes = 0
            m_mech%ParticleTypeNames = ""
            m_compNames = ""
            m_spnames = ""
        End If
    End Function

    ! -------------------------------------------------------

    Subroutine DeallocateTempMemory()
        ! DESCRIPTION.
        !   Deallocates temporary memory used to read a
        !   mechanism from a XML file.
        Implicit None
        Integer :: err
        Deallocate(m_mech%Inceptions, m_mech%Reactions, &
                   m_mech%Groups, m_mech%GroupNames, &
                   m_mech%DeferMask, m_mech%ParticleTypeNames, &
                   m_compNames, m_mech%CoagRules, STAT=err)
    End Subroutine

End Module

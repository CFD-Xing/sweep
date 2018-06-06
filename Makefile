
# Compiler options.

f77=ifort
f90=ifort
flag=-c -O2 -fpp -r8
flag77=-c -80 -fpp -r8

# Source directories.

srcdir=source
auxdir=auxilliary
xmldir=xmlf90

###########################################################################
# Set directory and file names:
###########################################################################

all:    profile_driver_program.o profile_driver_input.o profile_driver_output.o \
	strconv.o contable.o settings_file_reader.o list_t.o \
	swpensemble.o swpcoag.o swpcoag_model.o swpstats.o swperr.o swpchem.o swpmech.o swpmech_types.o mt19937.o \
	swppart.o swpchem_shared.o swpparams.o swpstep.o sweep.o swprng.o swpprocess.o swpsoln.o swpmech_reader.o \
	flib_sax.o m_buffer.o m_charset.o m_converters.o m_debug.o m_dictionary.o m_elstack.o m_entities.o m_fsm.o \
	m_io.o m_reader.o m_xml_error.o m_xml_parser.o

	$(f90) -r8 -o bin_sect/sweep_aggregation_rho_trung_dc.x profile_driver_program.o profile_driver_input.o profile_driver_output.o \
	strconv.o contable.o settings_file_reader.o list_t.o \
	swpensemble.o swpcoag.o swpcoag_model.o swpstats.o swperr.o swpchem.o swpmech.o swpmech_types.o mt19937.o \
	swppart.o swpchem_shared.o swpparams.o swpstep.o sweep.o swprng.o swpprocess.o swpsoln.o swpmech_reader.o \
	flib_sax.o m_buffer.o m_charset.o m_converters.o m_debug.o m_dictionary.o m_elstack.o m_entities.o m_fsm.o \
	m_io.o m_reader.o m_xml_error.o m_xml_parser.o

# Main Program Source.

profile_driver_program.o: $(srcdir)/profile_driver_program.f90 profile_driver_input.o profile_driver_output.o sweep.o
	$(f90) $(flag) $(srcdir)/profile_driver_program.f90

profile_driver_input.o: $(srcdir)/profile_driver_input.f90 strconv.o settings_file_reader.o sweep.o list_t.o
	$(f90) $(flag) $(srcdir)/profile_driver_input.f90

profile_driver_output.o: $(srcdir)/profile_driver_output.f90 contable.o strconv.o sweep.o swpchem.o swpstats.o
	$(f90) $(flag) $(srcdir)/profile_driver_output.f90

# Auxilliary Source.

strconv.o:    $(auxdir)/strconv.f90
	$(f90) $(flag) $(auxdir)/strconv.f90

contable.o:   $(auxdir)/contable.f90
	$(f90) $(flag) $(auxdir)/contable.f90

settings_file_reader.o:    $(auxdir)/settings_file_reader.f90 strconv.o
	$(f90) $(flag) $(auxdir)/settings_file_reader.f90

list_t.o: $(auxdir)/list_t.f90
	$(f90) $(flag) $(auxdir)/list_t.f90

# Sweep 2 Source.

swprng.o:  $(srcdir)/swprng.f90 mt19937.o
	$(f90) $(flag) $(srcdir)/swprng.f90

swpensemble.o:  $(srcdir)/swpensemble.f90  swpparams.o swppart.o swprng.o
	$(f90) $(flag) $(srcdir)/swpensemble.f90

swpcoag.o:   $(srcdir)/swpcoag.f90 swpchem_shared.o swpparams.o swppart.o swpcoag_model.o
	$(f90) $(flag) $(srcdir)/swpcoag.f90

swpcoag_model.o:   $(srcdir)/swpcoag_model.f90 swpchem_shared.o swpparams.o swppart.o
	$(f90) $(flag) $(srcdir)/swpcoag_model.f90

swpstats.o:   $(srcdir)/swpstats.f90 swpensemble.o swpparams.o swppart.o swpsoln.o
	$(f90) $(flag) $(srcdir)/swpstats.f90

swperr.o:   $(srcdir)/swperr.f90
	$(f90) $(flag) $(srcdir)/swperr.f90

swpchem.o:   $(srcdir)/swpchem_profile.f90 swpchem_shared.o swperr.o swpparams.o
	$(f90) $(flag) -o swpchem.o $(srcdir)/swpchem_profile.f90

swpprocess.o:   $(srcdir)/swpprocess.f90 swpcoag.o swpcoag_model.o swpparams.o swpchem.o swpensemble.o \
	swperr.o swppart.o swprng.o swpmech.o swpsoln.o
	$(f90) $(flag) $(srcdir)/swpprocess.f90

mt19937.o:   $(srcdir)/mt19937.f
	$(f77) $(flag77) $(srcdir)/mt19937.f

swppart.o:   $(srcdir)/swppart.f90 swpparams.o swpmech_types.o
	$(f90) $(flag) $(srcdir)/swppart.f90

swpchem_shared.o:   $(srcdir)/swpchem_shared.f90 strconv.o swperr.o swpparams.o swpmech_types.o
	$(f90) $(flag) $(srcdir)/swpchem_shared.f90

swpparams.o:   $(srcdir)/swpparams.f90
	$(f90) $(flag) $(srcdir)/swpparams.f90

swpstep.o:   $(srcdir)/swpstep.f90 swpchem.o swpensemble.o swperr.o swpmech.o swpparams.o swpprocess.o \
	swprng.o swpsoln.o
	$(f90) $(flag) $(srcdir)/swpstep.f90

sweep.o:   $(srcdir)/sweep.f90 swpparams.o swpchem.o swpensemble.o swperr.o swpmech.o swppart.o swpstep.o \
	swpstats.o swprng.o swpprocess.o swpsoln.o
	$(f90) $(flag) $(srcdir)/sweep.f90

swpmech_reader.o:   $(srcdir)/swpmech_reader.f90 swpparams.o swpcoag.o swpcoag_model.o swppart.o \
	strconv.o flib_sax.o swpmech_types.o 
	$(f90) $(flag) $(srcdir)/swpmech_reader.f90

swpmech.o:   $(srcdir)/swpmech.f90 swperr.o swpmech_reader.o swpmech_types.o swppart.o
	$(f90) $(flag) $(srcdir)/swpmech.f90

swpmech_types.o:   $(srcdir)/swpmech_types.f90 swpparams.o
	$(f90) $(flag) $(srcdir)/swpmech_types.f90

swpsoln.o:   $(srcdir)/swpsoln.f90 swpchem.o swpensemble.o swpmech_types.o 
	$(f90) $(flag) $(srcdir)/swpsoln.f90

# XMLF90 source.

flib_sax.o:   $(xmldir)/flib_sax.f90 m_converters.o m_dictionary.o m_xml_error.o m_xml_parser.o
	$(f90) $(flag) $(xmldir)/flib_sax.f90

m_buffer.o:   $(xmldir)/m_buffer.f90
	$(f90) $(flag) $(xmldir)/m_buffer.f90

m_charset.o:   $(xmldir)/m_charset.f90
	$(f90) $(flag) $(xmldir)/m_charset.f90

m_converters.o:   $(xmldir)/m_converters.f90 m_debug.o
	$(f90) $(flag) $(xmldir)/m_converters.f90

m_debug.o:   $(xmldir)/m_debug.f90
	$(f90) $(flag) $(xmldir)/m_debug.f90

m_dictionary.o:   $(xmldir)/m_dictionary.f90 m_buffer.o
	$(f90) $(flag) $(xmldir)/m_dictionary.f90

m_elstack.o:   $(xmldir)/m_elstack.f90 m_buffer.o
	$(f90) $(flag) $(xmldir)/m_elstack.f90

m_entities.o:   $(xmldir)/m_entities.f90 m_buffer.o
	$(f90) $(flag) $(xmldir)/m_entities.f90

m_fsm.o:   $(xmldir)/m_fsm.f90 m_buffer.o m_charset.o m_dictionary.o m_elstack.o m_entities.o
	$(f90) $(flag) $(xmldir)/m_fsm.f90

m_io.o:   $(xmldir)/m_io.f90
	$(f90) $(flag) $(xmldir)/m_io.f90

m_reader.o:   $(xmldir)/m_reader.f90 m_io.o
	$(f90) $(flag) $(xmldir)/m_reader.f90

m_xml_error.o:   $(xmldir)/m_xml_error.f90 m_elstack.o
	$(f90) $(flag) $(xmldir)/m_xml_error.f90

m_xml_parser.o:   $(xmldir)/m_xml_parser.f90 m_buffer.o m_debug.o m_dictionary.o m_elstack.o m_entities.o m_fsm.o m_reader.o m_xml_error.o
	$(f90) $(flag) $(xmldir)/m_xml_parser.f90

# TARGET clean: leave just source code.

clean:
	rm -rf *.o *.mod \#* *~

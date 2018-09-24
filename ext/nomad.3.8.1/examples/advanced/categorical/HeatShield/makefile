UNAME              := $(shell uname)

EXE                 = nomad_HS.exe
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -ansi -Wall -O3

LIB_DIR             = $(NOMAD_HOME)/lib
LIB_NOMAD           = libnomad.so


CXXFLAGS            =           

ifeq ($(UNAME), Darwin)
CXXFLAGS          += -headerpad_max_install_names
endif

ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif

LDLIBS              = -lm -lnomad
INCLUDE             = -I$(NOMAD_HOME)/src -I$(NOMAD_HOME)/ext/sgtelib/src -I.

COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = nomad_HS.o HS_Param.o HS_Material.o

ifndef NOMAD_HOME
define ECHO_NOMAD
	@echo Please set NOMAD_HOME environment variable!
	@false
endef
endif

$(EXE): $(OBJS)
	$(ECHO_NOMAD)
	@echo "   building the scalar version ..."
	@echo "   exe file : "$(EXE)
	@$(COMPILATOR) -o $(EXE) $(OBJS) $(LDLIBS) $(CXXFLAGS) -L$(LIB_DIR)
ifeq ($(UNAME), Darwin)
	@install_name_tool -change $(LIB_NOMAD) $(LIB_DIR)/$(LIB_NOMAD) $(EXE)
endif

nomad_HS.o: nomad_HS.cpp HS_Param.o
	$(ECHO_NOMAD)
	$(COMPILE) nomad_HS.cpp

HS_Param.o: HS_Param.cpp HS_Param.hpp HS_Material.o
	$(ECHO_NOMAD)
	$(COMPILE) HS_Param.cpp

HS_Material.o: HS_Material.cpp HS_Material.hpp
	$(ECHO_NOMAD)
	$(COMPILE) HS_Material.cpp

clean:
	@echo "   cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "   cleaning trash files"
	@rm -f core *~
	@echo "   cleaning obj files"
	@rm -f $(OBJS)
	@echo "   cleaning exe file"
	@rm -f $(EXE)

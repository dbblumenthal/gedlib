UNAME := $(shell uname)

EXE                    = nomad_for_dll.exe
COMPILATOR             = g++
COMPILATOR_OPTIONS     = -O3 -ansi

LIB_DIR                = $(NOMAD_HOME)/lib
LIB_NOMAD              = libnomad.so 

CXXFLAGS               =           
ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif

LDLIBS                 = -lm -lnomad

INCLUDE                = -I$(NOMAD_HOME)/src -I.

COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c

OBJS                   = nomad_for_dll.o

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
	@install_name_tool -change $(LIB_NOMAD) $(NOMAD_HOME)/lib/$(LIB_NOMAD) $(EXE)
endif


nomad_for_dll.o: nomad_for_dll.cpp
	$(ECHO_NOMAD)
	@$(COMPILE) nomad_for_dll.cpp

$(L1): ;
	$(ECHO_NOMAD)

clean:
	@echo cleaning obj files
	@del $(OBJS)

del:
	@echo cleaning obj files
	@del $(OBJS)
	@echo cleaning exe file
	@del $(EXE)

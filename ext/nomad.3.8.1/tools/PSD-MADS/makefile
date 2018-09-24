UNAME := $(shell uname)

EXE                    = psdmads.exe
COMPILATOR             = mpic++

COMPILATOR_OPTIONS     = -ansi -O2

LIB_DIR                = $(NOMAD_HOME)/lib
LIB_NOMAD              = libnomad.so 

CXXFLAGS               =           
ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif

ifeq ($(UNAME), Darwin)
CXXFLAGS              += -headerpad_max_install_names
endif

LDLIBS                 = -lm -lnomad -lmpi

INCLUDE                = -I$(NOMAD_HOME)/ext/sgtelib/src -I$(NOMAD_HOME)/src -I.
COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                   = main.o Cache_Server.o Master_Slaves.o

ifndef NOMAD_HOME
define ECHO_NOMAD
	@echo Please set NOMAD_HOME environment variable!
	@false
endef
endif

$(EXE): $(OBJS)
	$(ECHO_NOMAD)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LDLIBS) $(CXXFLAGS) -L$(LIB_DIR) 
ifeq ($(UNAME), Darwin)
	@install_name_tool -change $(LIB_NOMAD) $(NOMAD_HOME)/lib/$(LIB_NOMAD) $(EXE)
endif

main.o: main.cpp Master_Slaves.o
	$(ECHO_NOMAD)
	$(COMPILE) main.cpp

Master_Slaves.o: Master_Slaves.cpp Master_Slaves.hpp Cache_Server.o
	$(ECHO_NOMAD)
	$(COMPILE) Master_Slaves.cpp

Cache_Server.o: Cache_Server.cpp Cache_Server.hpp 
	$(ECHO_NOMAD)
	$(COMPILE) Cache_Server.cpp
	
clean:
	@echo "cleaning obj files"
	@rm -f $(OBJS)

del:
	@echo "cleaning trash files"
	@rm -f core *~
	@echo "cleaning obj files"
	@rm -f $(OBJS)
	@echo "cleaning exe file"
	@rm -f $(EXE)

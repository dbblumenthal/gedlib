UNAME := $(shell uname)

EXE                    = test.exe
COMPILATOR             = g++
FCOMPILATOR            = gfortran
COMPILATOR_OPTIONS     = -ansi -Wall -O3 

LIB_DIR                = $(NOMAD_HOME)/lib
LIB_NOMAD              = libnomad.so 

CXXFLAGS               =           
ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif


LDLIBS                 = -lm -lnomad -lstdc++

INCLUDE                = -I$(NOMAD_HOME)/src -I$(NOMAD_HOME)/ext/sgtelib/src -I.
COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
FCOMPILE               = $(FCOMPILATOR) -O3 -c

OBJS                   = nomad.o test.o

ifndef NOMAD_HOME
define ECHO_NOMAD
	@echo Please set NOMAD_HOME environment variable!
	@false
endef
endif

$(EXE): $(OBJS)
	$(ECHO_NOMAD)
	$(FCOMPILATOR) -o $(EXE) $(OBJS) $(LDLIBS) -L$(LIB_DIR) $(CXXFLAGS)

nomad.o: nomad.cpp
	$(ECHO_NOMAD)
	$(COMPILE) nomad.cpp

test.o: test.f
	$(ECHO_NOMAD)
	$(FCOMPILE) test.f
	
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

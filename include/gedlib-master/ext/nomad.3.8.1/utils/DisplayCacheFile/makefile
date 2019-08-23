UNAME := $(shell uname)

EXE                    = displayCache.exe

COMPILATOR             = g++
COMPILATOR_OPTIONS     = -O2 -ansi 

LIB_DIR                = $(NOMAD_HOME)/lib
LIB_NOMAD              = libnomad.so 

CXXFLAGS               =           
ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif

LDLIBS                 = -lm -lnomad

INCLUDE                = -I$(NOMAD_HOME)/src -I$(NOMAD_HOME)/ext/sgtelib/src -I.
COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c

OBJS                   = displayCache.o

ifndef NOMAD_HOME
define ECHO_NOMAD
	@echo Please set NOMAD_HOME environment variable!
	@false
endef
endif


$(EXE): $(OBJS)
	$(ECHO_NOMAD)
	@echo "   building the executable "$(EXE)
	@$(COMPILATOR) -o $(EXE) $(OBJS) $(LDLIBS) $(CXXFLAGS) -L$(LIB_DIR) 
ifeq ($(UNAME), Darwin)
	@install_name_tool -change $(LIB_NOMAD) $(NOMAD_HOME)/lib/$(LIB_NOMAD) $(EXE)
endif


displayCache.o: displayCache.cpp
	$(ECHO_NOMAD)
	@$(COMPILE) displayCache.cpp


all: $(EXE)

clean: ;
	@echo "   cleaning obj files"
	@rm -f $(OBJS)

del: ;
	@echo "   cleaning trash files"
	@rm -f core *~
	@echo "   cleaning obj files"
	@rm -f $(OBJS)
	@echo "   cleaning exe file"
	@rm -f $(EXE) 



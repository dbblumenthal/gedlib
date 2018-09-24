UNAME := $(shell uname)

EXE                    = cache_inspect.exe
COMPILATOR             = g++
COMPILATOR_OPTIONS     = -ansi -Wall -O2

LIB_DIR                = $(NOMAD_HOME)/lib
LIB_NOMAD              = libnomad.so 

CXXFLAGS               =           
ifeq ($(UNAME), Linux)
CXXFLAGS              += -Wl,-rpath,'$(LIB_DIR)'
CXXFLAGS              += -ansi
endif

LDLIBS                 = -lm -lnomad

INCLUDE                = -I$(NOMAD_HOME)/src -I. -I$(NOMAD_HOME)/ext/sgtelib/src
COMPILE                = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                   = cache_inspect.o

$(EXE): $(OBJS)
	$(ECHO_NOMAD)
	@$(COMPILATOR) -o $(EXE) $(OBJS) $(LDLIBS) $(CXXFLAGS) -L$(LIB_DIR) 
ifeq ($(UNAME), Darwin)
	@install_name_tool -change $(LIB_NOMAD) $(NOMAD_HOME)/lib/$(LIB_NOMAD) $(EXE)
endif


cache_inspect.o: cache_inspect.cpp
	$(ECHO_NOMAD)
	$(COMPILE) cache_inspect.cpp

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

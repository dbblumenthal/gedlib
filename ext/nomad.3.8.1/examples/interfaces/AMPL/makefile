LIBAMPL_DIR = $(LIBAMPLDIR)

CXX = gcc -std=c99
INCLUDE_DIRS = $(LIBAMPL_DIR)/Src/solvers
CXX_OPTS = $(addprefix -I,$(INCLUDE_DIRS))
LIBS = -L$(LIBAMPL_DIR)/Lib -lampl -lfuncadd0 -lm -ldl

%.o: %.c
	$(CXX) -c $(CXX_OPTS) $?

bb.exe: bb.o
	$(CXX) -o $@ $? $(LIBS)

clean:
	rm *.o miniampl

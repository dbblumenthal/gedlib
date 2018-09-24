EXE                 = surrogate.exe
COMPILATOR          = g++
COMPILATOR_OPTIONS  = -ansi -Wall -O3
LIBS                = -lc -lm
INCLUDE             = -I.
COMPILE             = $(COMPILATOR) $(COMPILATOR_OPTIONS) $(INCLUDE) -c
OBJS                = bb.o chemical.o thermolib.o secant.o stream.o servor.o mix.o bissection.o pump.o heatx.o reactor.o pfr.o reaction.o RungeKutta.o flash.o column.o split.o burner.o combrx.o cashflow.o profitability.o

$(EXE): $(OBJS)
	$(COMPILATOR) -o $(EXE) $(OBJS) $(LIBS) $(COMPILATOR_OPTIONS)

bb.o: bb.cpp servor.o profitability.o
	$(COMPILE) bb.cpp

profitability.o: profitability.cpp profitability.hpp cashflow.o secant.o
	$(COMPILE) profitability.cpp

cashflow.o: cashflow.cpp cashflow.hpp defines.hpp
	$(COMPILE) cashflow.cpp

servor.o: servor.cpp servor.hpp mix.o pump.o reactor.o column.o split.o burner.o
	$(COMPILE) servor.cpp

burner.o: burner.cpp burner.hpp combrx.o stream.o
	$(COMPILE) burner.cpp

combrx.o : combrx.cpp combrx.hpp chemical.o
	$(COMPILE) combrx.cpp

column.o: column.cpp column.hpp flash.o
	$(COMPILE) column.cpp

reactor.o: reactor.cpp reactor.hpp pfr.o
	$(COMPILE) reactor.cpp

pfr.o: pfr.cpp pfr.hpp RungeKutta.o stream.o
	$(COMPILE) pfr.cpp

reaction.o: reaction.cpp reaction.hpp chemical.o
	$(COMPILE) reaction.cpp

RungeKutta.o: RungeKutta.cpp RungeKutta.hpp defines.hpp
	$(COMPILE) RungeKutta.cpp

pump.o: pump.cpp pump.hpp stream.o
	$(COMPILE) pump.cpp

flash.o: flash.cpp flash.hpp stream.o bissection.o
	$(COMPILE) flash.cpp

split.o: split.cpp split.hpp stream.o
	$(COMPILE) split.cpp

mix.o: mix.cpp mix.hpp stream.o bissection.o
	$(COMPILE) mix.cpp

heatx.o: heatx.cpp heatx.hpp stream.o bissection.o
	$(COMPILE) heatx.cpp

bissection.o: bissection.cpp bissection.hpp defines.hpp
	$(COMPILE) bissection.cpp

stream.o: stream.cpp stream.hpp chemical.o
	$(COMPILE) stream.cpp

chemical.o: chemical.cpp chemical.hpp thermolib.o
	$(COMPILE) chemical.cpp

thermolib.o: thermolib.cpp thermolib.hpp secant.cpp secant.o
	$(COMPILE) thermolib.cpp

secant.o: secant.cpp secant.hpp
	$(COMPILE) secant.cpp

clean:
	@echo "cleaning obj files"
	@rm -f $(OBJS)
	@echo "cleaning trash files"
	@\rm -f core *~
del:
	@echo "cleaning trash files"
	@\rm -f core *~
	@echo "cleaning obj files"
	@\rm -f $(OBJS)
	@echo "cleaning exe file"
	@\rm -f $(EXE)

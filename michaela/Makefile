# vars
EXE     = calcE
CPP     = g++
LL      = g++
WARN    = -Wall
OPTI    = -O2
DEFS    = 
FLAGS   = -g
INCS    = -I/usr/include/openbabel-2.0
CCFLAGS = $(FLAGS) $(WARN) $(OPTI) $(DEFS) $(INCS)
LDFLAGS = -L/usr/lib
LIBS    = -lopenbabel
CPPSRC  = energy_calculation.cpp
CPPOBJ  = $(foreach file, $(CPPSRC:.cpp=.o), $(file))
XTRA    = Makefile
# rules
all: $(EXE)
$(CPPOBJ):  %.o: %.cpp
	$(CPP) -c $(CCFLAGS) $< -o $@
$(EXE): $(CPPOBJ)
	$(LL) -o $(EXE) $(CPPOBJ) $(LDFLAGS) $(LIBS)
clean:
	rm -f *~ $(CPPOBJ) $(EXE) 
PHONY: clean

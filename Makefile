include Makefile.common
TARGET		= bin/test_differential_equations.bin
OBJS		= bin/differential_equations.o bin/io_functions.o bin/Timer.o bin/test_differential_equations.o bin/Trace.o
CPPFLAGS	= $(DEFAULT_CPPFLAGS)
LDFLAGS		= $(DEFAULT_LDFLAGS)
PROJECTDIR = .

all:
	@mkdir -p bin
	@echo make $(TARGET)
	@$(MAKE) --silent $(OBJS)
	@$(MAKE) --silent $(TARGET)
	@echo Done.

# object files:
bin/differential_equations.o: include/automatic_derivatives/Fad_construct.h include/Trace.h include/automatic_derivatives/Fad.h include/io_functions.h src/differential_equations.cpp include/differential_equations.h include/constants.h include/automatic_derivatives/Fad_operators.h include/automatic_derivatives/Fad_functions.h
	@echo "   $(CCC) bin/differential_equations.o"
	@$(CCC) -o bin/differential_equations.o -c src/differential_equations.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion
bin/io_functions.o: include/io_functions.h include/constants.h src/io_functions.cpp
	@echo "   $(CCC) bin/io_functions.o"
	@$(CCC) -o bin/io_functions.o -c src/io_functions.cpp $(CPPFLAGS) 
bin/Timer.o: include/io_functions.h src/Timer.cpp include/Timer.h
	@echo "   $(CCC) bin/Timer.o"
	@$(CCC) -o bin/Timer.o -c src/Timer.cpp $(CPPFLAGS) 
bin/test_differential_equations.o: src/test_differential_equations.cpp include/automatic_derivatives/Fad_construct.h include/Trace.h include/automatic_derivatives/Fad.h include/Timer.h include/io_functions.h include/differential_equations.h include/constants.h include/automatic_derivatives/Fad_operators.h include/automatic_derivatives/Fad_functions.h
	@echo "   $(CCC) bin/test_differential_equations.o"
	@$(CCC) -o bin/test_differential_equations.o -c src/test_differential_equations.cpp $(CPPFLAGS) -Wno-shadow -I${PROJECTDIR}/include -Wno-conversion
bin/Trace.o: include/io_functions.h include/Trace.h src/Trace.cpp
	@echo "   $(CCC) bin/Trace.o"
	@$(CCC) -o bin/Trace.o -c src/Trace.cpp $(CPPFLAGS) 

# target files (executable or shared object):
bin/test_differential_equations.bin: bin/differential_equations.o bin/io_functions.o bin/Timer.o bin/test_differential_equations.o bin/Trace.o
	@echo "   $(CCC) bin/test_differential_equations.bin"
	@$(CCC) -o bin/test_differential_equations.bin bin/differential_equations.o bin/io_functions.o bin/Timer.o bin/test_differential_equations.o bin/Trace.o $(LDFLAGS) 

clean:
	@rm -f $(OBJS)
	@rm -f $(TARGET)
	@if [ -d bin ]; then rm -r bin ; fi

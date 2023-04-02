C = gcc
CCFLAGS = -Wall -Wextra -Wpedantic -O1
GFLAGS = -Wall -Wextra -Wpedantic -Og -g
LFLAGS = -lgsl -lgslcblas -lm
SDLFLAGS = -lSDL2

default: output

output: orbit.c
	$(C) orbit.c $(GFLAGS) $(LFLAGS) -o output.o
	
output2: hamiltonian.c
	$(C) hamiltonian.c $(GFLAGS) $(LFLAGS) -o output2.o

output3: poincare_section.c
	$(C) poincare_section.c $(GFLAGS) $(LFLAGS) -o output3.o
C = gcc
CCFLAGS = -Wall -Wextra -Wpedantic -O1
GFLAGS = -Wall -Wextra -Wpedantic -Og -g
LFLAGS = -lgsl -lgslcblas -lm

default: output
	
output: hamiltonian.c
	$(C) $^ $(GFLAGS) $(LFLAGS) -o $@.o

output2: poincare_section.c
	$(C) $^ $(GFLAGS) $(LFLAGS) -o $@.o
	
ouput3: sec_poincare.c
	$(C) $^ $(CCFLAGS) $(LFLAGS) -o $@.o
	
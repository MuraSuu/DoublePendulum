CC = g++
C = gcc
CCFLAGS = -Wall -Wextra -Wpedantic -O1
GFLAGS = -Wall -Wextra -Wpedantic -Og -g
LFLAGS = -lgsl -lgslcblas -lm
SDLFLAGS = -lSDL2

default: output

output: lagrangian.cc
	$(CC) lagrangian.cc $(GFLAGS) $(LFLAGS) -o output.o
	
output2: hamiltonian.c
	$(C) hamiltonian.c $(GFLAGS) $(LFLAGS) -o output2.o
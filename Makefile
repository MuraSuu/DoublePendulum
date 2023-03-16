CC = g++
CCFLAGS = -Wall -Wextra -Wpedantic -O1
GFLAGS = -Wall -Wextra -Wpedantic -Og -g
LFLAGS = -lgsl -lgslcblas -lm
SDLFLAGS = -lSDL2

default: output

output: main.cc
	$(CC) main.cc $(GFLAGS) $(LFLAGS) -o output.o
	

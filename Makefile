CC = g++
CFLAGS = -c -IBasisClasses/ -Wall
LFLAGS = -lm -IBasisClasses/ -Wall

OBJECTS = main.o gamma.o coupling.o basis.o

all: main

main: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main

main.o: main.cpp gamma.h
	$(CC) $(CFLAGS) main.cpp 

gamma.o: gamma.h gamma.cpp
	$(CC) $(CFLAGS) gamma.cpp 

coupling.o: coupling.h coupling.cpp gamma.h
	$(CC) $(CFLAGS) coupling.cpp

basis.o: basis.h HOstate.h 
	$(CC) $(CFLAGS) basis.cpp

clean: 
	rm *.o BasisClasses/*.o main

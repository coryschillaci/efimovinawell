CC = g++
CFLAGS = -c -I../eigen/ -Wall
LFLAGS = -lm -Wall

OBJECTS = main.o gamma.o coupling.o basis.o HOstate.o

all: main

main: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main

main.o: main.cpp gamma.h coupling.h basis.h
	$(CC) $(CFLAGS) main.cpp 

gamma.o: gamma.h gamma.cpp
	$(CC) $(CFLAGS) gamma.cpp 

coupling.o: coupling.h coupling.cpp gamma.h
	$(CC) $(CFLAGS) coupling.cpp

basis.o: basis.h basis.cpp HOstate.h 
	$(CC) $(CFLAGS) basis.cpp

HOstate.o: HOstate.h HOstate.cpp
	$(CC) $(CFLAGS) HOstate.cpp

clean: 
	rm *.o main

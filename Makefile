CC = g++
CFLAGS = -c -Wall
LFLAGS = -lm -Wall

OBJECTS = main.o gamma.o coupling.o

all: main

main: $(OBJECTS)
	$(CC) $(LFLAGS) $(OBJECTS) -o main

main.o: main.cpp gamma.h
	$(CC) $(CFLAGS) main.cpp 

gamma.o: gamma.h gamma.cpp
	$(CC) $(CFLAGS) gamma.cpp 

coupling.o: coupling.h coupling.cpp gamma.h
	$(CC) $(CFLAGS) coupling.cpp

clean: 
	rm *.o main

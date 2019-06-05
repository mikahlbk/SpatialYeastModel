CC=g++ -std=c++11

CFLAGS=-c -Wall

all: program

program: main.o coord.o cell.o colony.o
		$(CC) main.o coord.o cell.o colony.o -o program

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

colony.o: colony.cpp
		$(CC) $(CFLAGS) colony.cpp

clean: wipe
		rm -rf *o program

wipe:

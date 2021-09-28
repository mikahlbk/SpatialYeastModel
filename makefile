CC=-g++ -fopenmp -static-libstdc++ -std=c++11

CFLAGS=-c -Wall -O3

all: program

program: main.o coord.o colony.o mesh.o mesh_pt.o cell.o
		$(CC) main.o coord.o colony.o mesh.o mesh_pt.o cell.o -o program

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

colony.o: colony.cpp
		$(CC) $(CFLAGS) colony.cpp

mesh.o: mesh.cpp
		$(CC) $(CFLAGS) mesh.cpp

mesh_pt.o: mesh_pt.cpp
		$(CC) $(CFLAGS) mesh_pt.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

clean: wipe
		rm -rf *o program
wipe:

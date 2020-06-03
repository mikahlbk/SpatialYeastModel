CC=-g++ -fopenmp -static-libstdc++

CFLAGS=-c -Wall -O3

all: program

program: main.o coord.o cell.o colony.o mesh_pt.o mesh.o
		$(CC) main.o coord.o cell.o colony.o mesh_pt.o mesh.o -o program

main.o: main.cpp
		$(CC) $(CFLAGS) main.cpp

coord.o: coord.cpp
		$(CC) $(CFLAGS) coord.cpp

cell.o: cell.cpp
		$(CC) $(CFLAGS) cell.cpp

colony.o: colony.cpp
		$(CC) $(CFLAGS) colony.cpp

mesh_pt.o: mesh_pt.cpp
		$(CC) $(CFLAGS) mesh_pt.cpp

mesh.o: mesh.cpp
		$(CC) $(CFLAGS) mesh.cpp

clean: wipe
		rm -rf *o program

wipe:

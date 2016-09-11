CC=g++
OUT=licz

program: RozwRown.h RozwRown.o main.cpp
	$(CC) RozwRown.o main.cpp -o $(OUT)

RozwRown.o: RozwRown.h RozwRown.cpp
	$(CC) -c RozwRown.cpp

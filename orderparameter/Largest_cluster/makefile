
CC=g++
CX=mpicxx


lib:
	$(CC) -O3 -Wall -fPIC -c steinhardt.cpp
	$(CC) -O3 -Wall -shared -o steinhardt.so steinhardt.o

all:
	#$(CC) -O3 -Wall -fPIC -c steinhardt.cpp main.cpp -o struct
	$(CC) -g -lm -Wall -o struct steinhardt.cpp main.cpp 

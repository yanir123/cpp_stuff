CC = g++
CFLAGS = -lm -fopenmp -O3
TARGET = main.cpp

dafault: main

main: ann.o $(TARGET)
	$(CC) $(CFLAGS) $(TARGET) ann.o -o main.out

ann.o: ann.cpp ann.hpp
	$(CC) $(CFLAGS) -c ann.cpp

clean: 
	$(RM) main.out *.o
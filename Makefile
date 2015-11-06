CC = g++
FLAGS = -ansi -std=c++11 -lm 
MAIN_SRC = matrix.cpp

matrix: matrix.cpp
	g++ matrix.cpp -o matrix $(FLAGS)
clean:
	rm -r matrix
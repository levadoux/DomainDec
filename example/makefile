GCC=g++-13 -std=c++23 -Wall
FLAGS=-O3
INC=-I../eigen/ -I../femtool/

example: example.o
	$(GCC) $(FLAGS) example.o -o example.exe


example.o: example.cpp
	$(GCC) $(FLAGS) $(INC) -c example.cpp -o example.o

vizirhere:  
	../../vizir4.2025.01.22.windows/vizir4.exe -in

clean:
	rm *.o *.exe output*



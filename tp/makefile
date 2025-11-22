GCC=g++-13 -std=c++23 -Wall
FLAGS=-O3
INC=-I../eigen/ -I../femtool/

tp1-exo1: tp1-exo1.o
	$(GCC) $(FLAGS) tp1-exo1.o -o tp1-exo1.exe

tp1-exo1.o: tp1-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo1.cpp -o tp1-exo1.o

tp1-exo2: tp1-exo2.o
	$(GCC) $(FLAGS) tp1-exo2.o -o tp1-exo2.exe

tp1-exo2.o: tp1-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo2.cpp -o tp1-exo2.o

tp1-exo3: tp1-exo3.o
	$(GCC) $(FLAGS) tp1-exo3.o -o tp1-exo3.exe

tp1-exo3.o: tp1-exo3.cpp
	$(GCC) $(FLAGS) $(INC) -c tp1-exo3.cpp -o tp1-exo3.o

tp2-exo1: tp2-exo1.o
	$(GCC) $(FLAGS) tp2-exo1.o -o tp2-exo1.exe

tp2-exo1.o: tp2-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo1.cpp -o tp2-exo1.o

tp2-exo2: tp2-exo2.o
	$(GCC) $(FLAGS) tp2-exo2.o -o tp2-exo2.exe

tp2-exo2.o: tp2-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp2-exo2.cpp -o tp2-exo2.o

tp3-exo1: tp3-exo1.o
	$(GCC) $(FLAGS) tp3-exo1.o -o tp3-exo1.exe

tp3-exo1.o: tp3-exo1.cpp
	$(GCC) $(FLAGS) $(INC) -c tp3-exo1.cpp -o tp3-exo1.o

tp3-exo2: tp3-exo2.o
	$(GCC) $(FLAGS) tp3-exo2.o -o tp3-exo2.exe

tp3-exo2.o: tp3-exo2.cpp
	$(GCC) $(FLAGS) $(INC) -c tp3-exo2.cpp -o tp3-exo2.o

tp3-exo3: tp3-exo3.o
	$(GCC) $(FLAGS) tp3-exo3.o -o tp3-exo3.exe

tp3-exo3.o: tp3-exo3.cpp
	$(GCC) $(FLAGS) $(INC) -c tp3-exo3.cpp -o tp3-exo3.o

H_LIST = 0.025 0.01 0.005 0.0025

iter:
	
	@for h in $(H_LIST); do \
		echo "Lancement avec h = $$h"; \
		nb_it=$$(./tp1-exo3.exe $$h | tail -n 1); \
		echo "$$h ; $$nb_it" >> dat/iter.dat; \
	done

iter_precond:
	
	@for h in $(H_LIST); do \
		echo "Lancement avec h = $$h"; \
		nb_it=$$(./tp2-exo1.exe $$h | tail -n 1); \
		echo "$$h ; $$nb_it" >> dat/iter_precond.dat; \
	done

iter_precond_Schwarz_add:
	
	@for h in $(H_LIST); do \
		echo "Lancement avec h = $$h"; \
		nb_it=$$(./tp3-exo3.exe $$h | tail -n 1); \
		echo "$$h ; $$nb_it" >> dat/iter_precond_Schwarz_add.dat; \
	done

vizirhere:
	../../vizir4.2025.01.22.windows/vizir4.exe &

plot_conv:
	python3 plot_conv.py

plot_conv_precond:
	python3 plot_conv_precond.py

plot_iter:
	python3 plot_iter.py

clean:
	rm *.o *.exe

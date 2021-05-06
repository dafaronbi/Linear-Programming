#build executable of two variable
all: three\ variable.c three\ variable.h pairing.c main.c random_problem.c dyer_algo.c median.c
		cc -g -Wall -O -c three\ variable.c pairing.c main.c random_problem.c dyer_algo.c median.c
		cc  three\ variable.o pairing.o main.o random_problem.o dyer_algo.o median.o -o exec
		
clean:
	$(RM) exec
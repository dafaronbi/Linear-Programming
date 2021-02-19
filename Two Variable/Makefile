#build executable of two variable
all: main.c pairing.c median.c random_problem.c  dyer_algo.c two\ variable.h
		cc -g -Wall -O -c main.c pairing.c median.c random_problem.c dyer_algo.c 
		cc  main.o pairing.o median.o random_problem.o  dyer_algo.o -o tv 
		
clean:
	$(RM) myprog
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "defines.h"
#include <fstream>
 
using namespace std;

#define ROWS 10
#define COLUMNS 10   
#define STARTING_SEED 10000
#define INCREMENTAL_SEED 2000
#define NUMBER_PROBLEMS 10
#define A_MULTIPLE_FACTOR 1000
#define A_SUM_FACTOR 0
#define B_TIGHTNESS_RATIO 0.5
#define C_MULTIPLE_FACTOR 200
#define NEG_NUMBERS 0   // 0 = No	| 1 = Yes
#define PROB_NEGATIVES 0.1
#define SPARSE 0        // 0 = No Sparse Matrix		| 1 = Sparse Matrix
#define SPARSE_PERCENTAGE 0.90

struct Rows_Double {
    double *cols;
};
 
struct Rows_Integer {
    int *cols;
};
 
struct Instance {
    double *cost_vector;
    double *rhs_vector;
    Rows_Double *a_matrix;
    int number_variables;
    int number_constraints;
};
 
typedef struct Instance INSTANCE, *INSTANCEptr;
 
double function_random(double random_number);
void create_instance(Instance *instance);
void create_sparse(Instance *instance);

static double z;
long long int seed_value = 0;
 
int main(void) {
    int i, j, k, inc_seed;
	char fn[256];

    FILE *fp = NULL;
    Instance *instance;
	instance = (Instance*)calloc(1, sizeof(Instance));
	instance->cost_vector = (double*)calloc(COLUMNS, sizeof(double));
    instance->rhs_vector = (double*)calloc(ROWS, sizeof(double));
    instance->a_matrix = (Rows_Double*)calloc(ROWS, sizeof(Rows_Double));
	for (i = 0; i < ROWS; i++) {
		instance->a_matrix[i].cols = (double*)calloc(COLUMNS, sizeof(double));
    }
    instance->number_variables = COLUMNS;
    instance->number_constraints = ROWS;

	inc_seed = 0;
	for (i = 0; i < NUMBER_PROBLEMS; i++) {
		seed_value = STARTING_SEED + inc_seed;
        inc_seed += INCREMENTAL_SEED;
		sprintf(fn, "Instances/%d_%d_%d.txt", instance->number_constraints, instance->number_variables, seed_value);
		if (!(fp = fopen(fn, "w"))) {
			fprintf(stderr, "There was a problem opening the output file!\n");
			goto TERMINATE;
		}
		if (SPARSE == 1) {
			create_sparse(instance);
		}
		else {
			create_instance(instance);
		}
		for (j = 0; j < instance->number_variables; j++) {
			fprintf(fp, "%f ", instance->cost_vector[j]);
		}
		fprintf(fp, "\n\n");		
		for (j = 0; j < instance->number_constraints; j++) {
			for (k = 0; k < instance->number_variables; k++) {
				fprintf(fp, "%f ", instance->a_matrix[j].cols[k]);
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n");
		for (j = 0; j < instance->number_constraints; j++) {
			fprintf(fp, "%f\n", instance->rhs_vector[j]);
		}
		fclose(fp);
    }
	free(instance);
 
TERMINATE:

    system("PAUSE");
    return 0;
}
 
double function_random(double random_number) {
    const int a = 100801;
    const int c = 103319;
    const int m = 193723;
    seed_value = (a*seed_value*c) % m;
    random_number = (double)(seed_value) / (double)(m);
    return random_number;
}
 
void create_instance(Instance *instance) {
    int i, j;
    double random_number = 0, auxiliar, sum, store;
    for (i = 0; i < instance->number_constraints; i++) {
        sum = 0;
        for (j = 0; j < instance->number_variables; j++) {
            auxiliar = function_random(random_number);
            store = ceil(auxiliar * A_MULTIPLE_FACTOR + A_SUM_FACTOR);
            if (store == 0) {
                store = 1;
            }
            if (NEG_NUMBERS == 1) {
                auxiliar = function_random(random_number);
                if (auxiliar <= PROB_NEGATIVES) {
                    store = -store;
                }
            }
            instance->a_matrix[i].cols[j] = store;
            sum += instance->a_matrix[i].cols[j];
        }
		instance->rhs_vector[i] = ceil(sum * B_TIGHTNESS_RATIO);
    }
    for (j = 0; j < instance->number_variables; j++) {
        sum = 0;
        for (i = 0; i < instance->number_constraints; i++) {
            sum += instance->a_matrix[i].cols[j];
        }
        auxiliar = function_random(random_number);
        store = ceil(auxiliar * C_MULTIPLE_FACTOR);
        instance->cost_vector[j] = sum + store;
    }
}
 
void create_sparse(Instance *instance) {
    int i, j;
    double random_number = 0, auxiliar, sum, store;
 
    for (i = 0; i < instance->number_constraints; i++) {
        for (j = 0; j < (instance->number_variables + instance->number_constraints); j++) {
            instance->a_matrix[i].cols[j] = 0;
            instance->cost_vector[i] = 0;
        }
        instance->rhs_vector[i] = 0;
    }
    for (i = 0; i < instance->number_constraints; i++) {
        sum = 0;
        for (j = 0; j < instance->number_variables; j++) {
            auxiliar = function_random(random_number);
            if (auxiliar > SPARSE_PERCENTAGE) {
                if (NEG_NUMBERS == 1) {
                    auxiliar = function_random(random_number);
                    if (auxiliar <= PROB_NEGATIVES) {
                        instance->a_matrix[i].cols[j] = -1;
                    }
                    else {
                        instance->a_matrix[i].cols[j] = 1;
                    }
                }
                else {
                    instance->a_matrix[i].cols[j] = 1;
                }
            }
            sum += instance->a_matrix[i].cols[j];
        }
        instance->rhs_vector[i] = ceil(sum * B_TIGHTNESS_RATIO);
    }
    for (j = 0; j < instance->number_variables; j++) {
        sum = 0;
        for (i = 0; i < instance->number_constraints; i++) {
            sum += instance->a_matrix[i].cols[j];
        }
        auxiliar = function_random(random_number);
        store = ceil(auxiliar * C_MULTIPLE_FACTOR);
        instance->cost_vector[j] = sum + store;
    }
}